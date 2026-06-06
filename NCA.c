
#include "udf.h"
#include <stdbool.h>
#include "hdfio.h"
#include <math.h>

static real V_f = 0.0; // Initialize flame spread rate variable, will be updated at end of each iteration in calc_FSR and used in inlet velocity profile and solid motion BCs
static real alpha = 1; // Under-relaxation factor for FSR update, adjust as needed for stability and convergence speed

// Temp based flame spread rate calculation variables
static real R;
static real V_f_old;
static real R_old;
static real dR_by_dVf;
static int eigen_face_zoneID = 0; // Zone ID of separated surface where temp is monitored
static real Delta_Vf = 0.000001; // Initial change in FSR
static const real T_infty = 300; // Ambient temperature [K]
static bool converged_flag = false;

// Eigenposition based FSR calculation, Bisection method
static real V_f_bracket[2] = {0.0, 0.00015}; //Initial Bracet for FSR 0 - 200 um/s
static real R_bracket[2] = {338.0, -60}; //Bracket for residual, lower bracket will be updated after first iteration. Uppwer 
static int N_UPDATE = 1000; // Number of iterations to wait between FSR updates in bisection method, adjust as needed for stability and convergence speed
static int N_FSR_ITER = 0; // Iteration counter at last FSR change, used to control update frequency in eigenvalue method

// Eigenposition PID controller method
static real e; // Error between eigen position temperature and target temperature
// Function to calculate if scaled residuslas are below a threshold
static real false_position_residual_tolerance = 1e-5; // Residual tolerance for accepting a given R(V_f) in false position method, adjust as needed
static int last_fsr_update_iter = -100000000;
static int min_iter_after_fsr_update = 1000;  

static bool residuals_below(Domain *d, real threshold)
{
    real scaled_res;
    real max_scaled_res = 0.0;
    int residuals_ok = 0;

    /* Loop over all available equations/species in Fluent */
#if !RP_HOST
    for (int nw = 0; nw < DOMAIN_NUMEQN(d); nw++)
    {
		int eqn = DOMAIN_EQNS(d, nw); // Get the equation index for the current equation/species
		const char* strEqnLabel = DOMAIN_EQN_LABEL(d, eqn);

		if ( strEqnLabel && (strlen(strEqnLabel) > 0) )
		{
			/* Check to prevent division by zero for unscaled residuals */
			if (0 == DOMAIN_RES_SCALE(d, nw)) DOMAIN_RES_SCALE(d,nw) = 1.0;
			
			scaled_res = DOMAIN_RES(d, nw) / DOMAIN_RES_SCALE(d,nw);

            #if RP_HOST
            #else
                max_scaled_res = MAX(max_scaled_res, fabs(scaled_res));
            #endif
		}
    }

	PRF_GRHIGH1(max_scaled_res); // Get the maximum scaled residual across all compute nodes
#endif

	node_to_host_real_1(max_scaled_res); // Update max_scaled_res on host process so it can be used in report definition and for convergence checking

    /* Evaluate against threshold */
	

#if RP_HOST
		residuals_ok = (max_scaled_res < threshold);
#endif 

	host_to_node_int_1(residuals_ok); // Update residuals_ok on compute nodes so it can be used in UDF logic
	
	return residuals_ok ? true : false;
}



/* typedef enum {
	FSR_RAMP_FROM_BELOW = 0;
	FSR_FALSE_POSITION = 1;
	FSR_DONE = 2;
} fsr_state_t;

fsr_state_t fsr_state = FSR_RAMP_FROM_BELOW; // Initial state for bisection method */


/* CODE SECTION */
/* FD INLET VELOCITY PROFILE */
/*
	Modified version of pNCA.c to use RP variable for mean velocity.
	Put your velocity in the brackets (remove the brackets)
	Define the variable using the following commands in Fluent TUI:
	To define: (rp-var-define `user/U_mean [velocity] 'real #f)
		   ex: (rp-var-define `user/U_mean 0.146073469 'real #f)
	To get   : (%rpgetvar `user/U_mean)
	To update: (rpsetvar `user/U_mean 0.1)
	To access in UDF: RP_Get_Real("user/U_mean");
				  ex: U_mean = RP_Get_Real("user/U_mean");
*/
DEFINE_PROFILE(inlet_x_vel_rpvar, thread, position)
{
	real x[ND_ND]; /* this will hold the position vector */
	real y, h, U_mean, U_max, m, n;
	face_t f;

	h = 0.00495; /* m; inlet height, do not change */
	m = 27.59596236; /* constant, do not change */
	n = 2.0; /* constant, do not change */

	// Get U_mean from RP var
	U_mean = 0.082; /* m/sec; inlet mean velocity, update with geom, default value */
	bool U_mean_exists = RP_Variable_Exists_P("user/u_mean"); // Check if user-defined parameter exists
	Message0("Checking for user-defined parameter 'user/u_mean': %d\n", U_mean_exists);

	if (RP_Variable_Exists_P("user/u_mean"))
	{
		U_mean = RP_Get_Real("user/u_mean"); // Get mean velocity from user-defined parameter if it exists
	}
	else
	{
		Message0("Warning: User-defined parameter 'user/u_mean' not found. Using default value of %f m/s.\n", U_mean);
	}

	U_mean += V_f; // Add calculated FSR to mean velocity for inlet profile

	U_max = U_mean * ((m + 1) / m) * ((n + 1) / n); /* m/sec; max velocity, at centerline... calc */

	begin_f_loop(f, thread)
	{
		F_CENTROID(x, f, thread);
		y = 2. * (x[1] - 0.5 * h) / h; /* non-dimensional y coordinate, b/c coord sys is at bottom of geom not centerline... calc */

		F_PROFILE(f, thread, position) = U_max * (1.0 - (y * y)); /* m/sec; velocity as f(y) at centerline... calc */
	}
	end_f_loop(f, thread)
}

DEFINE_PROFILE(inlet_x_vel_8cms_fsr_preset, thread, position)
{
	real x[ND_ND]; /* this will hold the position vector */
	real y, h, U_mean, U_max, m, n;
	face_t f;

	h = 0.00495; /* m; inlet height, do not change */
	m = 27.59596236; /* constant, do not change */
	n = 2.0; /* constant, do not change */

	U_mean = 0.082 + 0.0000617; /* m/sec; inlet mean velocity, update with geom */
	U_max = U_mean * ((m + 1) / m) * ((n + 1) / n); /* m/sec; max velocity, at centerline... calc */

	begin_f_loop(f, thread)
	{
		F_CENTROID(x, f, thread);
		y = 2. * (x[1] - 0.5 * h) / h; /* non-dimensional y coordinate, b/c coord sys is at bottom of geom not centerline... calc */

		F_PROFILE(f, thread, position) = U_max * (1.0 - y * y); /* m/sec; velocity as f(y) at centerline... calc */
	}
	end_f_loop(f, thread)
}

// Calculate Flame Spread Rate
DEFINE_EXECUTE_AT_END(calc_FSR)
{

	// Integral of the x-component of the normal vector f'(x)/sqrt(f'(x)^2+1
	// Units of I are m^2. Even though the integral is over one spatial direction, x
	// it obtains units of m^2 from multiplication by the reference length (depth) of 1 m

	const real I = 0.002243881756148; // [m^2] Computed in Matlab, changes with different surface profiles. Corresponds to Hossain V_g = 8.2 cm/s curve
	const real rho = 1190; // [kg/m^3] Density of solid phase



	real mdot_chem = 0.; //Mass flux from chemical reaction at surface [kg/s]
	//real V_f;

	face_t f; // Face along surface

#if !RP_HOST
	// Find wall_mass_flux thread
	Domain* d = Get_Domain(1); // Get domain pointer, update if different
	int zone_ID = 5; // ID of surface zone where chemical reaction occurs, update if different. Zone is shown in Boundary conditions tab
	Thread* t = Lookup_Thread(d, zone_ID); // Get thread pointer for surface zone where chemical reaction occurs

	// Loop through faces along surface and sum mass flux from chemical reaction
	begin_f_loop(f, t)
		if PRINCIPAL_FACE_P(f, t)
		{
			// Evem though this macro is labeled as "FLUX", it is actually a mass flow rate through a face 
			// according to section 3.2.2.4 of Fluent Customization manual. So we can sum this value across
			// all faces along the surface to get total mass flow rate from chemical reaction at surface.

			mdot_chem += F_FLUX(f, t); // Sum mass fluxes on each face from chemical reaction at surface
		}
	end_f_loop(f, t)

		// Sum mdot over all compute nodes
		mdot_chem = PRF_GRSUM1(mdot_chem);
	// Calculate corrected FSR and print result
	V_f = fabs(mdot_chem / (rho * I)); // Calculate flame spread rate [m/s]
#endif
	node_to_host_real_1(V_f); // update V_f on host process so report is correct.

	// Print Every 25 iterations to avoid excessive printing, update with different frequency if desired.
	// Count figure out how to automatically pass the profile update interval (count find a macro or rp var for it)
	if (N_ITER % 25 == 0)
	{
		Message0("Mass flux from chemical reaction at surface: %g kg/s\n", mdot_chem);
		Message0("Calculated Flame Spread Rate: %g m/s\n", V_f);
	}

}

// Calculate FSR using the least squares method
/*
V_f = \int_{A_{face}} \dot{m}''  \hat{n_x} dA_{face} / \rho \int_{A_{face}} \hat{n_x}^2 dA_{face}

*/
DEFINE_EXECUTE_AT_END(update_FSR_LSQ)
{
	const real rho = 1190; // [kg/m^3] Density of solid phase

	real mdot = 0.; //Mass flux at surface [kg/s]
	real mdot_face = 0; // Mass flux at a face [kg/s]

	//real I = 0; // Integral of x-component of face area vector, sum of A_face_x = sum of A_face * f'(x)/sqrt(f'(x)^2+1) over all faces along surface
	real A_face[ND_ND]; // Face area vector
	real A_face_mag; // Face area magnitude
	real nhat_x; //x-component of face normal vector
	//real V_f;

	// Integrals in the numerator and deniminator
	real I_num = 0;
	real I_denom = 0;

	face_t f; // Face along surface

#if !RP_HOST
	// Find wall_mass_flux thread
	Domain* d = Get_Domain(1); // Get domain pointer, update if different
	int zone_ID = 5; // ID of surface zone where chemical reaction occurs, update if different. Zone is shown in Boundary conditions tab
	int zone_fixed_ID = 65; // ID of surface of zone with fixed 668 temperature to pick eigenvalue

	Thread* t = Lookup_Thread(d, zone_ID); // Get thread pointer for surface zone where chemical reaction occurs
	Thread* t_fixed = Lookup_Thread(d, zone_fixed_ID); //pointer to fixed temp surface 

	// Loop through faces along surface and sum mass flux from chemical reaction
	begin_f_loop(f, t)
		if PRINCIPAL_FACE_P(f, t)
		{
			// Evem though this macro is labeled as "FLUX", it is actually a mass flow rate through a face 
			// according to section 3.2.2.4 of Fluent Customization manual. So we can sum this value across
			// all faces along the surface to get total mass flow rate from chemical reaction at surface.
			mdot_face = F_FLUX(f, t); // Mass flow from single face
			mdot += mdot_face; // Sum mass fluxes on each face from chemical reaction at surface

			F_AREA(A_face, f, t); // Get face area vector for current face
			A_face_mag = NV_MAG(A_face); // Get face area magnitude for current face

			//mass_flux = mdot_face / A_face_mag; //average mass flux at face [kg/m^2-s]

			nhat_x = A_face[0] / A_face_mag; // Get x-component of face normal vector for current face

			I_num += mdot_face * nhat_x; // Increment numerator integral with contribution from current face, mdot'' * n_x * A_face = mdot * n_x
			I_denom += nhat_x * nhat_x * A_face_mag; // Increment denominator integral with contribution from current face, n_x^2 * A_face

		}
	end_f_loop(f, t)

		// Add flux other reaction thread surface
		begin_f_loop(f, t_fixed)
		if PRINCIPAL_FACE_P(f, t_fixed)
		{
			mdot_face = F_FLUX(f, t_fixed); // Mass flow from single face
			mdot += mdot_face; // Sum mass fluxes on each face from chemical reaction at surface

			F_AREA(A_face, f, t_fixed); // Get face area vector for current face
			A_face_mag = NV_MAG(A_face); // Get face area magnitude for current face

			//mass_flux = mdot_face / A_face_mag; //average mass flux at face [kg/m^2-s]

			nhat_x = A_face[0] / A_face_mag; // Get x-component of face normal vector for current face

			I_num += mdot_face * nhat_x; // Increment numerator integral with contribution from current face, mdot'' * n_x * A_face = mdot * n_x
			I_denom += nhat_x * nhat_x * A_face_mag; // Increment denominator integral with contribution from current face, n_x^2 * A_face
		}
	end_f_loop(f, t_fixed)

		// Sum mdot over all compute nodes
		mdot = PRF_GRSUM1(mdot);

	// Sum integrals over all compute nodes
	I_num = PRF_GRSUM1(I_num);
	I_denom = PRF_GRSUM1(I_denom);

	// Calculate corrected FSR and print result
	//V_f = mdot_chem / (rho * I); // Calculate flame spread rate [m/s]
	double V_f_new = I_num / (rho * I_denom); // Calculate flame spread rate using least squares method [m/s]

	// Apply under-relaxation
	V_f = V_f + alpha * (V_f_new - V_f); // Update flame spread rate with under-relaxation factor of 0.5, adjust factor as needed for stability and convergence speed
#endif

	node_to_host_real_1(V_f); // update V_f on host process so report is correct.

	// Print Every 25 iterations to avoid excessive printing, update with different frequency if desired.
	// Count figure out how to automatically pass the profile update interval (count find a macro or rp var for it)
	if (N_ITER % 25 == 0)
	{
		Message0("Mass flow at surface: %g kg/s\n", mdot);
		Message0("LSQ Calculated Flame Spread Rate: %g m/s\n", V_f);
	}

}

/* // Flame spread rate calculation based on eigenposition temperature
DEFINE_EXECUTE_AT_END(calc_FSR_eigen)
{
	real T_eig = 0.0; // initialize temperatuere as zero on all nodes
	real dR;
	real dVf;

	// Find eigen position surface 
	Domain* d = Get_Domain(1); // Get domain pointer, update if different
	Thread* t_fixed = Lookup_Thread(d, eigen_face_zoneID); //pointer to fixed temp surface 

	face_t f;

	//real R_old = R; // Store old residual value 
	//real V_f_old = V_f; // Store old FSR value 

#if !RP_HOST
	begin_f_loop(f, t_fixed)
		if PRINCIPAL_FACE_P(f, t_fixed)
		{
			T_eig = F_T(f, t_fixed); //Temperture at x(eig)
		}
	end_f_loop(f, t_fixed)


	//R = T_eig - 1.2 * T_infty; //Residual
#endif
	
	// Sum temperature over compute nodes 
	// T_eig will be zero on all compute nodes except on the one that the face belongs to
	T_eig = PRF_GRSUM1(T_eig);

	// Calculate Residual
	R_old = R; // Set R_old to previous value
	R = T_eig - 1.2 * T_infty; //calculate new R

	//node_to_host_real_1(V_f); // update V_f on host process so report is correct.
	//node_to_host_real_1(Delta_Vf);
	node_to_host_real_1(R);// update R on host process so report is correct for residual report definition

	//Calculate Derivative if iteration counter is greater than zero 
	// otherwise old value is not valid
	if (N_ITER != 0)
	{
		dR = R - R_old;
		dVf = V_f - V_f_old;

		dR_by_dVf = dR / dVf;
		node_to_host_real_1(dR_by_dVf);
	}

	// Only converged if solution has been iterating enough and the temperature is within 1/1000th of a degree of the desired
	if (N_ITER > 20000 && fabs(dR) < 0.001)
	{
		converged_flag = true;
	}
	// Print Every 25 iterations to avoid excessive printing, update with different frequency if desired.
	// Count figure out how to automatically pass the profile update interval (count find a macro or rp var for it)
	if (N_ITER % 25 == 0)
	{
		Message0("Eigenvalue method Flame Spread Rate: %g m/s\n", V_f);
		Message0("Eigenvalue method temperture Residual R: %g K\n", R);
	}
} */

/* DEFINE_ADJUST(update_Vf_eigen, d)
{
#if !RP_HOST

    real dVf;
    real dR;
    real dV_update;

    
    //   Do not start eigen update immediately.
    //   Wait until the thermal solution has reached the eigen point.
    //   You can use either condition:
    //   - N_ITER >= eigen_start_iter
    //   - R > 0
    //   Here I use R > 0 as the main trigger, with the 1000-iteration condition
    //   as an additional safety delay.
    
    if (!eigen_update_started)
    {
        if ((N_ITER >= eigen_start_iter) && (R > 0.0))
        {
            eigen_update_started = true;
            waiting_for_response = true;

            R_sample = R;
            V_f_sample = V_f;
            V_f_old = V_f;

            V_f = V_f + Delta_Vf;

            last_vf_change_iter = N_ITER;

            Message0("Eigen update started at N_ITER=%d: R_sample=%g, V_f_sample=%g, perturbed V_f=%g\n",
                     N_ITER, R_sample, V_f_sample, V_f);
        }
        else
        {
            if (N_ITER % 100 == 0)
            {
                Message0("Eigen update not started: N_ITER=%d, R=%g, V_f=%g\n",
                         N_ITER, R, V_f);
            }
        }

        node_to_host_real_1(V_f);
        return;
    }

    
    //   If we just changed V_f, wait for the CFD/thermal field to respond.
    //   Do not calculate dR/dVf immediately.
    
    if (waiting_for_response)
    {
        if ((N_ITER - last_vf_change_iter) < vf_update_interval)
        {
            node_to_host_real_1(V_f);
            return;
        }

        dR = R - R_sample;
        dVf = V_f - V_f_sample;

        if (fabs(dVf) > 1.0e-20)
        {
            dR_by_dVf = dR / dVf;
            waiting_for_response = false;

            Message0("Eigen derivative sampled at N_ITER=%d: R=%g, R_sample=%g, dR=%g, dVf=%g, dR_by_dVf=%g\n",
                     N_ITER, R, R_sample, dR, dVf, dR_by_dVf);
        }
        else
        {
            Message0("Bad derivative sample: dVf=%g. Re-perturbing V_f.\n", dVf);

            R_sample = R;
            V_f_sample = V_f;
            V_f = V_f + Delta_Vf;

            last_vf_change_iter = N_ITER;
            waiting_for_response = true;
        }

        node_to_host_real_1(V_f);
        return;
    }

    
    //   Main damped Newton/secant update.
    //   This only happens after a valid derivative sample exists.
    
    if (fabs(dR_by_dVf) > 1.0e-12 && isfinite(dR_by_dVf))
    {
        V_f_old = V_f;

        dV_update = -eigen_relax * R / dR_by_dVf;

        if (dV_update > max_dVf)
            dV_update = max_dVf;
        else if (dV_update < -max_dVf)
            dV_update = -max_dVf;

        
        //   Store the current state before changing V_f.
        //  The next derivative will be computed after the flow responds.
        
        R_sample = R;
        V_f_sample = V_f;

        V_f = V_f + dV_update;

        
        //   Keep this floor. Negative flame spread rate is not useful unless you
        //  are explicitly allowing reverse spread.
        
        V_f = MAX(V_f, 0.0);

        last_vf_change_iter = N_ITER;
        waiting_for_response = true;

        Message0("Eigen update at N_ITER=%d: R=%g, dR_by_dVf=%g, dV_update=%g, new V_f=%g\n",
                 N_ITER, R, dR_by_dVf, dV_update, V_f);
    }
    else
    {
        Message0("Bad dR_by_dVf=%g. Taking new perturbation sample instead.\n",
                 dR_by_dVf);

        R_sample = R;
        V_f_sample = V_f;
        V_f = V_f + Delta_Vf;

        last_vf_change_iter = N_ITER;
        waiting_for_response = true;
    }

#endif

    node_to_host_real_1(V_f);
} */

// Update Solid Motion 
DEFINE_ZONE_MOTION(update_solid_motion, omega, axis, origin, velocity, current_time, dtime)
{
	NV_D(velocity, =, V_f, 0.0, 0.0); // Update solid motion velocity with calculated FSR
	//Message("Updated solid motion to %g m/s\n", V_f);
	return;
}

// Update Moving Wall BCs
DEFINE_PROFILE(update_wall_motion, thread, position)
{
	face_t f;

	begin_f_loop(f, thread)
	{
		F_PROFILE(f, thread, position) = V_f;// Update inlet velocity profile with calculated FSR
	}
	end_f_loop(f, thread)
		//Message("Updated wall motion to %g m/s\n", V_f);
}

// Create Report Definition for FSR
DEFINE_REPORT_DEFINITION_FN(flame_spread_rate)
{
	return V_f; // Return calculated flame spread rate for report definition
}


DEFINE_ON_DEMAND(check_rp_vars)
{
	bool U_mean_exists = RP_Variable_Exists_P("user/u_mean"); // Check if user-defined parameter for mean velocity exists
	bool alpha_exists = RP_Variable_Exists_P("user/alpha"); // Check if user-defined parameter for under-relaxation factor exists

	Message0("Checking for user-defined parameter 'user/u_mean': %d\n", U_mean_exists);
	Message0("Checking for user-defined parameter 'user/alpha': %d\n", alpha_exists);

	if (U_mean_exists)
	{
		real U_mean = RP_Get_Real("user/u_mean"); // Get mean velocity from user-defined parameter if it exists
		Message0("User-defined parameter 'user/u_mean' found with value: %f m/s\n", U_mean);
	}
	else
	{
		Message0("Warning: User-defined parameter 'user/u_mean' not found.\n");
	}

	if (alpha_exists)
	{
		alpha = RP_Get_Real("user/alpha"); // Get under-relaxation factor from user-defined parameter if it exists
		Message0("User-defined parameter 'user/alpha' found with value: %f\n", alpha);
	}
	else
	{
		Message0("Warning: User-defined parameter 'user/alpha' not found. Using default value: %f\n", alpha);
	}

}

DEFINE_ON_DEMAND(set_FSR)
{

	bool V_f_init_exists = RP_Variable_Exists_P("user/v_f_init"); // Check if user-defined parameter for initial FSR exists)

	Message0("Checking for user-defined parameter 'user/v_f_init': %d\n", V_f_init_exists);

	if (V_f_init_exists)
	{
		V_f = RP_Get_Real("user/v_f_init"); // Get initial FSR from user-defined parameter if it exists
		Message0("User-defined parameter 'user/v_f_init' found with value: %g m/s\n", V_f);
	}
	else
	{
		V_f = 0.0; // Default initial FSR value if user-defined parameter does not exist
		Message0("Warning: User-defined parameter 'user/v_f_init' not found. Using default value of 0 m/s.\n");
	}
	node_to_host_real_1(V_f); // update V_f on host process so report is correct.
}

/*=================================================================================
* Functions for calculating the FSR using the eigenposition-based method with 
* bisection based updates. This method uncouples the perscribed surface regression
* profile from the flame spread rate calculation. 
* To use this method, follow this process:
* 1. Define the zone ID of the face that serves as the Eigen position. Define a
*    cell register in fluent with a single cell that has a face on the fuel surface
*    Go to Domain > Zones > Separate > Faces > Mark > select the cell register and the
*    fuel surface zone. Then click separate. 
* 2. Set the user-defined rp variable "user/eigen_zone_id" to the zone ID of the 
*    separated face zone created in step 1.
* 3. Set the interval of the flame spread rate. The minimum should be zero. The maximum
*    interval should be chosen with care. Too high and the flame could be put out. Set 
*    the interval using the user-defined rp variables "user/vf_bracket_lower" and
*    "user/vf_bracket_upper". 
* 4. Set the number of iterations between FSR updates using the user-defined rp variable
*     "user/bisection_update_interval". 
* 5. Guess an initial flame spread rate if desired. (Default is zero). Set the rp variable
*    "user/v_f_init" to the desired initial FSR.
* 6. Run the on demand functions to set the desired values internally in the UDF. 
* 7. Set the DEFINE_EXECUTE_AT_END function hook to be calc_FSR_bisection.
* 8. Create a report definition for R and monitor it. 
* 9. Run
=================================================================================*/
DEFINE_ON_DEMAND(set_eigen_face_zoneID)
{
	//real T_eig;

	bool zone_ID_exists = RP_Variable_Exists_P("user/eigen_zone_id"); // Check if user-defined parameter for eigen face zone ID exists
	Message0("Checking for user-defined parameter 'user/eigen_zone_id': %d\n", zone_ID_exists);
	if (zone_ID_exists)
	{
		eigen_face_zoneID = RP_Get_Integer("user/eigen_zone_id"); // Get eigen face zone ID from user-defined parameter if it exists
		Message0("User-defined parameter 'user/eigen_zone_id' found with value: %d\n", eigen_face_zoneID);
	}
	else
	{
		eigen_face_zoneID = 0; // Default value if user-defined parameter does not exist, update with different default if desired
		Message0("Warning: User-defined parameter 'user/eigen_zone_id' not found. Using default value of 0.\n");
	}
	node_to_host_int_1(eigen_face_zoneID); // update eigen_face_zoneID on host process so it can be used in calc_FSR_eigen

}

// Function used to setup false position Vf solution procedure. Sets the initial bracket and the residual tolerance for accelpting a
// given R(V_f)
DEFINE_ON_DEMAND(setup_solve_Vf_false_position)
{
	bool lower_exists = RP_Variable_Exists_P("user/vf_bracket_lower"); // Check if user-defined parameter for lower FSR bracket exists
	Message0("Checking for user-defined parameter 'user/vf_bracket_lower': %d\n", lower_exists);
	bool upper_exists = RP_Variable_Exists_P("user/vf_bracket_upper"); // Check if user-defined parameter for upper FSR bracket exists
	Message0("Checking for user-defined parameter 'user/vf_bracket_upper': %d\n", upper_exists);

	bool R_bracket_lower_exists = RP_Variable_Exists_P("user/r_bracket_lower"); // Check if user-defined parameter for R bracket exists
	Message0("Checking for user-defined parameter 'user/r_bracket_lower': %d\n", R_bracket_lower_exists);
	bool R_bracket_upper_exists = RP_Variable_Exists_P("user/r_bracket_upper"); // Check if user-defined parameter for R bracket exists
	Message0("Checking for user-defined parameter 'user/r_bracket_upper': %d\n", R_bracket_upper_exists);

	bool residual_tolerance_exists = RP_Variable_Exists_P("user/false_position_residual_tolerance"); // Check if user-defined parameter for residual tolerance exists
	Message0("Checking for user-defined parameter 'user/false_position_residual_tolerance': %d\n", residual_tolerance_exists);

	// Check V_f bracket
	if (lower_exists)
	{
		V_f_bracket[0] = RP_Get_Real("user/vf_bracket_lower"); // Get lower FSR bracket from user-defined parameter if it exists
		Message0("User-defined parameter 'user/vf_bracket_lower' found with value: %g m/s\n", V_f_bracket[0]);
	}
	else
	{
		Message0("Warning: User-defined parameter 'user/vf_bracket_lower' not found. Using default value of %g m/s.\n", V_f_bracket[0]);
	}

	if (upper_exists)
	{
		V_f_bracket[1] = RP_Get_Real("user/vf_bracket_upper"); // Get upper FSR bracket from user-defined parameter if it exists
		Message0("User-defined parameter 'user/vf_bracket_upper' found with value: %g m/s\n", V_f_bracket[1]);
	}
	else
	{
		Message0("Warning: User-defined parameter 'user/vf_bracket_upper' not found. Using default value of %g m/s.\n", V_f_bracket[1]);
	}

	// Check R bracket
	if(R_bracket_lower_exists)
	{
		R_bracket[0] = RP_Get_Real("user/r_bracket_lower"); // Get lower R bracket from user-defined parameter if it exists
		Message0("User-defined parameter 'user/r_bracket_lower' found with value: %g K\n", R_bracket[0]);
	}
	else
	{
		Message0("Warning: User-defined parameter 'user/r_bracket_lower' not found. Using default value of %g K.\n", R_bracket[0]);
	}

	if(R_bracket_upper_exists)
	{
		R_bracket[1] = RP_Get_Real("user/r_bracket_upper"); // Get upper R bracket from user-defined parameter if it exists
		Message0("User-defined parameter 'user/r_bracket_upper' found with value: %g K\n", R_bracket[1]);
	}
	else
	{
		Message0("Warning: User-defined parameter 'user/r_bracket_upper' not found. Using default value of %g K.\n", R_bracket[1]);
	}

	if (residual_tolerance_exists)
	{
		false_position_residual_tolerance = RP_Get_Real("user/false_position_residual_tolerance"); // Get residual tolerance from user-defined parameter if it exists
		Message0("User-defined parameter 'user/false_position_residual_tolerance' found with value: %g K\n", false_position_residual_tolerance);
	}
	else
	{
		Message0("Warning: User-defined parameter 'user/false_position_residual_tolerance' not found. Using default value of %g.\n", 1e-5);
	}

	node_to_host_real(V_f_bracket,2); // update V_f_bracket on host process so it can be used in calc_FSR_eigen
	node_to_host_real(R_bracket,2); // update R_bracket on host process so it can be used in calc_FSR_eigen
	node_to_host_real_1(false_position_residual_tolerance); // update false_position_residual_tolerance on host process so it can be used in calc_FSR_eigen

	// Calculate the initial V_f using false position
	V_f = V_f_bracket[0]*R_bracket[1] - V_f_bracket[1]*R_bracket[0]; // False position method formula for new guess
	V_f = V_f / (R_bracket[1] - R_bracket[0]); // Divide by difference in residuals at bracket endpoints
	node_to_host_real_1(V_f); // update V_f on host process so it can be used in calc_FSR_eigen
	Message0("Initial FSR guess from false position method: %g m/s\n", V_f);
}

DEFINE_ON_DEMAND(set_bisection_update_interval)
{
	bool interval_exists = RP_Variable_Exists_P("user/bisection_update_interval"); // Check if user-defined parameter for bisection update interval exists
	Message0("Checking for user-defined parameter 'user/bisection_update_interval': %d\n", interval_exists);

	if (interval_exists)
	{
		N_UPDATE = RP_Get_Integer("user/bisection_update_interval"); // Get bisection update interval from user-defined parameter if it exists
		Message0("User-defined parameter 'user/bisection_update_interval' found with value: %d iterations\n", N_UPDATE);
	}
	else
	{
		Message0("Warning: User-defined parameter 'user/bisection_update_interval' not found. Using default value of %d iterations.\n", N_UPDATE);
	}

	node_to_host_int_1(N_UPDATE); // update N_UPDATE on host process so it can be used in calc_FSR_eigen

}

DEFINE_EXECUTE_AT_END(calc_FSR_false_position)
{
	// Access residuals and check to see if its time to update
	/* bool residuals_low =  */
	Domain* d = Get_Domain(1); // Get domain pointer, update if different
	if ((N_ITER - last_fsr_update_iter) < min_iter_after_fsr_update)
	{
		Message0("Skipping FSR update: only %d iterations since last V_f change. Need %d.\n",
				N_ITER - last_fsr_update_iter, min_iter_after_fsr_update);
		return;
	}


	real residual_tol = false_position_residual_tolerance; // Use value set in setup_solve_Vf_false_position, update if desired for convergence criteria
	//real residual_tol = RP_Get_Real("user/false_position_residual_tolerance"); // Get residual tolerance for accepting R values from user-defined parameter

	bool residuals_low = residuals_below(d, residual_tol); // Check if all residuals are below the specified tolerance, update threshold as needed for convergence criteria

	if (residuals_low)
	{
		// Initialize T_eig
		real T_eig = 0.0; // initialize temperatuere as zero on all nodes

		Thread* t_fixed = Lookup_Thread(d, eigen_face_zoneID); //pointer to fixed temp surface 

		face_t f;

#if !RP_HOST
		begin_f_loop(f, t_fixed)
			if PRINCIPAL_FACE_P(f, t_fixed)
			{
				T_eig = F_T(f, t_fixed); //Temperture at x(eig)
			}
		end_f_loop(f, t_fixed)
#endif
	
		// Sum temperature over compute nodes 
		// T_eig will be zero on all compute nodes except on the one that the face belongs to
		T_eig = PRF_GRSUM1(T_eig);

		// Calculate Residual
		R = T_eig - 1.2 * T_infty; //calculate new R
		node_to_host_real_1(R);// update R on host process so report is correct for residual report definition
	} 
	else
	{
		return; // Skip update if residuals are not low enough for convergence
	}

	// Compute new V_f
	/* real V_f_new;
	V_f_new = V_f_bracket[0]*R_bracket[1] - V_f_bracket[1]*R_bracket[0]; // False position method formula for new guess
	V_f_new = V_f_new / (R_bracket[1] - R_bracket[0]); // Divide by difference in residuals at bracket endpoints
	//node_to_host_real_1(V_f); // update V_f on host process so report is
	Message0("False position method update at iteration %d: R = %g K, new FSR = %g m/s\n", N_ITER, R, V_f); */

	if ( (R*R_bracket[0]) < 0)
	{
		V_f_bracket[1] = V_f; // If R has opposite sign of R at lower bracket endpoint, new FSR becomes upper bracket endpoint
		R_bracket[1] = R; // Update residual at upper bracket endpoint
	}
	else if ( (R*R_bracket[1]) < 0)
	{
		V_f_bracket[0] = V_f; // If R has opposite sign of R at upper bracket endpoint, new FSR becomes lower bracket endpoint
		R_bracket[0] = R; // Update residual at lower bracket endpoint
	}
	else
	{
		Message0("Warning: Residual does not change sign across bracket endpoints. Check for convergence or consider updating brackets.\n");
	}

	// Compute new V_f
	//real V_f_new;
	V_f = V_f_bracket[0]*R_bracket[1] - V_f_bracket[1]*R_bracket[0]; // False position method formula for new guess
	V_f = V_f / (R_bracket[1] - R_bracket[0]); // Divide by difference in residuals at bracket endpoints
	node_to_host_real_1(V_f); // update V_f on host process so report is

	last_fsr_update_iter = N_ITER;
	node_to_host_int_1(last_fsr_update_iter); // update last_fsr_update_iter on host process so it can be used in convergence checks or other UDFs if needed
	
	Message0("False Position ITERATION %d: \n", N_FSR_ITER);
	Message0("False position method update at iteration %d: R = %g K, new FSR = %g m/s\n", N_ITER, R, V_f);
	Message0("Current FSR bracket: [%g, %g] m/s\n", V_f_bracket[0], V_f_bracket[1]);
	Message0("Current R bracket:   [%g, %g] K\n", R_bracket[0], R_bracket[1]);

	N_FSR_ITER++;
	node_to_host_int_1(N_FSR_ITER); // update N_FSR_ITER on host process so it can be used in reports or other UDFs if needed
}

DEFINE_EXECUTE_AT_END(calc_FSR_control)
{
	if (N_ITER % N_UPDATE ==0)
	{
		// Initialize T_eig
		real T_eig = 0.0; // initialize temperatuere as zero on all nodes

		// Find eigen position surface 
		Domain* d = Get_Domain(1); // Get domain pointer, update if different
		Thread* t_fixed = Lookup_Thread(d, eigen_face_zoneID); //pointer to fixed temp surface 

		face_t f;

#if !RP_HOST
		begin_f_loop(f, t_fixed)
			if PRINCIPAL_FACE_P(f, t_fixed)
			{
				T_eig = F_T(f, t_fixed); //Temperture at x(eig)
			}
		end_f_loop(f, t_fixed)
#endif
	
		// Sum temperature over compute nodes 
		// T_eig will be zero on all compute nodes except on the one that the face belongs to
		T_eig = PRF_GRSUM1(T_eig);

		// Calculate Residual
		R = T_eig - 1.2 * T_infty; //calculate new R
		node_to_host_real_1(R);// update R on host process so report is correct for residual report definition
	}
	else 
	{
		return; // Skip update if not at the correct iteration
	}

	real R_tolerance = 0.001; // Set residual tolerance for convergence, update as needed
	real delta_Vf = 0.00000001; // Set FSR adjustment step size, update as needed 0.1 um/s

	// increment V_f based on sign of R
	if (fabs(R) < R_tolerance)
	{
		Message0("Control method converged at iteration %d with FSR = %g m/s and Residual R = %g K\n", N_ITER, V_f, R);
		return; // Converged, no update needed
	}

	if (R > 0)
	{
		V_f = V_f + delta_Vf; // If too hot (positive residual), FSR is too small and more cooling is needed. Increase FSR.
	}
	else	
	{
		V_f = V_f - delta_Vf; // If too cold (negative residual), FSR is too large and less cooling is needed. Decrease FSR.

		V_f = MAX(V_f, 0.0); // Ensure FSR does not become negative
	}
	node_to_host_real_1(V_f); // update V_f on host process so report is correct for next iteration
	Message0("Control method update at iteration %d: R = %g K, new FSR = %g m/s\n", N_ITER, R, V_f);
}

DEFINE_ON_DEMAND(residual_list)
{
	int nw;
	real scaled_res;
	Domain *domain=Get_Domain(1);

	for(nw=0; nw<DOMAIN_NUMEQN(domain); ++nw)
	{
		int eqn = DOMAIN_EQNS(domain, nw); // Get the equation index for the current equation/species
		
		if(strlen(DOMAIN_EQN_LABEL(domain, eqn))>0)
		{
			if(0==DOMAIN_RES_SCALE(domain,nw)) DOMAIN_RES_SCALE(domain,nw)=1.0;

			scaled_res=DOMAIN_RES(domain,nw)/DOMAIN_RES_SCALE(domain,nw);
			Message0("%s equation,residual=%g\n",DOMAIN_EQN_LABEL(domain,nw ),scaled_res);
		}
	}
}

// Report Definition for Residual of eigenvalue-based FSR calculation
DEFINE_REPORT_DEFINITION_FN(residual_R)
{
	return R; // Return residual for report definition
}
/*======================================================================================*/
// I think this is old and is unused 
DEFINE_INIT(init_R, d)
{
	// Find fixed temp thread
	real T_eig;

	if (eigen_face_zoneID == 0)
	{
		Message0("Error: eigen_face_zoneID is not set. Please set it using the 'user/eigen_zone_id' RP variable before initializing.\n");
		Message0("Warning: Residual R for calc_FSR_eigen not initialized!\n");
		return;
	}

	Thread* t_fixed = Lookup_Thread(d, eigen_face_zoneID); //pointer to fixed temp surface 

	face_t f;

#if !RP_HOST
	begin_f_loop(f, t_fixed)
		if PRINCIPAL_FACE_P(f, t_fixed)
		{
			T_eig = F_T(f, t_fixed); //Temperture at x(eig)
			Message0("Initial temperature at eigenposition: %g K\n", T_eig);
		}
	end_f_loop(f, t_fixed)

		R = T_eig - 1.2 * T_infty; //Residual

#endif
	node_to_host_real_1(R);	// update R on host process so it can be used in calc_FSR_eigen
	Message0("Initialized and passed residual R = %g K to host node\n", R);
	
}


DEFINE_RW_FILE(write_FSR, fp)
{
	Message0("Writing FSR to file: %g m/s\n", V_f);
#if !RP_NODE
	fprintf(fp, "%g", V_f); // Write FSR value to file
#endif
}

DEFINE_RW_FILE(read_FSR, fp)
{
	Message0("Reading FSR from file...\n");
#if !RP_NODE
	fscanf(fp, "%g", &V_f);
#endif
}

DEFINE_RW_HDF_FILE(write_FSR_hdf, filename)
{

	char* path = "/FSR_data"; // HDF5 dataset path for FSR data
	real* data_ptr = &V_f; // Pointer to FSR data to write
	size_t nelems = 1; // Number of elements to write (1 in this case since we're writing a single value)

	Write_Complete_User_Dataset(filename, path, data_ptr, nelems);

}

DEFINE_RW_HDF_FILE(read_FSR_hdf, filename)
{
	char* path = "/FSR_data"; // HDF5 dataset path for flame spread rate data
	real* data_ptr = &V_f; // Pointer to FSR variable to read into
	size_t nelems = 1; // Number of elements to read (1 in this case since we're reading a single value)

	Read_Complete_User_Dataset(filename, path, data_ptr, nelems);
}