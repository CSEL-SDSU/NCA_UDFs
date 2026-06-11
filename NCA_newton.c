
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
static real Delta_Vf = 0.000010; // Initial change in FSR
static const real T_infty = 300; // Ambient temperature [K]

static real residual_tolerance = 1e-6; // Residual tolerance for accepting a given R(V_f) in false position method, adjust as needed

static real relative_change_V_f = 0.000001; // Stop changing V_f is the step is only 1 um/s
static int N_NEWTON_ITER = 0;
static int N_LAST_UPDATE = -1;
static const int MIN_ITER_BEFORE_UPDATE = 1000;

// State flags
static int LSQ_mode = 1;
static int newton_mode = 0;

#define R_NP 1000
#define R_HISTORY_SIZE (R_NP + 1)

static real R_relative_tolerance = 1e-6;
static real R_denominator_floor = 1e-12;

static real R_history[R_HISTORY_SIZE];
static int R_history_count = 0;
static int R_history_head = 0;

static real R_relative_residual_max = 1e20;



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

// Function to clear the storage of R
static void reset_R_history(void)
{
#if !RP_HOST
	R_history_count = 0;
	R_history_head = 0;
	R_relative_residual_max = 1e20;
#endif
	node_to_host_int_2(R_history_count,R_history_head);
	node_to_host_real_1(R_relative_residual_max);
}

// Function to write current R into R_history array
static void store_R_in_history(real R_new)
{
#if !RP_HOST
	R_history[R_history_head] = R_new;

	R_history_head++;

	// If we exceed the numeber of stored iterations, write the next one to the
	// start of the array
	if (R_history_head >= R_HISTORY_SIZE)
	{
		R_history_head = 0;
	}

	// if we have stored less than the storage size increase the count
	// of stored R's 
	if (R_history_count < R_HISTORY_SIZE)
	{
		R_history_count++;
	}
#endif

	node_to_host_int_2(R_history_head,R_history_count);
}

// function to check R relative change. Determines whether the largest R
// residual over the last 1000 iterations is less than the tolerance 
static bool R_relative_change_below(real threshold)
{
	int R_converged = 0;

#if !RP_HOST
	int k; 
	int idx_curr;
	int idx_old;
	real R_curr;
	real denom;
	real res_k;

	R_relative_residual_max = 1e20;

	// if we've counted enough values check for convergence 
	if (R_history_count >= R_HISTORY_SIZE)
	{

		// Get the index to the current R. R_history_head always points to
		// the next place to write so the one before it will be the current
		idx_curr = R_history_head - 1;

		//if were at the end of storage and the next place to write is idx =0,
		// the current is the other end.
		if (idx_curr < 0)
		{
			idx_curr += R_HISTORY_SIZE;
		}

		R_curr = R_history[idx_curr];

		denom = MAX(fabs(R_curr), R_denominator_floor);

		R_relative_residual_max = 0.0;

		for (k = 1; k <= R_NP; k++)
		{
			idx_old = idx_curr - k;
			while (idx_old < 0)
			{
				idx_old += R_HISTORY_SIZE;
			}

			res_k = fabs(R_curr - R_history[idx_old]) / denom;
			R_relative_residual_max = MAX(R_relative_residual_max, res_k);
		}

		R_converged = (R_relative_residual_max < threshold);
	}
#endif

	node_to_host_int_1(R_converged);
	node_to_host_real_1(R_relative_residual_max);

	return R_converged ? true : false;
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
// Use LSQ Method for the initial guess for Newton's method
//Should be above update_FSR_newton in function handle list
DEFINE_EXECUTE_AT_END(update_FSR_LSQ)
{
	if (!LSQ_mode)
	{
		return; //Dont do LSQ update if not in LSQ mode.
	}

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

//Make sure this function is above 'update_FSR_newton' in define at end list
//Will be executed every iteration
DEFINE_EXECUTE_AT_END(update_R)
{
	Domain* d = Get_Domain(1); // Get domain pointer, update if different

	// Calculate R(V_f)
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
	R = T_eig - 1.2 * T_infty; // calculate new R
	node_to_host_real_1(R);    // update R on host process so report is correct

	store_R_in_history(R);

	if (N_ITER % 25 == 0)
	{
		Message0("NEWTON_ITER=%d, V_f = %g m/s, R = %f K, R_rel_res_max = %g, R_hist_count = %d\n",
				N_NEWTON_ITER, V_f, R, R_relative_residual_max, R_history_count);
	}
}

DEFINE_EXECUTE_AT_END(update_FSR_newton)
{
	
	
	// Check residuals. Only do a newton update if R(V_f) has been safely evaluated.
	//Domain* d = Get_Domain(1); // Get domain pointer, update if different

	//if residuals are not below tolance or if enough iterations have not passed since last update, return
	//if (!residuals_below(d, residual_tolerance) || !((N_ITER - N_LAST_UPDATE) > MIN_ITER_BEFORE_UPDATE))
	//{
	//	return; // Skip FSR update 
	//}
	// Only do a Newton update after R(V_f) has become stationary
	// over the last R_NP iterations.
	if (!R_relative_change_below(R_relative_tolerance) ||
		!((N_ITER - N_LAST_UPDATE) > MIN_ITER_BEFORE_UPDATE))
	{
		return; // Skip FSR update
	}

	// When residuals are low enough using LSQ mode switch to newton mode
	if (N_NEWTON_ITER == 0 && LSQ_mode)
	{
		LSQ_mode = 0; // Turn off LSQ mode so it doesn't interfere with Newton updates
		newton_mode = 1; // Turn on Newton mode to allow for Newton updates in this function

		// Update state flages on host nodes. This function runs on compute nodes
		node_to_host_int_1(LSQ_mode); 
		node_to_host_int_1(newton_mode);

		Message0("Switching from LSQ to Newton. V_f = %g m/s\n", V_f);
	}

	//V_f(1) and R(1) stored in V_f and R.
	//V_f(0) and R(0) stored in V_f_old and R_old


	if (newton_mode)
	{
		//V_f_old = V_f;
		//R_old = R; 
		//node_to_host_real_1(V_f_old);
		//node_to_host_real_1(R_old);

		// Now we have a set of (V_f, R(V_f)) points that are a valid evaluation
		if (N_NEWTON_ITER == 0)
		{
			// We only have V_f and R at N_NEWTON_ITER = 0. We do not have old values.
			// This is sort of a half iteration to change V_f and get a second evaluation R
			// First Newton iteration perturb with Delta_Vf. 
			// In this branch we basically take the initial V_f from LSQ and calculate what R is for it
			
			V_f_old = V_f;
			R_old = R; 
			node_to_host_real_1(V_f_old);
			node_to_host_real_1(R_old);

			// V_f(1) = V_f(0) + deltaV_f
			V_f += Delta_Vf;
			N_NEWTON_ITER++; 

			N_LAST_UPDATE = N_ITER; //Save iteration number of the update

			reset_R_history();

			Message0("NEWTON_ITER=%d, V_f = %g m/s. Updated V_f, solving R(V_f).\n ", N_NEWTON_ITER, V_f);
			node_to_host_real_1(V_f);
			node_to_host_int_1(N_NEWTON_ITER);
			node_to_host_int_1(N_LAST_UPDATE);
			return;
		}
		else
		{
			//V_f(n) and R(n) are stored in V_f and R. V_f(n-1) and R(n-1) are stored in V_f_old and R_old
			//Now we have the current evaluation V_f(1),R(V_f(1)) or V_f(N), R(V_f(N)) and the evaluation at the previous
			//iteration or the initial guess (V_f_LSQ, R(V_f_LSQ))

			//Calculate V_f
			dR_by_dVf = (R - R_old) / (V_f - V_f_old); // Numerical derivative dR/dVf using finite difference from previous iteration

			// make sure derivative is not zero (in that case procedure fails)
			if (fabs(dR_by_dVf) < 1e-12) { 
				Message0("Warning: Newton procedure failed dR_by_dVf = %g \n",dR_by_dVf);
				newton_mode=0;  
				node_to_host_int_1(newton_mode);
				return;
			}
			
			// Save old values for the next iteration
			V_f_old = V_f;
			R_old = R;
			node_to_host_real_1(V_f_old);
			node_to_host_real_1(R_old);

			// Update V_f 
			// V_f(n+1) = R(n)*( (V_f(n)) - V_f(n-1))/(R(n) - R(n-1)) )
			V_f = V_f_old - R / dR_by_dVf; 
			N_NEWTON_ITER++; 
			
			N_LAST_UPDATE = N_ITER; //Save iteration number of the update

			reset_R_history();

			node_to_host_real_1(V_f);
			node_to_host_int_1(N_NEWTON_ITER);
			node_to_host_int_1(N_LAST_UPDATE);

			Message0("NEWTON_ITER=%d, V_f = %g m/s. Updated V_f, solving R(V_f).\n ", N_NEWTON_ITER, V_f);

			if (fabs(V_f - V_f_old) <= relative_change_V_f)
			{
				Message0("Change in V_f is only 1 um/s. Procedure completed successfully. NEWTON_ITER=%d, V_f = %g m/s \n", N_NEWTON_ITER, V_f);
				newton_mode = 0;
				node_to_host_int_1(newton_mode);
				return;
			}

			Message0("NEWTON_ITER=%d, V_f = %g m/s. Updated V_f, solving R(V_f).\n ", N_NEWTON_ITER, V_f);

		}
	}
}

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

DEFINE_REPORT_DEFINITION_FN(residual_R_relative_change)
{
	return R_relative_residual_max;
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

#if !RP_NODE // get the V_f value on from the host
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
#endif

	//node_to_host_real_1(V_f); // update V_f on host process so report is correct.
	host_to_node_real_1(V_f); // broadcast V_f to nodes
	Message0("FSR initialized to V_f = %g m/s \n", V_f);
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
#if !RP_NODE //get desired value from host
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
#endif
	//node_to_host_int_1(eigen_face_zoneID); // update eigen_face_zoneID on host process so it can be used in calc_FSR_eigen
	host_to_node_int_1(eigen_face_zoneID); //broadcast variable to nodes

	Message0("eigen_face_zoneID initialized to %d \n", eigen_face_zoneID);
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
/*=====================================================================================*/



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