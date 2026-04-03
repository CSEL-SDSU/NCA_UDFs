
#include "udf.h"

static real V_f = 0; // Initialize flame spread rate variable, will be updated at end of each iteration in calc_FSR and used in inlet velocity profile and solid motion BCs

/* CODE SECTION */
/* FD INLET VELOCITY PROFILE */

DEFINE_PROFILE(inlet_x_vel_8cms, thread, position)
{
	real x[ND_ND]; /* this will hold the position vector */
	real y, h, U_mean, U_max, m, n;
	face_t f;

	h = 0.00495; /* m; inlet height, do not change */
	m = 27.59596236; /* constant, do not change */
	n = 2.0; /* constant, do not change */

	U_mean = 0.082 + V_f; /* m/sec; inlet mean velocity, update with geom */
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
	V_f = mdot_chem / (rho * I); // Calculate flame spread rate [m/s]
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
	Thread* t = Lookup_Thread(d, zone_ID); // Get thread pointer for surface zone where chemical reaction occurs

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

	// Sum mdot over all compute nodes
	mdot = PRF_GRSUM1(mdot);

	// Sum integrals over all compute nodes
	I_num = PRF_GRSUM1(I_num);
	I_denom = PRF_GRSUM1(I_denom);

	// Calculate corrected FSR and print result
	//V_f = mdot_chem / (rho * I); // Calculate flame spread rate [m/s]
	V_f = I_num / (rho * I_denom); // Calculate flame spread rate using least squares method [m/s]
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