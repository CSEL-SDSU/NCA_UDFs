
 #include "udf.h"
#include <stdbool.h>


/* CODE SECTION */
/* FD INLET VELOCITY PROFILE */ 
/*
	Modified version of parabolic.c to use RP variable for mean velocity.
    Put you velocity in the brackets (remove the brackets)
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
    bool U_mean_exists = RP_Variable_Exists_P("user/U_mean"); // Check if user-defined parameter exists
    Message("Checking for user-defined parameter 'user/U_mean': %d\n", U_mean_exists);

    if (RP_Variable_Exists_P("user/U_mean"))
    {
        U_mean = RP_Get_Real("user/U_mean"); // Get mean velocity from user-defined parameter if it exists
    }
    else
    {
        Message("Warning: User-defined parameter 'user/U_mean' not found. Using default value of %f m/s.\n", U_mean);
    }


	U_max = U_mean*((m+1)/m)*((n+1)/n); /* m/sec; max velocity, at centerline... calc */ 

    begin_f_loop(f,thread)
    {
      F_CENTROID(x, f, thread);
      y = 2. * (x[1] - 0.5 * h) / h; /* non-dimensional y coordinate, b/c coord sys is at bottom of geom not centerline... calc */

      F_PROFILE(f, thread, position) = U_max*(1.0 - (y * y)); /* m/sec; velocity as f(y) at centerline... calc */
    }
    end_f_loop(f, thread)
 }

 // Calculate Flame Spread Rate
 DEFINE_ADJUST(calc_FSR, d)
 {

     // Integral of the x-component of the normal vector f'(x)/sqrt(f'(x)^2+1
     // Units of I are m^2. Even though the integral is over one spatial direction, x
     // it obtains units of m^2 from multiplication by the reference length (depth) of 1 m

     const real I = 0.002243881756148; // [m^2] Computed in Matlab, changes with different surface profiles. Corresponds to Hossain V_g = 8.2 cm/s curve
     const real rho = 1190; // [kg/m^3] Density of solid phase

     // Find wall_mass_flux thread
	 int zone_ID = 5; // ID of surface zone where chemical reaction occurs, update if different. Zone is shown in Boundary conditions tab
	 Thread* t = Lookup_Thread(d, zone_ID); // Get thread pointer for surface zone where chemical reaction occurs

	 real mdot_chem = 0.; //Mass flux from chemical reaction at surface [kg/s]
     real V_f;

     face_t f; // Face along surface

	 // Loop through faces along surface and sum mass flux from chemical reaction
     begin_f_loop(f, t)
     {
         // Evem though this macro is labeled as "FLUX", it is actually a mass flow rate through a face 
         // according to section 3.2.2.4 of Fluent Customization manual. So we can sum this value across
         // all faces along the surface to get total mass flow rate from chemical reaction at surface.

		 mdot_chem += F_FLUX(f, t); // Sum mass fluxes on each face from chemical reaction at surface
     }
     end_f_loop(f, t)

	 // Calculate corrected FSR and print result
	 V_f = mdot_chem / (rho * I); // Calculate flame spread rate [m/s]

     printf("Mass flux from chemical reaction at surface: %g kg/s\n", mdot_chem);
	 printf("Calculated Flame Spread Rate: %g m/s\n", V_f);

	 // Update Solid Motion and Moving Wall BCs with calculated FSR
	 // Update U_mean in inlet velocity profile with calculated FSR
 }