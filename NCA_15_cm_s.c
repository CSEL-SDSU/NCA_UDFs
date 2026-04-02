
 #include "udf.h"



/* CODE SECTION */
/* FD INLET VELOCITY PROFILE */ 

 DEFINE_PROFILE(inlet_x_vel_15cms, thread, position) 
 {
    real x[ND_ND]; /* this will hold the position vector */
    real y, h, U_mean, U_max, m, n;
    face_t f;

	h = 0.00495; /* m; inlet height, do not change */
	m = 27.59596236; /* constant, do not change */ 
	n = 2.0; /* constant, do not change */

	U_mean = 0.082; /* m/sec; inlet mean velocity, update with geom */
	U_max = U_mean*((m+1)/m)*((n+1)/n); /* m/sec; max velocity, at centerline... calc */ 

    begin_f_loop(f,thread)
    {
      F_CENTROID(x, f, thread);
      y = 2.*(x[1]-0.5*h)/h; /* non-dimensional y coordinate, b/c coord sys is at bottom of geom not centerline... calc */

      F_PROFILE(f, thread, position) = U_max*(1.0-y*y); /* m/sec; velocity as f(y) at centerline... calc */
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

	 Thread* t;         

	 real mdot_chem = 0.; //Mass flux from chemical reaction at surface [kg/s]
     real V_f;

     face_t f; // Face along surface

	 // Loop through faces along surface and sum mass flux from chemical reaction
     begin_f_loop(f, t)
     {
         // Evem though this macro is labeled as "FLUX", it is actually a mass flow rate through a face 
         // according to section 3.2.2.4 of Fluent Customization manual. So we can sum this value across
         // all faces along the surface to get total mass flow rate from chemical reaction at surface.

		 mdot_chem += F_FLUX(f, thread); // Sum mass fluxes on each face from chemical reaction at surface
     }
     end_f_loop(f, t)

	 // Calculate corrected FSR and print result
	 V_f = mdot_chem / (rho * I); // Calculate flame spread rate [m/s]

     printf("Mass flux from chemical reaction at surface: %g kg/s\n", mdot_chem);
	 printf("Calculated Flame Spread Rate: %g m/s\n", V_f);

	 // Update Solid Motion and Moving Wall BCs with calculated FSR
	 // Update U_mean in inlet velocity profile with calculated FSR
 }