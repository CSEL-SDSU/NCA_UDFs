// This directive tells the compiler to expose gnu specific features, like qsort_r
#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include "udf.h"
#include "hdfio.h"

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define EPS 2.2204460492503131e-16
#define sgn(x)  ((x>0) - (x<0))
#define DISPLACEMENT_THRESHOLD 1e-5 // threashold for surface convergence [m]
/*--- PID-control FSR calculation globals---*/
static real V_f; //Flame Spread Rate 
static const int UPDATE_INTERVAL = 1; // Number of iterations between FSR updates
static real Kp = 1e-6; // Proportional gain for P-control update of FSR
static real Ki = 7.5e-9; // Integral gain control
static real Kd = 1e-7; // Derivative gain control

/*
 * Controller bias and integral state.
 *
 * At steady state:
 *     R = 0
 *     dRdn = 0
 *     V_f = V_f_bias + I_term
 */
static real V_f_bias = 40e-6;
static real I_term = 0.0; 

static const real R_integral_enable = 10; //abs of R to accumulate integral to prevent windup, prob make rp variable

//Variables required for evaluating F(X) - The residual temperature
// And determining when an evaluation is accepted
static real R;
static real R_old; 
static real dRdn = 0;
static real dRdn_old = 0;
static const real alpha = 0.1667; //Derivative term under-relaxation
static int eigen_face_zoneID = -1; // Zone ID of separated surface where temp is monitored
static const real T_infty = 300; 
static real T_eig_setpoint = 360; // Pseudo fixed temperature of eigen face

static const real I_min = -55e-6;
static const real I_max =  35e-6;

// Global regressed surface coordinates. Updated when calc_regressed_surf is called
// DO NOT USE ON HOST NODE NOT ALLOCATED
// Before install is called calc_regression owns x_f and y_f_new. After, the memory is owned by
// the global state x_f_g and y_f_g and x_f and y_f_new are NULL.

static size_t size_g; 
static real *x_f_g = NULL, *y_f_new_g = NULL; 
static bool valid_profile = false;
static real r_eig_g[ND_ND] = {0.0, 0.0};

// Dynamic grid update variables
static int N_MESH_UPDATE_IGNORE = 5000; //# of iterations to ignore before adjusting grid
static int N_MESH_UPDATES = 0; // Counter for number of mesh updates performed

// Dynamic Grid Under-relaxation
static real alpha_grid = 0.5; // Under-relaxation factor for grid motion, 0.5 is a good starting point

// Dynamic grid largest displacement residual. Should get smaller with each update
static real largest_displacement = 0.0; 
static int largest_displacement_conv = 0;

static void cleanup_profile(void)
{
	free(x_f_g);
	free(y_f_new_g);

	x_f_g = NULL;
	y_f_new_g = NULL;
	size_g = 0;
	valid_profile = false;
}

// This function takes the x_f and y_f_new arrays and assignes the global variables to point to their memory instead
// This allows all functions to access those values now without having to allocate more memory.
static int install_profile(real **x_src, real **y_src, size_t n)
{
	// Pass arguemtns as double pointers so that we may operate on the callers copy
	if (x_src == NULL || y_src == NULL || *x_src == NULL || *y_src == NULL)
	{
		// If the inputs are null or point to null, fail the install
		return 1;
	}

	// Otherwise clear the profile and install the new one
	cleanup_profile();

	x_f_g = *x_src;
	y_f_new_g = *y_src; 
	size_g = n;
	valid_profile = true;

	*x_src = NULL;
	*y_src = NULL;

	return 0; //Success
}



/*--- Comparator function for sorting face centroid coordinates ---*/
// Comparator function for qsort. 
// Should return:
// < 0 (negative) if a should be before b
// 0 if a and b are equal 
// > 0 (positive) if a should be after b)
static int compare_xf(void *context, const void* a, const void* b)
{

	// Cast arguemnts as integer pointers and dereference to get value
	int idx_a = *(const int*)a; 
	int idx_b = *(const int*)b;

	// Cast contest array as a real pointer
	real *x_f = (real*)context;

	// Get x_f values at the indicies
	real x_a = x_f[idx_a];
	real x_b = x_f[idx_b];

	return (x_a > x_b) - (x_b > x_a);
}

// If in linux (on HPC cluster) use qsort_r friendly wrapper.
#ifdef __linux__
static int compare_xf_linux(const void *a, const void *b, void *context)
{
	return compare_xf(context, a, b);
}
#endif

/*-----------------------------------------------------------------*/

/*----1D Interpolation/Table lookup function ---
	Used for interpolating nodal coordinates from new face centroid 
 	coordinates. Assumes that the data arrays, x and y are properly
 	sorted in ascending order.
 	Arguments:
	x: array of x data points
	y: array of y data points
	count: number of elements in x and y
	xq: query point to interpolate to
	extrap: boolean to flag linear extrapolation or not
	Returns: real value yq, interpolation of the table data at xq
--------------------------------------------------------------------*/
static real interp1d(const real *x, const real *y, size_t count, real xq, bool extrap)
{
	real yq, f0, f1, x0, x1, a1;

	// Check table bounds, linear extrapolation if out of bounds 
	if (xq <= x[0]) 
	{	
		if (extrap)
		{
			x0 = x[0];
			x1 = x[1];
			f0 = y[0];
			f1 = y[1];
			a1 = (f1 - f0) / (x1 - x0);

			yq = f0 + (xq - x0) * a1;
			return yq;
		}
		else
		{
			return y[0];
		}
	}
	if (xq >= x[count - 1]) 
	{
		if (extrap)
		{
			x0 = x[count - 2];
			x1 = x[count - 1];
			f0 = y[count - 2];
			f1 = y[count - 1];
			a1 = (f1 - f0) / (x1 - x0);

			yq = f0 + (xq - x0) * a1;
			return yq;
		}
		else
		{
			return y[count - 1];
		}
		
	} 

	// Binary search for the x interval that xq lies in
	size_t idx_upper = count - 1;
	size_t idx_lower = 0;
	// size_t idx_mid;

	while (idx_upper - idx_lower > 1)
	{
		size_t idx_mid = idx_lower + (idx_upper - idx_lower) / 2; // Middle index of bracket

		if (xq > x[idx_mid])
		{
			idx_lower = idx_mid;
		}
		else
		{
			idx_upper = idx_mid;
		}
	}

	// do linear interpolation
	x0 = x[idx_lower];
	x1 = x[idx_upper];
	f0 = y[idx_lower];
	f1 = y[idx_upper];
	a1 = (f1 - f0) / (x1 - x0);

	yq = f0 + (xq - x0) * a1;

	// return interpolated value
	return yq;
}
// Function to clear the storage of R

/* FD INLET VELOCITY PROFILE */
/*
	This funciton was originally written by Patterson, modified by Del Rio Blanco,
	and further modified by Volponi. Original correlations come from Shah and 
	London p. 197.
	Modified version of NCA.c to use RP variable for mean velocity.
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

	// U_mean += V_f; // Add calculated FSR to mean velocity for inlet profile

	U_max = U_mean * ((m + 1) / m) * ((n + 1) / n); /* m/sec; max velocity, at centerline... calc */

	begin_f_loop(f, thread)
	{
		F_CENTROID(x, f, thread);
		y = 2. * (x[1] - 0.5 * h) / h; /* non-dimensional y coordinate, b/c coord sys is at bottom of geom not centerline... calc */

		F_PROFILE(f, thread, position) = U_max * (1.0 - (y * y)) + V_f; /* m/sec; velocity as f(y) at centerline... calc */
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

//Make sure this function is above 'update_FSR_Pcontrol' in define at end list
//Will be executed every iteration
DEFINE_EXECUTE_AT_END(update_R)
{
	//real T_eig_sum = 0.0;
	real T_eig = 0.0;
	int n_eig_faces = 0;

#if !RP_HOST //will run on all compute nodes or in serial

	// Domain, thread, and face variables are not availible on the host
	Domain* d = Get_Domain(1); // Get domain pointer, update if different

	// Calculate R(V_f)
	//real T_eig = 0.0; // initialize temperatuere as zero on all nodes

	Thread* t_fixed = Lookup_Thread(d, eigen_face_zoneID); //pointer to fixed temp surface 

	face_t f;

	begin_f_loop(f, t_fixed)
		if PRINCIPAL_FACE_P(f, t_fixed)
		{
			T_eig = F_T(f, t_fixed); //Temperture at x(eig)
			n_eig_faces++;
		}
	end_f_loop(f, t_fixed)
#endif
	
	// Sum temperature over compute nodes 
	// T_eig will be zero on all compute nodes except on the one that the face belongs to
#if RP_NODE //only needed in parallel. Will only run in parallel
	T_eig = PRF_GRSUM1(T_eig);
	n_eig_faces = PRF_GISUM1(n_eig_faces);
#endif

	// Calculate Residual. Use synced T_eig value if parallel. Do not run on host if in parallel b/c host has T_eig = 0
#if !RP_HOST
	if (N_ITER > 2)
	{
		R_old = R;
		dRdn_old = dRdn; // old derivative

		// New R and dRdn
		R = T_eig - T_eig_setpoint;
		//R = T_eig - 1.2 * T_infty;
		dRdn = R - R_old; // raw derivative
		dRdn = alpha * (dRdn) + (1 - alpha)*dRdn_old;
	}
	else if (N_ITER > 1)
	{
		R_old = R;
		//R = T_eig - 1.2 * T_infty;
		R = T_eig - T_eig_setpoint;

		dRdn = R - R_old; // raw derivative
		
		dRdn_old = dRdn; // when we dont have dRdn_old, use raw derivatibe
		dRdn = alpha * (dRdn) + (1 - alpha) * dRdn_old;
	}
	else 
	{
		//R = T_eig - 1.2 * T_infty;
		R = T_eig - T_eig_setpoint;
	}
#endif

	// update R, Tigno_eig, and face count on host process so report is correct
	node_to_host_real_5(R,T_eig,R_old, dRdn, dRdn_old);    
	node_to_host_int_1(n_eig_faces);

	// Store R
	// store_R_in_history(R);

	// Report using host process or single serial process 
#if !RP_NODE
	if (n_eig_faces != 1 && N_ITER % 25 == 0)
	{
		Message("Warning: eigen face zone %d has %d principal faces. T_eig = %f K\n",
			eigen_face_zoneID, n_eig_faces, T_eig);
	}

	if (N_ITER % 25 == 0)
	{
		Message("Iteration %d: R = %g K, dR/dn = %g, T_eig = %g K, V_f = %g m/s\n", N_ITER, R, dRdn, T_eig, V_f);
	}
#endif
}

DEFINE_EXECUTE_AT_END(update_FSR_PIDcontrol)
{

	int return_flag = 0;

#if !RP_NODE //run on host in parallel or on single process in serial
	// Check if it is time for an update
	if (N_ITER % UPDATE_INTERVAL != 0)
	{
		return_flag = 1;
	}
#endif

	host_to_node_int_1(return_flag);	
	if(return_flag) {return;}

#if !RP_NODE
	
	// Flags to mark integral accumulation
	bool above_floor;
	bool R_move_away_floor; 

	real I_candidate = I_term; 

	// Prevent integral from accumulating when R is negative but V_f would not change b/c of floor
	// Ex in ignition before flame is close to eigen position, R is negative, but flame is actually
	// moving forward, so integral should not accumulate
	if (fabs(R) <= R_integral_enable)
	{
		I_candidate = I_term + Ki * R;

		// Clamp integral term
		I_candidate = MAX(I_min, MIN(I_max, I_candidate));

		// Evaulate candidate output
		real V_f_candidate = V_f_bias + Kp * R + I_candidate + Kd * dRdn;


		if (V_f_candidate > 0.0 || R > 0.0)
		{
			I_term = I_candidate; 
		}
	}

	// Update FSR using P-control based on the residual temperature R
	// R = T(x_eig) - 1.2*T_infty. 
	// If R > 0, then flame is spreading too fast and is spreading faster than frame -> increase V_f (flame moving towards inlet)
	// If R < 0, then flame is spreading too slow and frame is outpacing flame -> decrease V_f (flame moving towards outlet)
	real V_f_new = V_f_bias + Kp * R + I_term + Kd * dRdn;
	//real V_f_new = V_f + Kp * R + Kd * dRdn; 

	V_f = MAX(V_f_new, 0.0); //Floor V_f at 0 m/s
#endif

	//broadcast new V_f to nodes
	host_to_node_real_2(V_f, I_term);

#if !RP_NODE
	Message("Updated FSR to %g m/s using P-control with Kp = %g, Kd = %g, Ki = %g, and R = %g K, dRdn=%g, I_term=%g\n", V_f, Kp, Kd, Ki, R, dRdn, I_term);
#endif
}

// Function to execute scheme commands and update the grid
// When it is time to update, this function sets the string rp varibale that is the command that gets executed every iteration
// user define iterations to execute the grid motion 
// DEFINaE_ADJUST(update_grid, d)
// {
// 	// Check if it is time for update
// 	if (N_ITER < N_MESH_UPDATE_IGNORE)
// 	{
// 		Message0("Not enough iterations, skipping update. N_ITER=%d\n",N_ITER);
		
// 		return;
// 	}

// 	// Check to make sure flame spread it is greater than zero.
// 	if (V_f <= 0.0)
// 	{
// 		Message0("V_f too small. V_f = %g m/s\n",V_f);
// 		return;
// 	}

// 	RP_Set_String("update_grid_command",
//     "/define/user-defined/execute-on-demand \"calc_regression::lib_inlet_fsr\"\n"
//     "/solve/mesh-motion\n"
//     "yes\n");
// }

// Check variables that could be uninitialized at start of each iteration.
DEFINE_ADJUST(check_eigen_face_zoneID, d)
{
	if(eigen_face_zoneID == -1)
	{
		Message0("Error: eigen_face_zoneID = %d, STOPPING CALCULATION",eigen_face_zoneID);
		RP_Set_Integer("sol/iterations",0);

		Error("check_eigen_zone: eigen_face_zoneID not properly set. Please set the rp variable 'user/eigen_zone_id' and run the 'set_eigen_face_zone' ON_DEMAND UDF.\n");
	}

}

// This did not stop the calculation, why?
DEFINE_ADJUST(check_update_grid_command, d)
{
	bool scheme_command_exists = RP_Variable_Exists_P("user/update_grid_command");

	if (!scheme_command_exists)
	{
		Message0("Error: rp variable 'user/update_grid_command' does not exist. Cannot execute solve grid. STOPPING CALCULATION\n");
		RP_Set_Integer("sol/iterations",0);
	}
}

DEFINE_INIT(init_RP_vars, d)
{
#if !RP_NODE
	/* Check for V_f initializer value */
	bool V_f_init_exists =
		RP_Variable_Exists_P("user/v_f_init");

	Message(
		"Checking for user-defined parameter "
		"'user/v_f_init': %d\n",
		V_f_init_exists);

	if (V_f_init_exists)
	{
		V_f = RP_Get_Real("user/v_f_init");

		Message(
			"User-defined parameter 'user/v_f_init' "
			"found with value: %g m/s\n",
			V_f);
	}
	else
	{
		V_f = 0.0;

		Message(
			"Warning: User-defined parameter 'user/v_f_init' "
			"not found. Using default value of %g m/s.\n",
			V_f);
	}


	/* Check for proportional gain */
	bool Kp_exists =
		RP_Variable_Exists_P("user/kp");

	Message(
		"Checking for user-defined parameter "
		"'user/kp': %d\n",
		Kp_exists);

	if (Kp_exists)
	{
		Kp = RP_Get_Real("user/kp");

		Message(
			"User-defined parameter 'user/kp' "
			"found with value: %g\n",
			Kp);
	}
	else
	{
		Message(
			"Warning: User-defined parameter 'user/kp' "
			"not found. Using current value of %g.\n",
			Kp);
	}


	/* Check for integral gain */
	bool Ki_exists =
		RP_Variable_Exists_P("user/ki");

	Message(
		"Checking for user-defined parameter "
		"'user/ki': %d\n",
		Ki_exists);

	if (Ki_exists)
	{
		Ki = RP_Get_Real("user/ki");

		Message(
			"User-defined parameter 'user/ki' "
			"found with value: %g\n",
			Ki);
	}
	else
	{
		Message(
			"Warning: User-defined parameter 'user/ki' "
			"not found. Using current value of %g.\n",
			Ki);
	}


	/* Check for derivative gain */
	bool Kd_exists =
		RP_Variable_Exists_P("user/kd");

	Message(
		"Checking for user-defined parameter "
		"'user/kd': %d\n",
		Kd_exists);

	if (Kd_exists)
	{
		Kd = RP_Get_Real("user/kd");

		Message(
			"User-defined parameter 'user/kd' "
			"found with value: %g\n",
			Kd);
	}
	else
	{
		Message(
			"Warning: User-defined parameter 'user/kd' "
			"not found. Using current value of %g.\n",
			Kd);
	}


	/* Check for controller bias */
	bool V_f_bias_exists =
		RP_Variable_Exists_P("user/v_f_bias");

	Message(
		"Checking for user-defined parameter "
		"'user/v_f_bias': %d\n",
		V_f_bias_exists);

	if (V_f_bias_exists)
	{
		V_f_bias = RP_Get_Real("user/v_f_bias");

		Message(
			"User-defined parameter 'user/v_f_bias' "
			"found with value: %g m/s\n",
			V_f_bias);
	}
	else
	{
		Message(
			"Warning: User-defined parameter 'user/v_f_bias' "
			"not found. Using current value of %g m/s.\n",
			V_f_bias);
	}

	//I_term = V_f - V_f_bias;
	I_term = 0.0;

	Message(
		"PID controller initialized: "
		"V_f=%g m/s, V_f_bias=%g m/s, I_term=%g m/s, "
		"Kp=%g, Ki=%g, Kd=%g\n",
		V_f,
		V_f_bias,
		I_term,
		Kp,
		Ki,
		Kd);
#endif

	/*
	 * V_f is used by profile and motion UDFs on the compute nodes.
	 */
	host_to_node_real_1(V_f);
	Message0("V_f initialized to %g m/s\n", V_f);

	// broadcast to nodes
	host_to_node_real_1(Kp);
	Message0("Kp initialized to %g\n", Kp);

	host_to_node_real_1(Ki);
	Message0("Ki initialized to %g\n", Ki);

	host_to_node_real_1(Kd);
	Message0("Kd initialized to %g\n", Kd);

	host_to_node_real_1(V_f_bias);
	Message0("V_f_bias initialized to %g m/s\n", V_f_bias);

	host_to_node_real_1(I_term);
	Message0("I_term initialized to %g m/s\n", I_term);
}


// Update Solid Motion 
DEFINE_ZONE_MOTION(update_solid_motion, omega, axis, origin, velocity, current_time, dtime)
{
	NV_D(velocity, =, V_f, 0.0, 0.0); // Update solid motion velocity with calculated FSR
	//Message("Updated solid motion to %g m/s\n", V_f);
	return;
}

// Calculate new regressed surface coordinates (move to DEFINE_GRID_MOTION when debugged)
// This function calculates and populates the globals x_f_g and y_f_new_g for the DEFINE_GRID-MOTION
// UDFs to access. It also allocates them, so the variables need to be unallocated before this gets called.
// The process is:
//	1. calc_regression (populate x_f_g and y_f_new_g)
//	2. regress_surface for wall and its shadow (uses x_f_g and y_f_g to move nodes)
//	3. free x_f_g and y_f_new_g

// New process for this:
/*
1. Setup a calculation activity with these commands with the frequency you want the mesh to update at: 
	/define/user-defined/execute-on-demand/inlet_fsr:calc_regression
	(ti-menu-load-string (%rpgetvar 'user/update_grid_command)
2. calc_regression will always run and will calculate and populate the globals x_f_g and y_f_new_g for the DEFINE_GRID-MOTION
	UDFs to access. If the calculated profile is valid (doesn't contain nans/infinite slope) the command that gets executed in the activity
	'user/update_grid_command' is set to 'solve/mesh-motion yes' so that the grid is actually moved, otherwise the globals are populated
	with place holders and the 'user/update_grid_command' is set to a non-action, skipping the mesh update. 
3. When the profile is okay and the command is 'solve/mesh-motion yes', regress_surface is executed. This function checks again that 
	the profile is not null and moves the nodes accordingly. 
*/
DEFINE_ON_DEMAND(calc_regression)
{
	int profile_ok = 1;

#if !RP_HOST // Run on nodes in parallel or single process in serial host does nothing for this UDF
	// Data lookup variables
	Domain *d = Get_Domain(1); // Get domain pointer, update if different
	int regression_zone_id = 5; 
	Thread *t = Lookup_Thread(d, regression_zone_id); // get pointer to regressed surface thread

	real *y_f_new; // Create vectors to store face centroid coordinates.
	real *x_f, *y_f, *mdot_f, *A_f, *dx_f; //Current centroid coordinates vector
	int *idx_f; // Array to store indicies for f 

	real r[ND_ND]; // vector to store face centroid coordinates
	//FaceData *faces_data;

	face_t f;

	Thread *t_eigen = Lookup_Thread(d, eigen_face_zoneID); // get pointer to eigen face thread
	real r_eig[ND_ND] = {0.0, 0.0}; // coordinate of eigen face centroid, 
	real iwork[ND_ND];  // i work is temp array storage
	

	// Node **eigen_nodes; // create vector of node pointers for eigen face nodes
	// Node *v; // node pointer for looping throuh nodes

	// Integration variables
	const real rho = 1190; // [kg/m^3] Density of solid phase

	// Integrals in the numerator and deniminator
	real I_num = 0;
	real I_denom = 0;
	real I = 0; // Integral of x-component of face area vector, sum of A_face_x = sum of A_face * f'(x)/sqrt(f'(x)^2+1) over all faces along surface

	real A_face[ND_ND]; // Face area vector
	real A_face_mag; // Face area magnitude

	// Get the eigen face nodes from the eigen face zone ID
	begin_f_loop(f, t_eigen)
	{
		if (PRINCIPAL_FACE_P(f, t_eigen))
		{
			F_CENTROID(r_eig, f, t_eigen); // get eigen face centroid coordinates
		}
	}
	end_f_loop(f, t_eigen)

	// Synchronize centroid coordinates
	PRF_GRSUM(r_eig, ND_ND, iwork);
	r_eig_g[0] = r_eig[0];
	r_eig_g[1] = r_eig[1];

	Message0("Eigenface centroid found. x_eig = %g m, y_eig = %g m \n", r_eig[0], r_eig[1]);

	// Get number of interior faces on thread on each compute node
	int size = 0;
	begin_f_loop(f, t)
	{
		if PRINCIPAL_FACE_P(f, t)
		{
			size++;
		}
	}
	end_f_loop(f, t)
	
	Message("faces:%d found on partition/node: %d \n", size, myid);

	// Send the size (number of faces) to node zero from each node
	if (! I_AM_NODE_ZERO_P)
	{
		PRF_CSEND_INT(node_zero, &size, 1, myid);
	}

	// Allocate memory for faces_data and y_f_new on each compute node
	x_f = malloc(size * sizeof(real));
	y_f = malloc(size * sizeof(real));
	mdot_f = malloc(size * sizeof(real));
	A_f = malloc(size * sizeof(real));
	dx_f = malloc(size * sizeof(real));
	// idx_f = malloc(size * sizeof(int));

	//faces_data = malloc(size * sizeof(FaceData));
	//y_f_new = malloc(size * sizeof(real));
	Message0("faces_data allocated on nodes\n");

	// Fill arrays on each node
	int i = 0;
	begin_f_loop(f, t)
	{
		if (PRINCIPAL_FACE_P(f, t))
		{
			//faces_data[i].f = f; 
		
			// faces[i] = f; 
			F_CENTROID(r, f, t); //get face centroid coordinates

			x_f[i] = r[0];
			y_f[i] = r[1];
			mdot_f[i] = F_FLUX(f, t);

			F_AREA(A_face, f, t); // Get face area vector for current face
			A_face_mag = NV_MAG(A_face);
			A_f[i] = A_face_mag; // Get face area magnitude for current face
			dx_f[i] = fabs(A_face[1]); // get dx which is y-normal component of face area vector
			//nhat_x[i] = A_face[0] / A_face_mag; // Get x-component of face normal vector for current face

			//idx_f[i] = i;
			//faces_data[i].x_f = r[0];
			//faces_data[i].y_f = r[1];

			i++;
		}
	}
	end_f_loop(f, t)

	// Send face coordinates array to node zero from node 1,2,...
	Message0("Starting Transfer from node 1,2,... to 0\n");
	if (! I_AM_NODE_ZERO_P)
	{
		PRF_CSEND_REAL(node_zero, x_f, size, myid);
		PRF_CSEND_REAL(node_zero, y_f, size, myid);
		PRF_CSEND_REAL(node_zero, mdot_f, size, myid);
		PRF_CSEND_REAL(node_zero, A_f, size, myid);
		PRF_CSEND_REAL(node_zero, dx_f, size, myid);
		
		//PRF_CSEND_INT(node_zero, idx_f, size, myid);
	}

	// Recieve data on node zero
	if (I_AM_NODE_ZERO_P)
	{
		// I think loop over all nodes except node zero (This is not well documented by fluent)
		// Supposedly this macro is in para.h, but it does not seem to be there anymore.
		compute_node_loop_not_zero(i)
		{
			int old_size = size;

			PRF_CRECV_INT(i, &size, 1, i); //Replace size with the value sent from node i earlier

			// Calculate new size to append data to single node zero owned array 
			int new_size = old_size + size; 

			Message0("Reallocating arrays from size:%d to new_size: %d \n", old_size, new_size);

			// Reallocate arrays
			x_f = realloc(x_f, new_size * sizeof(real));
			y_f = realloc(y_f, new_size * sizeof(real));
			mdot_f = realloc(mdot_f, new_size * sizeof(real));
			A_f = realloc(A_f, new_size * sizeof(real));
			dx_f = realloc(dx_f, new_size * sizeof(real));

			//idx_f = realloc(idx_f, new_size * sizeof(int));

			// Recieve data from arrays sent earlier from node 1,2,..
			// Store the data in the appended room to the x_f and y_f arrays
			PRF_CRECV_REAL(i, &x_f[old_size], size, i);
			PRF_CRECV_REAL(i, &y_f[old_size], size, i);
			PRF_CRECV_REAL(i, &mdot_f[old_size], size, i);
			PRF_CRECV_REAL(i, &A_f[old_size], size, i);
			PRF_CRECV_REAL(i, &dx_f[old_size], size, i);
			
			//PRF_CRECV_INT(i, &idx_f[old_size], size, i);

			// update size
			size = new_size; 

			Message0("Data recieved from node %d \n", i);
		}

		// Sort using qsort_s (may need to be switched to qsort_r when using linux)
		// Sort index array based on context array x_f

		// build index array
		idx_f = malloc(size * sizeof(int));
		for (int l = 0; l < size; l++) {idx_f[l] = l;}

		// Use correct sorting for platform
		#ifdef _WIN32
		qsort_s(idx_f, size, sizeof(int), compare_xf, x_f);
		#elif defined(__linux__)
		qsort_r(idx_f, size, sizeof(int), compare_xf_linux, x_f);
		#else
		#error Unsupported platform.
		#endif

		// Re order x_f and y_f
		real *temp_x_f = malloc(size * sizeof(real));
		real *temp_y_f = malloc(size * sizeof(real));
		real *temp_mdot_f =  malloc(size * sizeof(real));
		real *temp_A_f = malloc(size * sizeof(real));
		real *temp_dx_f = malloc(size * sizeof(real));

		memcpy(temp_x_f, x_f, size * sizeof(real));
		memcpy(temp_y_f, y_f, size * sizeof(real));
		memcpy(temp_mdot_f, mdot_f, size * sizeof(real));
		memcpy(temp_A_f, A_f, size * sizeof(real));
		memcpy(temp_dx_f, dx_f, size * sizeof(real));
	
		for (int j = 0; j < size; j++)
		{
			x_f[j] = temp_x_f[idx_f[j]];
			y_f[j] = temp_y_f[idx_f[j]];
			mdot_f[j] = temp_mdot_f[idx_f[j]];
			A_f[j] = temp_A_f[idx_f[j]];
			dx_f[j] = temp_dx_f[idx_f[j]];
		}
		
		free(temp_x_f); //free temporary arrays
		free(temp_y_f); 
		free(temp_mdot_f);
		free(temp_A_f);
		free(temp_dx_f);
		free(idx_f);

		// Perform cumulative integral
		// checks for infinite slope
		
		// allocate y_f_new on node zero 
		y_f_new = malloc(size * sizeof(real));
		
		for (int k = 0; k < size; k++)
		{
			// check for infinite slope
			real rad = pow(rho * V_f * A_f[k], 2) - pow(mdot_f[k], 2);

			if (rad <= 0)
			{
				profile_ok = 0;
				Message0("Warning: Infinite slope detected, regression profile invalid. mdot_f = %g kg/s, rhoV_fA = %g kg/s, x = %g m\n",
					 mdot_f[k], rho * V_f * A_f[k], x_f[k]);
				break;
			}

			I += mdot_f[k] * dx_f[k] / sqrt(pow(rho * V_f * A_f[k], 2) - pow(mdot_f[k], 2));
			//I += mdot_f[k]*nhat_x[k] / A_f[k] / sqrt(pow(rho * V_f * A_f[k], 2) + pow(mdot_f[k], 2));

			y_f_new[k] = r_eig[1] + I;

			// // Apply under-relaxation to new profile to help steady state convergence
			// y_f_new[k] = y_f[k] + alpha_grid * (y_f_new[k] - y_f[k]);
			
			// // Update the largest displacement
			// real displacement = fabs(y_f_new[k] - y_f[k]);
			// if (displacement > largest_displacement)
			// {
			// 	largest_displacement = displacement;
			// }

			Message0("y_f_new(x = %g m) = %g m \n", x_f[k], y_f_new[k]);
		}	

		/* If any point failed, overwrite the partially computed profile with a finite dummy profile.
		It should not be used because mycommand will be set to no-op, but this avoids sending NaNs. */
		if (!profile_ok)
		{
			for (int k = 0; k < size; k++)
			{
				y_f_new[k] = y_f[k];
			}

			Message0("Warning: Regression rejected. profile_ok = %d \n", profile_ok);
		}

		// Send a copy of the sorted x_f and y_f_new to all other nodes
		compute_node_loop_not_zero(i)
		{
			Message0("Sending size of y_f_new from node 0 to %d \n", i);
			PRF_CSEND_INT(i, &size, 1, myid); //send size of y_f new to other nodes 
		}

	}

	// Allocate memory for y_f_new on all other nodes
	// reallocate x_f since it should already exist on other nodes from initial array fill
	if (! I_AM_NODE_ZERO_P)
	{
		PRF_CRECV_INT(node_zero, &size, 1, node_zero);
		y_f_new = malloc(size * sizeof(real));
		x_f = realloc(x_f, size * sizeof(real));
	}
	Message0("y_f_new allocated and x_f reallocated on all nodes\n");

	// Wait and make sure y_f_new is allocated
	PRF_GSYNC();

	// Send y_f_new array from zero to 1,2...
	if (I_AM_NODE_ZERO_P)
	{
		compute_node_loop_not_zero(i)
		{
			PRF_CSEND_REAL(i, y_f_new, size, myid);
			PRF_CSEND_REAL(i, x_f, size, myid);
		}
	}

	// Recieve y_f_new and x_f from node zero
	if (! I_AM_NODE_ZERO_P)
	{
		PRF_CRECV_REAL(node_zero, y_f_new, size, node_zero);
		PRF_CRECV_REAL(node_zero, x_f, size, node_zero);
		Message("y_f_new and x_f recieved on all other node %d from node 0 \n", myid);
	}

	// Set global storage on compute nodes. If fail (install_profile returns 1), cleanup
	if (install_profile(&x_f, &y_f_new, (size_t)size))
	{
		free(x_f);
		free(y_f);
		free(mdot_f);
		free(A_f);
		free(dx_f);
		free(y_f_new);

		Error("calc_regression: install_profile failed.\n");
	}

	// x_f_g = malloc(size * sizeof(real));
	// y_f_new_g = malloc(size * sizeof(real));
	// size_g = size; 

	// memcpy(x_f_g, x_f, size * sizeof(real));
	// memcpy(y_f_new_g, y_f_new, size * sizeof(real));

	// Loop over nodes and interpolate their new position
	Node *v; 
	int n;

	profile_ok = PRF_GILOW1(profile_ok); //if profile_ok is zero on any node, set it to zero on all nodes
	if (profile_ok)
	{
		largest_displacement = 0.0;

		begin_f_loop(f, t)
		{
			f_node_loop(f, t, n)
			{
				v = F_NODE(f, t, n);
				real x_node = NODE_X(v);
				real y_node_new;
				real y_node_old = NODE_Y(v);
				real displacement_residual = 0.0;


				// if the node is left of the first centroid, it connects to the
				// eigen face and should not move. 
				if (x_node < x_f_g[0])
				{
					y_node_new = r_eig[1]; 
				}
				else
				{
					y_node_new = interp1d(x_f_g, y_f_new_g, size, x_node, true);
				}
				
				// Calcuate displacement residual based on target position and current position. 
				// Do not use the under-relaxes postion b/c that would change the residual with the under-relaxation
				displacement_residual = fabs(y_node_new - y_node_old);
				if (displacement_residual > largest_displacement)
				{
					largest_displacement = displacement_residual;
				}

				// Apply under-relaxation to new profile
				y_node_new = y_node_old + alpha_grid * (y_node_new - y_node_old);

				//displacement_residual = fabs(y_node_new - y_node_old);

				
				// Calculate 
				// Print new coordinates (may duplicate)
				Message("Node would be moved: y_node_new(x = %g m) = %g m \n", x_node, y_node_new);
			}
		}
		end_f_loop(f, t)

	}

	// Free allocated memory on all nodes. Do not free x_f or y_f_new after installing the profile
	//free(x_f);
	free(y_f);
	free(mdot_f);
	free(A_f);
	free(dx_f);
	//free(y_f_new); 
	
	//get the largest displacement on all nodes
	largest_displacement = PRF_GRHIGH1(largest_displacement); //get the largest displacement on all nodes
	Message0("calc_regression: largest_displacement = %g m \n", largest_displacement);
	
#endif 

	//update profile_ok on host
	node_to_host_int_1(profile_ok);

	// update largest_displacement on host
	node_to_host_real_1(largest_displacement);

	//PRF_GSYNC();
	// Set the proper mesh action command
#if !RP_NODE
	if (profile_ok && fabs(R) <= 0.001)
	{
		RP_Set_String(
			"user/update_grid_command",
			"/solve/mesh-motion yes\n"
			"/define/models/radiation/s2s-parameters/compute-write-vf \"SYS-4.s2s.h5\" yes\n");

		Message(
			"calc_regression: valid profile and |R| = %g <= 0.001 K. "
			"Mesh motion enabled.\n",
			fabs(R));

		N_MESH_UPDATES++;
		Message("N_MESH_UPDATES=%d\n", N_MESH_UPDATES);
	}
	else
	{
		RP_Set_String("user/update_grid_command", "\n");

		if (!profile_ok)
		{
			Message(
				"Warning: calc_regression: invalid profile. "
				"Mesh motion skipped.\n");
		}
		else
		{
			Message(
				"calc_regression: |R| = %g > 0.001 K. "
				"Mesh motion skipped.\n",
				fabs(R));
		}
	}
#endif
}

DEFINE_GRID_MOTION(regress_surface, d, dt, time, dtime)
{
#if !RP_HOST
	//Check for valid profile
	if (!valid_profile || x_f_g == NULL || y_f_new_g == NULL || size_g < 2)
	{
		Error("regress_surface: no valid regression profile. Run calc_regression first.\n");
	}

	// Convert dynamic thread into normal thread
	Thread *t = DT_THREAD(dt);
	face_t f;

	Node *v; 
	int n;
	real r_new[ND_ND];

	// Set deforming flag on adjacent cell zone/s?
	SET_DEFORMING_THREAD_FLAG(THREAD_T0(t));

	begin_f_loop(f, t)
	{
		f_node_loop(f, t, n)
		{
			v = F_NODE(f, t, n);
			real x_node = NODE_X(v);
			real y_node_old = NODE_Y(v);
			real y_node_new;

			if (NODE_POS_NEED_UPDATE (v))
			{
				NODE_POS_UPDATED(v); //indicate that node position was updated to prevent repeats
				// if the node is left of the first centroid, it connects to the
				// eigen face and should not move. 
				if (x_node < x_f_g[0])
				{
					y_node_new = r_eig_g[1]; 
				}
				else
				{
					y_node_new = interp1d(x_f_g, y_f_new_g, size_g, x_node, true);
				}

				// Apply under-relaxation to new profile
				y_node_new = y_node_old + alpha_grid * (y_node_new - y_node_old);
				
				r_new[0] = x_node;
				r_new[1] = y_node_new;

				// Set the new node coordinates
				NV_V(NODE_COORD(v), = , r_new);

				// Print new coordinates (may duplicate)
				Message("Node moved: y_node_new(x = %g m) = %g m \n", x_node, y_node_new);
			}			
		}
	}
	end_f_loop(f, t)

#endif
}

// Only run after calc_regression has been run
DEFINE_ON_DEMAND(free_xf_and_yf_new_globs)
{
#if !RP_HOST
	// free(x_f_g);
	// free(y_f_new_g);
	cleanup_profile();
	Message0("x_f_g and y_f_g freed \n");
#endif
}

DEFINE_GRID_MOTION(restore_exp_surface, d, dt, time, dtime)
{
	Thread *t = DT_THREAD(dt);
	face_t f;
	Node *v;
	int n;

	// Set deforming flag on the adjacent cell zone/s?
	SET_DEFORMING_THREAD_FLAG(THREAD_T0(t));

	real r_new[ND_ND];

	begin_f_loop(f, t)
	{
		f_node_loop(f, t, n)
		{
			v = F_NODE(f, t, n);

			if (NODE_POS_NEED_UPDATE (v))
			{
				NODE_POS_UPDATED(v);
				real x = NODE_X(v);
				real y;

				if (x < 0.020)
				{
					y = 0.0;
				}
				else
				{
					y = -0.52 * pow(x*1000 - 20 ,0.64)/1000;
				}

				r_new[0] = x;
				r_new[1] = y;

				NV_V(NODE_COORD(v), =, r_new);
			}
		}
	}
	end_f_loop(f, t)
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

//Limit massflux on surface to rho*V_f
// Outputs:
//		net molar reaction rates rr (in kmol/m^3-s for volumetric reaction and kmol/m^2-s for surface reaction)
//		the units of rr depend on the cell the UDF is being called from.
//		jacobian, the derivative of the surface net reaction rate wrt to species concentration		   
// This does not override the stiff chemistry solver, it just replaces the reaction rate the the 
// stiff chemistry solver integrates. The species are just evolved with a limited reaction rate.
// DEFINE_NET_REACTION_RATE(limit_surf, c, t, particle, pressure, temp, yi, rr, jac)
// {
// 	//Check for surface or volumetric cell
// 	if (BOUNDARY_FACE_THREAD_P(t))
// 	{
// 		// Surface call rr[kmol / (m^2-s)]
// 		// Pyrolysis only.
// 	}
// 	else
// 	{
// 		// Volumetric call rr[kmol / (m^3-s)]
// 	}
// 	RP_Get_List_Ref_Real()
// }

// Create Report Definition for FSR
DEFINE_REPORT_DEFINITION_FN(flame_spread_rate)
{
	return V_f; // Return calculated flame spread rate for report definition
}

// Report Definition for Residual of eigenvalue-based FSR calculation
DEFINE_REPORT_DEFINITION_FN(R_eigen_temp_diff)
{
	return R; // Return residual for report definition
}

// Stores the largest displacment of the last/most recent mesh update
DEFINE_REPORT_DEFINITION_FN(largest_displacement_report)
{
	return largest_displacement; // Return largest displacement for report definition
}

// Report definition to make a user defined absolute convergence condition.
// Returns the iteration number when largest_displacemt > threshold
// returns the iteration number when largest_displacement crossed the threshold when largest_displacement < threshold
DEFINE_REPORT_DEFINITION_FN(largest_displacement_conv_report)
{
	// Since largest_displacement = 0 at start, do not allow convergence during the initial 
	// surface calculation
	if (fabs(largest_displacement) > DISPLACEMENT_THRESHOLD || N_MESH_UPDATES < 1)
	{
		largest_displacement_conv = N_ITER;
	}

	return largest_displacement_conv;
}

DEFINE_ON_DEMAND(check_rp_vars)
{
	bool U_mean_exists = RP_Variable_Exists_P("user/u_mean"); // Check if user-defined parameter for mean velocity exists
	//bool alpha_exists = RP_Variable_Exists_P("user/alpha"); // Check if user-defined parameter for under-relaxation factor exists

	Message0("Checking for user-defined parameter 'user/u_mean': %d\n", U_mean_exists);
	//Message0("Checking for user-defined parameter 'user/alpha': %d\n", alpha_exists);

	if (U_mean_exists)
	{
		real U_mean = RP_Get_Real("user/u_mean"); // Get mean velocity from user-defined parameter if it exists
		Message0("User-defined parameter 'user/u_mean' found with value: %f m/s\n", U_mean);
	}
	else
	{
		Message0("Warning: User-defined parameter 'user/u_mean' not found.\n");
	}

	// if (alpha_exists)
	// {
	// 	alpha = RP_Get_Real("user/alpha"); // Get under-relaxation factor from user-defined parameter if it exists
	// 	Message0("User-defined parameter 'user/alpha' found with value: %f\n", alpha);
	// }
	// else
	// {
	// 	Message0("Warning: User-defined parameter 'user/alpha' not found. Using default value: %f\n", alpha);
	// }

}

DEFINE_ON_DEMAND(set_FSR)
{

#if !RP_NODE

	bool V_f_init_exists = RP_Variable_Exists_P("user/v_f_init"); // Check if user-defined parameter for initial FSR exists)

	Message("Checking for user-defined parameter 'user/v_f_init': %d\n", V_f_init_exists);

	if (V_f_init_exists)
	{
		V_f = RP_Get_Real("user/v_f_init"); // Get initial FSR from user-defined parameter if it exists
		Message("User-defined parameter 'user/v_f_init' found with value: %g m/s\n", V_f);
	}
	else
	{
		V_f = 0.0; // Default initial FSR value if user-defined parameter does not exist
		Message("Warning: User-defined parameter 'user/v_f_init' not found. Using default value of 0 m/s.\n");
	}

	// Set the I term when V_f is changed so there is no ''bump''
	I_term = V_f - V_f_bias;
#endif

	//node_to_host_real_1(V_f); // update V_f on host process so report is correct.
	host_to_node_real_2(V_f, I_term);
	Message("V_f initialized to %g m/s \n",V_f);
}

DEFINE_ON_DEMAND(set_Kp)
{
#if !RP_NODE
	bool Kp_exists = RP_Variable_Exists_P("user/kp");

	if (Kp_exists)
	{
		Kp = RP_Get_Real("user/kp");
		Message("User-defined parameter 'user/kp' found with value: %g\n", Kp);
	}
	else
	{
		Kp = 1e-6; // Default initial FSR value if user-defined parameter does not exist
		Message("Warning: User-defined parameter 'user/kp' not found. Using default value of 1e-6.\n");
	}

#endif

	// Pass parameters to nodes
	host_to_node_real_1(Kp);
	Message0("Kp initialized to %g \n",Kp);
}

DEFINE_ON_DEMAND(set_Kd)
{
#if !RP_NODE
	bool Kd_exists = RP_Variable_Exists_P("user/kd");

	if (Kd_exists)
	{
		Kd = RP_Get_Real("user/kd");
		Message("User-defined parameter 'user/kd' found with value: %g\n", Kd);
	}
	else
	{
		Kd = 1e-7; // Default initial FSR value if user-defined parameter does not exist
		Message("Warning: User-defined parameter 'user/kd' not found. Using default value of 1e-7.\n");
	}

#endif

	// Pass parameters to nodes
	host_to_node_real_1(Kd);
	Message0("Kd initialized to %g \n",Kd);
}

DEFINE_ON_DEMAND(set_Ki)
{
#if !RP_NODE
	bool Ki_exists = RP_Variable_Exists_P("user/ki");

	if (Ki_exists)
	{
		Ki = RP_Get_Real("user/ki");
		Message("User-defined parameter 'user/ki' found with value: %g\n", Ki);
	}
	else
	{
		Ki = 7.5e-9; // Default initial FSR value if user-defined parameter does not exist
		Message("Warning: User-defined parameter 'user/ki' not found. Using default value of %g\n", Ki);
	}

#endif

	// Pass parameters to nodes
	host_to_node_real_1(Ki);
	Message0("Ki initialized to %g \n",Ki);
}

DEFINE_ON_DEMAND(set_V_f_bias)
{
#if !RP_NODE
	bool V_f_bias_exists =
		RP_Variable_Exists_P("user/v_f_bias");

	if (V_f_bias_exists)
	{
		real V_f_bias_new =
			RP_Get_Real("user/v_f_bias");

		// Preserve V_f_bias + I_term so changing the bias
		// does not cause a sudden controller-output jump.
		 
		I_term += V_f_bias - V_f_bias_new;
		V_f_bias = V_f_bias_new;

		Message(
			"User-defined parameter 'user/v_f_bias' "
			"found with value: %g m/s\n",
			V_f_bias);
	}
	else
	{
		real V_f_bias_new = 40e-6;

		
		// Apply the default bias without changing the current
		//zero-error controller baseline.
		
		I_term += V_f_bias - V_f_bias_new;
		V_f_bias = V_f_bias_new;

		Message(
			"Warning: User-defined parameter 'user/v_f_bias' "
			"not found. Using default value of %g m/s\n",
			V_f_bias);
	}
#endif

	//Pass parameters to nodes 
	host_to_node_real_2(V_f_bias, I_term);

	Message0(
		"V_f_bias initialized to %g m/s; "
		"I_term adjusted to %g m/s\n",
		V_f_bias,
		I_term);
}

DEFINE_ON_DEMAND(set_T_eig_setpoint)
{
#if !RP_NODE
	bool t_eig_setpoint_exists = RP_Variable_Exists_P("user/t_eig_setpoint");

	if (t_eig_setpoint_exists)
	{
		T_eig_setpoint = RP_Get_Real("user/t_eig_setpoint");
		Message("User-defined parameter 'user/t_eig_setpoint' found with value %g\n", T_eig_setpoint);
	}
	else
	{
		T_eig_setpoint = 360; 
		Message("Warning: User-defined parameter 'user/t_eig_setpoint' not found. Using default value of %g.\n", T_eig_setpoint);
	}
#endif

	// Pass to nodes 
	host_to_node_real_1(T_eig_setpoint);
	Message0("T_eig_setpoint initialized to %g \n", T_eig_setpoint);
}

DEFINE_ON_DEMAND(set_N_MESH_UPDATE_IGNORE)
{
	#if !RP_NODE
	bool n_mesh_ignore_exists = RP_Variable_Exists_P("user/n_mesh_ignore");

	if (n_mesh_ignore_exists)
	{
		N_MESH_UPDATE_IGNORE = RP_Get_Integer("user/n_mesh_ignore");
		Message("User-defined parameter 'user/n_mesh_ignore' found with value: %d\n", N_MESH_UPDATE_IGNORE);
	}
	else
	{
		N_MESH_UPDATE_IGNORE = 5000;
		Message("Warning: User-defined parameter 'user/n_mesh_ignore' not found. Using default value of %d.\n", N_MESH_UPDATE_IGNORE);
	}

#endif

	// Pass parameters to nodes
	host_to_node_int_1(N_MESH_UPDATE_IGNORE);
	Message0("N_MESH_UPDATE_IGNORE initialized to %d \n",N_MESH_UPDATE_IGNORE);
}

// DEFINE_ON_DEMAND(set_UPDATE_INTERVAL)
// {
// #if !RP_NODE
// 	// Check for UPDATE INTERVAL value
// 	bool UPDATE_INTERVAL_exists = RP_Variable_Exists_P("user/update_interval");

// 	Message("Checking for user-defined parameter 'user/update_interval': %d\n", UPDATE_INTERVAL_exists);

// 	if (UPDATE_INTERVAL_exists)
// 	{
// 		UPDATE_INTERVAL = RP_Get_Integer("user/update_interval");
// 		Message("User-defined parameter 'user/update_interval' found with value: %d\n", UPDATE_INTERVAL);
// 	}
// 	else
// 	{
// 		UPDATE_INTERVAL = 100; // Default initial FSR value if user-defined parameter does not exist
// 		Message("Warning: User-defined parameter 'user/update_interval' not found. Using default value of 100.\n");
// 	}
// #endif

// 	// Pass parameters to nodes
// 	host_to_node_int_1(UPDATE_INTERVAL);
// 	Message0("UPDATE_INTERVAL initialized to %d iterations \n",UPDATE_INTERVAL);
// }

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
#if !RP_NODE
	bool zone_ID_exists = RP_Variable_Exists_P("user/eigen_zone_id"); // Check if user-defined parameter for eigen face zone ID exists
	Message0("Checking for user-defined parameter 'user/eigen_zone_id': %d\n", zone_ID_exists);
	if (zone_ID_exists)
	{
		eigen_face_zoneID = RP_Get_Integer("user/eigen_zone_id"); // Get eigen face zone ID from user-defined parameter if it exists
		Message("User-defined parameter 'user/eigen_zone_id' found with value: %d\n", eigen_face_zoneID);
	}
	else
	{
		eigen_face_zoneID = -1; // Default value if user-defined parameter does not exist, update with different default if desired
		Message("Warning: User-defined parameter 'user/eigen_zone_id' not found. Using default value of -1. Please set the rp variable. \n");
	}
#endif
	//node_to_host_int_1(eigen_face_zoneID); // update eigen_face_zoneID on host process so it can be used in calc_FSR_eigen
	
	host_to_node_int_1(eigen_face_zoneID); //broadcast to nodes
	Message("eigen_face_zoneID initialized to %d \n", eigen_face_zoneID);
}

DEFINE_ON_DEMAND(set_grid_urf)
{

#if !RP_NODE
	bool grid_urf_exists = RP_Variable_Exists_P("user/grid_urf");

	if (grid_urf_exists)
	{
		alpha_grid = RP_Get_Real("user/grid_urf");
		Message("User-defined parameter 'user/grid_urf' found with value: %g\n", alpha_grid);
	}
	else
	{
		alpha_grid = 0.5; // Default value if user-defined parameter does not exist
		Message("Warning: User-defined parameter 'user/grid_urf' not found. Using default value of %g.\n", alpha_grid);
	}
#endif

	// Pass parameters to nodes
	host_to_node_real_1(alpha_grid);
	Message0("alpha_grid initialized to %g \n",alpha_grid);
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
	host_to_node_real_1(V_f);
}

DEFINE_RW_HDF_FILE(write_FSR_hdf, filename)
{
	// Saves data to cas and dat file. 
	size_t nelems = 1;

	char *fsr_path = "/FSR_data";
	Write_Complete_User_Dataset(
		filename, fsr_path, &V_f, nelems);

	char *v_f_bias_path = "/V_f_bias_data";
	Write_Complete_User_Dataset(
		filename, v_f_bias_path, &V_f_bias, nelems);

	char *kp_path = "/Kp_data";
	Write_Complete_User_Dataset(
		filename, kp_path, &Kp, nelems);

	char *ki_path = "/Ki_data";
	Write_Complete_User_Dataset(
		filename, ki_path, &Ki, nelems);

	char *kd_path = "/Kd_data";
	Write_Complete_User_Dataset(
		filename, kd_path, &Kd, nelems);

	char *eigen_zone_path = "/eigen_zoneID_data";
	real eigen_face_zoneID_real = (real)eigen_face_zoneID;

	Write_Complete_User_Dataset(
		filename,
		eigen_zone_path,
		&eigen_face_zoneID_real,
		nelems);

	char *I_term_path = "/I_term";
	Write_Complete_User_Dataset(filename, I_term_path, &I_term, nelems);

	char* eigen_temp_setpoint_path = "/t_eig_setpoint";
	Write_Complete_User_Dataset(
		filename, eigen_temp_setpoint_path, &T_eig_setpoint, nelems);

	char* alpha_grid_path = "/alpha_grid";
	Write_Complete_User_Dataset(
		filename, alpha_grid_path, &alpha_grid, nelems);

	char* largest_displacement_path = "/largest_displacement";
	Write_Complete_User_Dataset(
		filename, largest_displacement_path, &largest_displacement, nelems);

	Message0(
		"FSR parameters written: "
		"V_f=%g m/s, V_f_bias=%g m/s, "
		"Kp=%g, Ki=%g, Kd=%g, eigen_zoneID=%d, I_term=%g, T_eig_setpoint=%g, alpha_grid=%g, largest_displacement=%g \n",
		V_f,
		V_f_bias,
		Kp,
		Ki,
		Kd,
		eigen_face_zoneID,
		I_term,
		T_eig_setpoint,
		alpha_grid,
		largest_displacement);
}

DEFINE_RW_HDF_FILE(read_FSR_hdf, filename)
{
	size_t nelems = 1;

	char *fsr_path = "/FSR_data";
	Read_Complete_User_Dataset(
		filename, fsr_path, &V_f, nelems);

	char *v_f_bias_path = "/V_f_bias_data";
	Read_Complete_User_Dataset(
		filename, v_f_bias_path, &V_f_bias, nelems);

	char *kp_path = "/Kp_data";
	Read_Complete_User_Dataset(
		filename, kp_path, &Kp, nelems);

	char *ki_path = "/Ki_data";
	Read_Complete_User_Dataset(
		filename, ki_path, &Ki, nelems);

	char *kd_path = "/Kd_data";
	Read_Complete_User_Dataset(
		filename, kd_path, &Kd, nelems);

	char *eigen_zone_path = "/eigen_zoneID_data";
	real eigen_face_zoneID_real = 0.0;
	Read_Complete_User_Dataset(
		filename,
		eigen_zone_path,
		&eigen_face_zoneID_real,
		nelems);

	eigen_face_zoneID = (int)eigen_face_zoneID_real;

	char *I_term_path = "/I_term";
	Read_Complete_User_Dataset(filename, I_term_path, &I_term, nelems);

	char *eigen_temp_setpoint_path = "/t_eig_setpoint";
	Read_Complete_User_Dataset(
		filename, eigen_temp_setpoint_path, &T_eig_setpoint, nelems);

	char* alpha_grid_path = "/alpha_grid";
	Read_Complete_User_Dataset(
		filename, alpha_grid_path, &alpha_grid, nelems);

	char* largest_displacement_path = "/largest_displacement";
	Read_Complete_User_Dataset(
		filename, largest_displacement_path, &largest_displacement, nelems);

	/*
	 * Reset the controller state consistently with the restored
	 * flame-spread rate and bias.
	 */
	//I_term = V_f - V_f_bias;

	R = 0.0;
	R_old = 0.0;
	dRdn = 0.0;
	dRdn_old = 0.0;

	/* Send values used on compute nodes. */
	// Limit of values for this command is somewhere between 6 and 9
	node_to_host_real_6(V_f, V_f_bias, Kp, Ki, Kd, I_term);
	node_to_host_real_3(T_eig_setpoint, alpha_grid, largest_displacement);
	// host_to_node_real_1(V_f);
	// host_to_node_real_1(V_f_bias);
	// host_to_node_real_1(Kp);
	// host_to_node_real_1(Ki);
	// host_to_node_real_1(Kd);

	node_to_host_int_1(eigen_face_zoneID);
	//host_to_node_int_1(eigen_face_zoneID);

	Message0(
		"FSR parameters read: "
		"V_f=%g m/s, V_f_bias=%g m/s, "
		"Kp=%g, Ki=%g, Kd=%g, eigen_zoneID=%d, "
		"I_term=%g, T_eig_setpoint=%g, alpha_grid=%g, largest_displacement=%g\n",
		V_f,
		V_f_bias,
		Kp,
		Ki,
		Kd,
		eigen_face_zoneID,
		I_term,
		T_eig_setpoint,
		alpha_grid,
		largest_displacement);
}

