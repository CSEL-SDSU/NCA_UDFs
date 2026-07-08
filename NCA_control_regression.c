#include <stdlib.h>
#include "udf.h"
#include <stdbool.h>
#include "hdfio.h"
#include <math.h>

#define EPS 2.2204460492503131e-16
#define sgn(x)  ((x>0) - (x<0))

static real V_f; //Flame Spread Rate 
static int UPDATE_INTERVAL = 100; // Number of iterations between FSR updates
static real Kp = 5e-8; // Proportional gain for P-control update of FSR

//Variables required for evaluating F(X) - The residual temperature
// And determining when an evaluation is accepted
static real R;
static int eigen_face_zoneID = 0; // Zone ID of separated surface where temp is monitored
static const real T_infty = 300; 

typedef struct face_data_struct {
		face_t f;
		real x_f;
		real y_f;
	} FaceData;

// Comparator function for qsort. 
// Should return:
// < 0 (negative) if a should be before b
// 0 if a and b are equal 
// > 0 (positive) if a should be after b)
int compare_xf(void *context, const void* a, const void* b)
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

// Function to clear the storage of R

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
	R = T_eig - 1.2 * T_infty;
#endif

	// update R, T_eig, and face count on host process so report is correct
	node_to_host_real_2(R,T_eig);    
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
		Message("Iteration %d: R = %g K, T_eig = %g K, V_f = %g m/s\n", N_ITER, R, T_eig, V_f);
	}
#endif
}

DEFINE_EXECUTE_AT_END(update_FSR_Pcontrol)
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

	// Update FSR using P-control based on the residual temperature R
	// R = T(x_eig) - 1.2*T_infty. 
	// If R > 0, then flame is spreading too fast and is spreading faster than frame -> increase V_f (flame moving towards inlet)
	// If R < 0, then flame is spreading too slow and frame is outpacing flame -> decrease V_f (flame moving towards outlet)
	real V_f_new = V_f + Kp * R; 

	V_f = MAX(V_f_new, 0.0); //Floor V_f at 0 m/s
#endif

	//broadcast new V_f to nodes
	host_to_node_real_1(V_f);

#if !RP_NODE
	Message("Updated FSR to %g m/s using P-control with Kp = %g and R = %g K\n", V_f, Kp, R);
#endif
}


DEFINE_INIT(init_RP_vars, d)
{
#if !RP_NODE
	// Check for V_f initializer value
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

	// Check for UPDATE INTERVAL value

	bool UPDATE_INTERVAL_exists = RP_Variable_Exists_P("user/update_interval");

	Message("Checking for user-defined parameter 'user/update_interval': %d\n", UPDATE_INTERVAL_exists);

	if (UPDATE_INTERVAL_exists)
	{
		UPDATE_INTERVAL = RP_Get_Integer("user/update_interval");
		Message("User-defined parameter 'user/update_interval' found with value: %d\n", UPDATE_INTERVAL);
	}
	else
	{
		UPDATE_INTERVAL = 100; // Default initial FSR value if user-defined parameter does not exist
		Message("Warning: User-defined parameter 'user/update_interval' not found. Using default value of 100.\n");
	}

	bool Kp_exists = RP_Variable_Exists_P("user/kp");

	if (Kp_exists)
	{
		Kp = RP_Get_Real("user/kp");
		Message("User-defined parameter 'user/kp' found with value: %g\n", Kp);
	}
	else
	{
		Message("Warning: User-defined parameter 'user/kp' not found. Using current value of %g.\n", Kp);
	}

#endif

	// Pass parameters to nodes
	//node_to_host_real_1(V_f); // update V_f on host process so report is correct.
	host_to_node_real_1(V_f);
	Message0("V_f initialized to %g m/s \n",V_f);

	host_to_node_int_1(UPDATE_INTERVAL);
	Message0("UPDATE_INTERVAL initialized to %d iterations \n",UPDATE_INTERVAL);

	host_to_node_real_1(Kp);
	Message0("Kp initialized to %g \n",Kp);
}


// Update Solid Motion 
DEFINE_ZONE_MOTION(update_solid_motion, omega, axis, origin, velocity, current_time, dtime)
{
	NV_D(velocity, =, V_f, 0.0, 0.0); // Update solid motion velocity with calculated FSR
	//Message("Updated solid motion to %g m/s\n", V_f);
	return;
}

// Calculate new regressed surface coordinates (move to DEFINE_GRID_MOTION when debugged)
DEFINE_ON_DEMAND(calc_regression)
{
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
		// Supposedly this macro is in para.h, but it does not see to be there anymore.
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

		qsort_s(idx_f, size, sizeof(int), compare_xf, x_f);

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

		// allocate y_f_new on node zero 
		y_f_new = malloc(size * sizeof(real));

		for (int k = 0; k < size; k++)
		{
			I += mdot_f[k] * dx_f[k] / sqrt(pow(rho * V_f * A_f[k], 2) - pow(mdot_f[k], 2));
			//I += mdot_f[k]*nhat_x[k] / A_f[k] / sqrt(pow(rho * V_f * A_f[k], 2) + pow(mdot_f[k], 2));

			y_f_new[k] = r_eig[1] + I;

			Message0("y_f_new(x = %g m) = %g m \n", x_f[k], y_f_new[k]);
		}	

	}

	// Free allocated memory on all nodes
	free(x_f);
	free(y_f);
	free(mdot_f);
	free(A_f);
	free(dx_f);
	//free(idx_f);
	
#endif 
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

// Report Definition for Residual of eigenvalue-based FSR calculation
DEFINE_REPORT_DEFINITION_FN(R_eigen_temp_diff)
{
	return R; // Return residual for report definition
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
#endif

	//node_to_host_real_1(V_f); // update V_f on host process so report is correct.
	host_to_node_real_1(V_f);
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
		Kp = 5e-8; // Default initial FSR value if user-defined parameter does not exist
		Message("Warning: User-defined parameter 'user/kp' not found. Using default value of 5e-8.\n");
	}

#endif

	// Pass parameters to nodes
	host_to_node_real_1(Kp);
	Message0("Kp initialized to %g \n",Kp);
}

DEFINE_ON_DEMAND(set_UPDATE_INTERVAL)
{
#if !RP_NODE
	// Check for UPDATE INTERVAL value
	bool UPDATE_INTERVAL_exists = RP_Variable_Exists_P("user/update_interval");

	Message("Checking for user-defined parameter 'user/update_interval': %d\n", UPDATE_INTERVAL_exists);

	if (UPDATE_INTERVAL_exists)
	{
		UPDATE_INTERVAL = RP_Get_Integer("user/update_interval");
		Message("User-defined parameter 'user/update_interval' found with value: %d\n", UPDATE_INTERVAL);
	}
	else
	{
		UPDATE_INTERVAL = 100; // Default initial FSR value if user-defined parameter does not exist
		Message("Warning: User-defined parameter 'user/update_interval' not found. Using default value of 100.\n");
	}
#endif

	// Pass parameters to nodes
	host_to_node_int_1(UPDATE_INTERVAL);
	Message0("UPDATE_INTERVAL initialized to %d iterations \n",UPDATE_INTERVAL);
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
		eigen_face_zoneID = 0; // Default value if user-defined parameter does not exist, update with different default if desired
		Message("Warning: User-defined parameter 'user/eigen_zone_id' not found. Using default value of 0.\n");
	}
#endif
	//node_to_host_int_1(eigen_face_zoneID); // update eigen_face_zoneID on host process so it can be used in calc_FSR_eigen
	
	host_to_node_int_1(eigen_face_zoneID); //broadcast to nodes
	Message("eigen_face_zoneID initialized to %d \n", eigen_face_zoneID);
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

	host_to_node_real_1(V_f);
}

