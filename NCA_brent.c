
#include "udf.h"
#include <stdbool.h>
#include "hdfio.h"
#include <math.h>

#define EPS 2.2204460492503131e-16
#define sgn(x)  ((x>0) - (x<0))

// State flags
typedef enum {
	EVAL_FA = 0,
	EVAL_FB = 1,
	EVAL_ITER = 2,
	CONVERGED = 3,
	FAILED = 4
} state_t;

static state_t STATE = EVAL_FA;
static real V_f; //Flame Spread Rate 

//Variables required for evaluating F(X) - The residual temperature
// And determining when an evaluation is accepted
static real R;
static int eigen_face_zoneID = 0; // Zone ID of separated surface where temp is monitored
static const real T_infty = 300; 

#define R_NP 1000
#define R_HISTORY_SIZE (R_NP + 1)
static int N_NEWTON_ITER = 0;
static int N_LAST_UPDATE = -1;
static const int MIN_ITER_BEFORE_UPDATE = 1000;

static real R_relative_tolerance = 1e-5; // Residual tolerance for accepting a given evaluation of the function R(V_f), adjust as needed
static real R_denominator_floor = 1e-12;

static real R_history[R_HISTORY_SIZE];
static int R_history_count = 0;
static int R_history_head = 0;

static real R_relative_residual_max = 1;

// Temp based flame spread rate calculation variables
// Starting bracket [AX,BX], Assumed that F(AX) and F(BX) have opposite signs 
// B is latest iterate and closest approximation to zero
// A is previus iterate
// C is the previous iterate or an older iterate so that F(B) and F(C) have opposite signs
static real A = 0.0,B = 0.000500, C, D, E, FA, FB, FC;
static real TOL = 1e-6; // Tolerance for accepting the final V_f (1 um/s)
static int N_BRENT = 0;

// Function to clear the storage of R
static void reset_R_history(void)
{
#if !RP_HOST // run on compute nodes in parallel, and single process in serial
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
#if !RP_HOST //Store R in history. Stored on compute nodes in parallel
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

	//Update node values (should be same on all nodes to host)
	node_to_host_int_2(R_history_head,R_history_count);
	node_to_host_real(R_history, R_HISTORY_SIZE);

}

// function to check R relative change. Determines whether the largest R
// residual over the last 1000 iterations is less than the tolerance 
static bool R_relative_change_below(real threshold)
{
	int R_converged = 0;

#if !RP_HOST // Run this on compute nodes in parallel or single process in serial
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

	// Pass node values to host. Should be the same on all nodes if in parallel
	node_to_host_int_1(R_converged);
	node_to_host_real_1(R_relative_residual_max);

	return R_converged ? true : false;
}


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

//Make sure this function is above 'update_FSR_newton' in define at end list
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
	store_R_in_history(R);

	// Report using host process or single serial process 
#if !RP_NODE
	if (n_eig_faces != 1 && N_ITER % 25 == 0)
	{
		Message("Warning: eigen face zone %d has %d principal faces. T_eig = %f K\n",
			eigen_face_zoneID, n_eig_faces, T_eig);
	}

	if (N_ITER % 25 == 0)
	{
		Message("NEWTON_ITER=%d, V_f = %g m/s, T_eig = %g K, R = %f K, "
			"R_rel_res_max = %g, R_hist_count = %d\n",
			N_NEWTON_ITER, V_f, T_eig, R,
			R_relative_residual_max, R_history_count);
	}
#endif
}


DEFINE_EXECUTE_AT_END(update_FSR_brent)
{
	int R_is_stationary = 0;
	int return_flag = 0;
	int did_update = 0;

	// Check if the evaluation of the function is changing relative to the last 1000 iterations
	R_is_stationary = R_relative_change_below(R_relative_tolerance);

	// Check if minimum number of iterations has passed. 
#if !RP_NODE //run on host in parallel or on single process in serial

	if (STATE == CONVERGED || STATE == FAILED)
	{
		return_flag = 1;
		//return;
	}
#endif

	host_to_node_int_1(return_flag);	
	if(return_flag) {return;}

#if !RP_NODE
	// Flags, only needed on host
	int enough_iters_since_update = 0;	
	int interpolation_flag = 0;
	

	enough_iters_since_update = ((N_ITER - N_LAST_UPDATE) > MIN_ITER_BEFORE_UPDATE);

	if (!R_is_stationary ||	!enough_iters_since_update)
	{
		return_flag = 1;
		//return;
	}
#endif

	host_to_node_int_1(return_flag);	
	if(return_flag) {return;}

#if !RP_NODE
	// Initialization of bracket (getting the first two evaluations)
	if (STATE == EVAL_FA)
	{
		FA = R;
		V_f = B; // Change V_f from V_f = A to V_f = B
		STATE = EVAL_FB;

		N_LAST_UPDATE = N_ITER;
		did_update = 1;
	}
	else if(STATE == EVAL_FB)
	{
		FB = R;
		STATE = EVAL_ITER;
	}
	else if(STATE == EVAL_ITER)
	{
		FB = R;
	}
#endif

	// Update State on nodes
	host_to_node_int_3(STATE,did_update, N_LAST_UPDATE);
	host_to_node_real_3(V_f, FA, FB);

	if (did_update)
	{
		reset_R_history();
	}

	if (STATE != EVAL_ITER)
	{
		return;
	}

	// If we are here, we are in STATE=EVAL_ITER, we will now compute new values

	// If FA and FB do not have differeing signs, a zero is not 
	// guaranteed by the IVT. Stop the run
#if !RP_NODE //run on host in parallel or on single process in serial

	real P, Q, Rb, S, XM;

	if (FA*FB >= 0 && N_BRENT == 0)
	{
		STATE = FAILED;
		RP_Set_Integer("sol/iterations",0);
		Message("Root is not bracketed between F(A=%g) = %g and F(B=%g) = %g\n", A, FA, B, FB);
		return_flag = 1;
	}

#endif

	host_to_node_int_2(return_flag,STATE);	
	if(return_flag) {return;}

#if !RP_NODE

	if ((FB > 0.0 && FC > 0.0) || (FB < 0.0 && FC < 0.0) || N_BRENT == 0)
	{
		C = A;
		FC = FA;
		D = B - A;
		E = D; //Size of bracket?
	}	

	// Check if sides need to be swapped. F(C) should be larger value and FB should be smallervalue
	if (fabs(FC) < fabs(FB))
	{
		A = B;
		B = C;
		C = A;
		FA = FB;
		FB = FC;
		FC = FA;
	}

	//--- Check for convergence ---
	real TOL1 = TOL + 4*EPS*fabs(B);

	// Calculate bisection point
	XM = 0.5*(C - B);

	// XM = B + 0.5C - 0.5B = )
	if (fabs(XM) <= TOL1 || FB == 0.0)
	{
		// Convergence reached
		STATE = CONVERGED;
		RP_Set_Integer("sol/iterations",0);
		V_f = B;
		Message("FSR Converged! V_f = %g m/s \n", V_f);
		return_flag = 1;
	}
#endif

	host_to_node_real_1(V_f);
	host_to_node_int_2(return_flag, STATE);	
	if(return_flag) {return;}

#if !RP_NODE
	// ---Is bisection neccesary---
	// if current point is is inside interval and not too close to endpoints
	if (fabs(E) < TOL1 || fabs(FA) <= fabs(FB))
	{
		D = XM;
		E = D;
	}
	else if( A != C) // is quadratic interpolation possible?
	{
		Q = FA/FC;
		Rb = FB/FC;
		S = FB/FA;
		P = S * (2.0 * XM * Q * (Q - Rb) - (B - A)*(Rb - 1.0));
		Q = (Q - 1.0) * (Rb - 1.0) * (S - 1.0);
		interpolation_flag = 1;
	}
	else // Linear interpolation (Secant method)
	{
		S = FB/FA;
		P = 2.0 * XM * S;
		Q = 1.0 - S;		
		interpolation_flag = 1;
	}

	// Is interpolation acceptable
	if (interpolation_flag)
	{		

		// Adjust signs 
		if (P > 0.0)
		{
			Q = -Q;
		}
		P = ABS(P);

		if ((2.0 * P) >= (3.0 * XM * Q - fabs(TOL1 * Q)) || P >= fabs(0.5 * E * Q))
		{
			// Interpolation is not acceptable, do bisection
			D = XM;
			E = D;
		}
		else
		{
			E = D;
			D = P/Q;
		}
	}

	// Complete Step 
	A = B;
	FA = FB;
	if (fabs(D) > TOL1)
	{
		B = B + D;
	}
	else if (fabs(D) <= TOL1)
	{
		B = B + fabs(TOL1)*sgn(XM);
	}

	// Report FSR
	N_BRENT++;
	V_f = B;

	Message("N_BRENT=%d, V_f = %g m/s \n",N_BRENT, V_f);

	// Dont think this section is needed since it is already taken care of at the entrance of this UDF
	// checks if FB and FC have same sign
	// // Function evaluation at new B needed
	// //FB = F(B)

	// // if the new evaluation at B 
	// if ((FB * (FC/fabs(FC)) > 0.0))
	// {
	// 	//start over
	// }

#endif 		

	N_LAST_UPDATE = N_ITER;
	reset_R_history();

	//Update A and B on the nodes
	//Update FSR and Brent iterations
	host_to_node_real_7(A, B, C, FA, FB, FC, V_f);
	host_to_node_int_1(N_BRENT);

}

DEFINE_INIT(set_Vf_at_start, d)
{
	// Start with V_f = A
	V_f = A;
	Message0("Set V_f to lower bracket start value. V_f = %g m/s \n",V_f);
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

DEFINE_ON_DEMAND(set_brent_bracket)
{
#if !RP_NODE
	bool AV_f_exists = RP_Variable_Exists_P("user/av_f"); // Check if user-defined parameter exists
	bool BV_f_exists = RP_Variable_Exists_P("user/bv_f");

	bool both_exist = AV_f_exists && BV_f_exists;

	Message0("Checking for user-defined parameters 'user/av_f' and 'user/bv_f'. Both exist: %d\n", both_exist);

	if (both_exist)
	{
		A = RP_Get_Real("user/av_f");
		B = RP_Get_Real("user/bv_f");

		Message("User-defined parameter 'user/avf' found with value: %g\n", A);
		Message("User-defined parameter 'user/bv_f' found with value: %g\n", B);
	}
	else
	{
		//eigen_face_zoneID = 0; // Default value if user-defined parameter does not exist, update with different default if desired
		Message("Warning: User-defined parameter 'user/av_f' or 'user/bv_f' not found. Using default value of A = 0.0, B = 0.000500 \n");
		//Message("Warning: User-defined parameter 'user/eigen_zone_id' not found. Using default value of 0.\n");
	}
#endif

	//node_to_host_int_1(eigen_face_zoneID); // update eigen_face_zoneID on host process so it can be used in calc_FSR_eigen
	host_to_node_real_2(A,B);
	Message("Bracket Initialized. A = %g, B = %g \n", A,B);

	STATE = EVAL_FA;
	N_BRENT = 0;
	V_f = A;
	N_LAST_UPDATE = N_ITER;

	FA = 0.0;
	FB = 0.0;
	FC = 0.0;
	C = A;
	D = B - A;
	E = D;

	host_to_node_real_9(A, B, C, D, E, FA, FB, FC, V_f);
	host_to_node_int_3(STATE, N_BRENT, N_LAST_UPDATE);

	reset_R_history();
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