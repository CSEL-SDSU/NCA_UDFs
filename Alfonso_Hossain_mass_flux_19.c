
 #include "udf.h"



/* CODE SECTION */
/* FD INLET VELOCITY PROFILE */ 

 DEFINE_PROFILE(inlet_x_vel_10cms, thread, position) 
 {
    real x[ND_ND]; /* this will hold the position vector */
    real y, h, U_mean, U_max, m, n;
    face_t f;

	h = 0.00495; /* m; inlet height, do not change */
	m = 27.59596236; /* constant, do not change */ 
	n = 2.0; /* constant, do not change */

	U_mean = 0.10; /* m/sec; inlet mean velocity, update with geom */
	U_max = U_mean*((m+1)/m)*((n+1)/n); /* m/sec; max velocity, at centerline... calc */ 

    begin_f_loop(f,thread)
    {
      F_CENTROID(x, f, thread);
      y = 2.*(x[1]-0.5*h)/h; /* non-dimensional y coordinate, b/c coord sys is at bottom of geom not centerline... calc */

      F_PROFILE(f, thread, position) = U_max*(1.0-y*y); /* m/sec; velocity as f(y) at centerline... calc */
    }
    end_f_loop(f, thread)
 }


int counter = 0;

 DEFINE_PROFILE(fuel_mass_flux, thread, position)
 {

	real K_store;

	face_t f;
	Thread *t0 = NULL; 
	cell_t c0 = -1;

	real T1, T2; 
	real xc[ND_ND];
	real xf[ND_ND];
	real m_dot=0.0;
	real m_dot_mid=0.0; 
	real q_net=0.0;
	real x1,y1,x2,y2,dist;
	real q_manual,grad;
	real par_C;
	real x_flame;

if (counter < 2) //for first 2 cells celt m'' = 0.02
{
	begin_f_loop(f,thread)
    {
	F_PROFILE(f, thread, position) = .02;
	

    }
    end_f_loop(f, thread)

	/* print iteration count to command window */ 
	counter += 1.0; 
	printf("Counter: %d\n", counter); /* \n denotes a new line */
}

else  // for the remaining cells use profile
{
    begin_f_loop(f,thread)
    {

	c0 = F_C0(f,thread); // get cell ID adjacent to face
	t0 = THREAD_T0(thread); // get thread that cell belongs to (cell center may not be in same thread as face)
	K_store = C_K_L(c0,t0); // get thermal conductivity of current cell
	C_UDMI(c0,t0,0) = K_store; // save k as a field in fluent managed memory (make k a field function)
	

	T1 = C_T(c0, t0); // get temperature at the cell centroid
	C_UDMI(c0,t0,1) = T1; //store cell temperature in user defined memory
	T2 = F_T(f,thread); // get the temperature at the cell face
	C_UDMI(c0,t0,2) = T2; //store the temperature at each face in user defined memory

	
		
	
	C_CENTROID(xc, c0, t0); //get the cell centroid coordinates of cell w/ id = c0 on thread/boundary t0
	F_CENTROID(xf, f, thread); //get the face centroid coordinates

	// compute temperature gradient 
	x1=xc[0];
	y1=xc[1];
	x2=xf[0];
	y2=xf[1];
	dist = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
	grad=(T1-T2)/dist;

	// compute heat flux 
	q_manual=K_store*grad;
	
	// not sure what this is (some sort of unit conversion to mm probably) 
	x_flame=1000*x2-19;

	if(x_flame <= 0)
	{
		m_dot=0;
		Message("x menor que 0\n"); //x menor que 0 means x less than zero
	}
	else
	{
		m_dot=0.01998 * pow(x_flame, -0.59) / sqrt(1 + 0.06672 		* pow(x_flame, -1.18));
	}

	C_UDMI(c0,t0,3) = dist;
	C_UDMI(c0,t0,4) = grad;
	C_UDMI(c0,t0,5) = q_manual;
	C_UDMI(c0,t0,6) = m_dot;
	C_UDMI(c0,t0,7) = par_C;
	C_UDMI(c0,t0,8) = counter;
	F_PROFILE(f, thread, position) = m_dot;

	
    }
    end_f_loop(f, thread)

	/* print iteration count to command window */ 
	counter += 1.0; 
	Message("Counter: %d\n", counter); /* \n denotes a new line */

 }
 }







int counter_2 = 0;
 DEFINE_PROFILE(fuel_heat_loose, thread, position)
 {
		
	real K_store;

	face_t f;
	Thread *t0 = NULL; 
	cell_t c0 = -1;

	real T1, T2; 
	real xc[ND_ND];
	real xf[ND_ND]; 
	real q_manual=0.0;
	real x1,y1,x2,y2,dist;
	real grad;
	real par_C;
	real q_real_out;

	real x_flame;
	real esp = 0.00001;
	real q_less=0.0;
	real q_less_mid=0.0;
		
if (counter_2 < 2)
{
	begin_f_loop(f,thread)
	{
	F_PROFILE(f, thread, position) = 0;
	
 	}
	end_f_loop(f, thread)
counter_2 += 1.0; 
printf("Counter_2: %d\n", counter_2); /* \n denotes a new line */
}

else  
{
	begin_f_loop(f,thread)
	{
	c0 = F_C0(f,thread); 
	t0 = THREAD_T0(thread); 
	K_store = C_K_L(c0,t0); 
	T1 = C_T(c0, t0);
	T2 = F_T(f,thread);	
	C_CENTROID(xc, c0, t0);
	F_CENTROID(xf, f, thread);
	x1=xc[0];
	y1=xc[1];
	x2=xf[0];
	y2=xf[1];
	dist = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
	grad=(T1-T2)/dist;
	q_manual=K_store*grad;

	x_flame=1000*x2-19;
	if(x_flame <= 0)
	{
		q_less=0;
		Message("x menor que 0\n");
	}
	else
	{
		q_real_out=16058.4 * pow(x_flame, -0.59) / sqrt(1 + 0.06672 * pow(x_flame, -1.18));
		q_less=-q_real_out/esp;

	}
	
	C_UDMI(c0,t0,9) = q_less;
	C_UDMI(c0,t0,10) = q_real_out;
	F_PROFILE(f, thread, position) = q_less;
	}
end_f_loop(f, thread)


counter_2 += 1.0; 

}
}

