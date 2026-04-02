
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

	U_mean = 0.15; /* m/sec; inlet mean velocity, update with geom */
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
	real heat_vap=803000;

if (counter < 2)
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

else  
{
    begin_f_loop(f,thread)
    {

	c0 = F_C0(f,thread); 
	t0 = THREAD_T0(thread); 
	K_store = C_K_L(c0,t0); 
	C_UDMI(c0,t0,0) = K_store;
	

	T1 = C_T(c0, t0);
	C_UDMI(c0,t0,1) = T1;
	T2 = F_T(f,thread);	
	C_UDMI(c0,t0,2) = T2;

	
		
	
	C_CENTROID(xc, c0, t0);
	F_CENTROID(xf, f, thread);
	x1=xc[0];
	y1=xc[1];
	x2=xf[0];
	y2=xf[1];
	dist = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
	grad=(T1-T2)/dist;
	q_manual=K_store*grad;
	

	if (q_manual <= 20)
	{
		m_dot=0.02;
	}

	else
	{
			if (counter < 8)
		{
			par_C = 1.7 - 35 * x2;
			if (par_C < 0)
			{
				par_C=0;
 			}
 		}
		else
		{
			par_C = C_UDMI(c0, t0, 7);
			if (T2 > 669.0)
			{
    				par_C += 0.001;
    				if (par_C > 1.0)
    				{
    					par_C = 1.0;
				}
			}
			else if (T2 < 667.0)
			{
   				par_C -= 0.001;  
				if (par_C < 0.0)
   				{
					par_C = 0.0;
				}
			}

		}
		
		
		m_dot_mid= par_C*q_manual/heat_vap;

		if (m_dot_mid > 0.2)
		{
			m_dot = 0.2;
		}
		else
		{
			m_dot = m_dot_mid;
		}
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

	if (q_manual <= 20)
	{
		q_less=0;
	}

	else
	{
		par_C = C_UDMI(c0, t0, 7);
		q_less_mid= -par_C*q_manual/esp;
		q_real_out=-par_C*q_manual;
		if (q_less_mid < -15000000000)
		{
			q_less= -15000000000;
		}
		else
		{
			q_less= q_less_mid ;
		}
	}

	
	C_UDMI(c0,t0,9) = q_less;
	C_UDMI(c0,t0,10) = q_real_out;
	F_PROFILE(f, thread, position) = q_less;
	}
end_f_loop(f, thread)


counter_2 += 1.0; 

}
}

