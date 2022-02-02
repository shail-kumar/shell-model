/*====================================================================================
Compute V(t+dt) using RK4 scheme

====================================================================================*/ 
#include "shell.h"
#include "evolution.h"
using namespace blitz;
typedef complex<double> cdouble;

Array<cdouble,1> *nlin;
void Time_evolution(Array<cdouble,1>* V,Array<cdouble,1>* W, Array<cdouble,1>* force, Array<double,1>* k, double E_input_rate) {

	static Array<cdouble,1> Vcopy(N+4);
	static Array<cdouble,1> tot_Vrhs(N+4);
	nlin=new Array<cdouble,1>(N+4);

	Vcopy = *V;
	tot_Vrhs = 0.0;		
	// int flag;	
	// cout<<"Initial\n"<<*V<<endl;		

	// Step 1
	Compute_nlin_helical(V, W, nlin, k); 
	Compute_force(V,force, E_input_rate);
	Compute_rhs(nlin,force);
	// cout<<"1\n"<<*nlin<<endl;
	Single_time_step_EULER(V, nlin, k, dt/2);						// u1: Go to the mid point
	// cout<<"euler\n"<<*V<<endl;
	Mult_nlin_exp_ksqr_dt(nlin, k, dt);								// *nlin = *nlin x exp(-nu k^2 dt)		
	tot_Vrhs = tot_Vrhs + *nlin;
	// cout<<"2\n"<<*nlin<<endl;
	
	// Step 2
	// Compute_force(V,force);

	Compute_force(V,force, E_input_rate);

	// else {Compute_force(V,force);}

	// Compute_nlin(V,nlin,k);											// Compute nlin using V(t+dt/2)
	Compute_nlin_helical(V, W, nlin, k); 
	// cout<<"3\n"<<*nlin<<endl;
	Compute_rhs(nlin,force);
	*V = Vcopy;
	Single_time_step_Semi_implicit(V, nlin, k, dt/2);							// u2: Goto the mid point
	// cout<<"semi\n"<<*V<<endl;
	Mult_nlin_exp_ksqr_dt(nlin, k, dt/2);										// rhs(t_mid, u1) in nlin	
	tot_Vrhs = tot_Vrhs + 2.0* (*nlin);
		
	// Step 3
		
	// Compute_force(V,force);
	Compute_force(V,force, E_input_rate);

	// else {Compute_force(V,force);}

	Compute_nlin_helical(V, W, nlin, k); 
	// cout<<"4\n"<<*nlin<<endl;
	Compute_rhs(nlin,force);
	*V = Vcopy;
	Single_time_step_RK2(V, nlin, k, dt);								// u3: Go to the end point
	// cout<<"RK2\n"<<*V<<endl;		
	Mult_nlin_exp_ksqr_dt(nlin, k, dt/2);									
	tot_Vrhs = tot_Vrhs + 2.0* (*nlin);						// rhs(t_mid, u2) in nlin
		
	// Step 4
		
	// Compute_force(V,force);
	Compute_force(V,force, E_input_rate);

	// else {Compute_force(V,force);}

	Compute_nlin_helical(V, W, nlin, k); 
	// cout<<"5\n"<<*nlin<<endl;
	Compute_rhs(nlin,force);
										
	tot_Vrhs = tot_Vrhs + (*nlin);							// rhs(t_end, u3) in nlin
		

	// Final result
				
	*V = Vcopy;
	Mult_field_exp_ksqr_dt(V, k, dt);
				
	 *V = *V + (dt/6) * (tot_Vrhs);
	// cout<<"Final\n"<<*V<<endl;
	delete nlin;		

}

//====================== fn time_advance ends ====================
   

