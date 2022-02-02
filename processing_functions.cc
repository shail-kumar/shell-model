#include "shell.h"
#include "evolution.h"
using namespace blitz;

typedef complex<double> cdouble;
const cdouble  minusI = cdouble(0,-1.0);			// minus -I
const cdouble  plusI = cdouble(0,1.0);			// plus +I
double a1 = 1.0; double a2=q-2.0; double a3=1.0-q;			
double alpha = 1.0;

/*double a31 = q + pow(q,2); double a32 =-(pow(q,-alpha) + pow(q,alpha));
double a33 = (pow(q,-alpha)-1)/q;
double normal_factor = 1;
*/
double a31 = 1.0; double a32 =-(pow(q,-alpha) + pow(q,alpha))/(pow(q,alpha) + 1);  // From Benzi's Helical shell models 
double a33 = (pow(q,-alpha)-1)/(pow(q,alpha) + 1);									//for three-dimensional turbulence
double normal_factor = (q+q*q);													// multiply for Lessiness

const cdouble f_0 = cdouble(1.0,1.0);

inline int Get_i(int n) {return n+2; }				// returns array index for shell n
inline int Get_n(int i) {return i-2; }				// returns shell index for array index i
	
void Compute_nlin(Array<cdouble,1> *V, Array<cdouble,1> *nlin, Array<double,1> *k) 
{

	(*V)(0) = (*V)(1) = (*V)(N+2) = (*V)(N+3) = 0.0;	
	*nlin = 0.0;
	for (int n=0; n<N; n++) 
		(*nlin)(Get_i(n)) = minusI *(a1* (*k)(Get_i(n))   * conj( (*V)(Get_i(n+2)) * (*V)(Get_i(n+1)) )
						 +  a2* (*k)(Get_i(n-1)) * conj( (*V)(Get_i(n+1)) * (*V)(Get_i(n-1)) )
						 +  a3* (*k)(Get_i(n-2)) * conj( (*V)(Get_i(n-1)) * (*V)(Get_i(n-2)) ));
}

void Compute_nlin_helical(Array<cdouble,1> *V, Array<cdouble,1> *W, Array<cdouble,1> *nlin, Array<double,1> *k) 
{

	(*V)(0) = (*V)(1) = (*V)(N+2) = (*V)(N+3) = 0.0;	// u^+ component
	(*W)(0) = (*W)(1) = (*W)(N+2) = (*W)(N+3) = 0.0;	// u^- component
	*nlin = 0.0;
	for (int n=0; n<N; n++)
	{
		/*(*nlin)(Get_i(n)) = plusI*(*k)(Get_i(n))*conj(a31*(*V)(Get_i(n+1))*(*W)(Get_i(n+2)) \
							+ a32*(*V)(Get_i(n-1))*(*W)(Get_i(n+1)) + a33*(*W)(Get_i(n-1)) * (*W)(Get_i(n-2)));
*/		/*(*nlin)(Get_i(n)) = normal_factor * plusI *(a31* (*k)(Get_i(n)) * conj( (*W)(Get_i(n+2)) * (*V)(Get_i(n+1)) )
						 +  a32* (*k)(Get_i(n-1)) * conj( (*W)(Get_i(n+1)) * (*V)(Get_i(n-1)) )
						 +  a33* (*k)(Get_i(n-2)) * conj( (*W)(Get_i(n-1)) * (*W)(Get_i(n-2)) )); 		// Helical shell model*/
		(*nlin)(Get_i(n)) = plusI * ( a1* (*k)(Get_i(n))   * conj( (*V)(Get_i(n+2)) * (*V)(Get_i(n+1)) )
								 	+  a2* (*k)(Get_i(n-1)) * conj( (*V)(Get_i(n+1)) * (*V)(Get_i(n-1)) )
								 	+  a3* (*k)(Get_i(n-2)) * conj( (*V)(Get_i(n-1)) * (*V)(Get_i(n-2)) ));		// GOY shell model

	} 
}

void Compute_force(Array<cdouble,1> *V, Array<cdouble,1> *force, double EIR)
{
	// EIR is Energy Input Rate
	// (*force)(n_f+1) = f_0*epsilon_u*(*V)(n_f+1)/pow(abs((*V)(n_f+1)),2);
	// (*force)(n_f+1) = epsilon_u*f_0/(real((*V)(n_f+1)) + imag((*V)(n_f+1)));
	*force=0.0;
	for (int n=0; n<N; n++){
		if (n==n_f-1)
		{
			(*force)(Get_i(n)) = f_0*EIR*(*V)(Get_i(n))/pow(abs((*V)(Get_i(n))),2) + minusI*Omega*(*V)(Get_i(n)); // Forcing with Coriolis force
		}
		else 
		{
			(*force)(Get_i(n)) = minusI*Omega*(*V)(Get_i(n));	// Coriolis force term
		}
	}
/*	Array<cdouble,1> *T;
	T = new Array<cdouble,1> (N+4);
	cout<<"Force_u: "<<*force;
	cout<<"V: "<<*V;
	*T = *V + *force;
	cout<<"total V: "<<*T;
	cout<<"V: "<<(*V)(n_f+1);*/

}

void Compute_force_u(Array<cdouble,1> *V, Array<cdouble,1> *force, double EIR)
{
	// EIR is Energy Input Rate
	// (*force)(n_f+1) = f_0*epsilon_u*(*V)(n_f+1)/pow(abs((*V)(n_f+1)),2);
	// (*force)(n_f+1) = epsilon_u*f_0/(real((*V)(n_f+1)) + imag((*V)(n_f+1)));
	*force=0.0;
	for (int n=0; n<N; n++){
		if (n==n_f-1)
		{
			(*force)(Get_i(n)) = f_0*EIR*(*V)(Get_i(n))/pow(abs((*V)(Get_i(n))),2) + minusI*Omega*(*V)(Get_i(n)); // Forcing with Coriolis force
		}
		else 
		{
			(*force)(Get_i(n)) = minusI*Omega*(*V)(Get_i(n));	// Coriolis force term
		}
	}
/*	Array<cdouble,1> *T;
	T = new Array<cdouble,1> (N+4);
	cout<<"Force_u: "<<*force;
	cout<<"V: "<<*V;
	*T = *V + *force;
	cout<<"total V: "<<*T;
	cout<<"V: "<<(*V)(n_f+1);*/

}

void Compute_force_v(Array<cdouble,1> *V, Array<cdouble,1> *force, double EIR)
{
	// (*force)(n_f+1) = f_0*epsilon_v*(*V)(n_f+1)/pow(abs((*V)(n_f+1)),2);
	// (*force)(n_f+1) = epsilon_v*f_0/(real((*V)(n_f+1)) + imag((*V)(n_f+1)));
	*force=0.0;
	for (int n=0; n<N; n++){
		if (n==n_f-1)
		{
			(*force)(Get_i(n)) = f_0*EIR*(*V)(Get_i(n))/pow(abs((*V)(Get_i(n))),2) + minusI*Omega*(*V)(Get_i(n)); // Forcing with Coriolis force
		}
		else 
		{
			(*force)(Get_i(n)) = minusI*Omega*(*V)(Get_i(n));	// Coriolis force term
		}
	}
	// cout<<"Force_v: "<<*force;
/*	Array<cdouble,1> *T;
	T = new Array<cdouble,1> (N+4);
	cout<<"Force_u: "<<*force;
	cout<<"V: "<<*V;
	*T = *V + *force;
	cout<<"total V: "<<*T;
	cout<<"V: "<<(*V)(n_f+1);*/
}

void Compute_rhs(Array<cdouble,1> *nlin, Array<cdouble,1> *force)
{
	*nlin = *nlin + *force;	
}


void Single_time_step_EULER(Array<cdouble,1> *V, Array<cdouble,1> *nlin, Array<double,1> *k, double dt)
{	
	// cout<<10<<*nlin<<endl;										
	*V = *V + dt* (*nlin);
	// cout<<11<<*V<<endl;										
	Mult_field_exp_ksqr_dt(V, k, dt);
	// cout<<13<<*V<<endl;										

}

void Single_time_step_Semi_implicit(Array<cdouble,1> *V, Array<cdouble,1> *nlin, Array<double,1> *k, double dt)
{	
	Mult_field_exp_ksqr_dt(V, k, dt);
	*V = *V + dt* (*nlin);	
	// cout<<"semi impljgd \t\n"<<*V<<endl;	

}
void Single_time_step_RK2(Array<cdouble,1> *V, Array<cdouble,1> *nlin, Array<double,1> *k, double dt)
{	
	Mult_field_exp_ksqr_dt(V, k, dt/2);
	*V = *V + dt* (*nlin);		
	Mult_field_exp_ksqr_dt(V, k, dt/2);
}



void Mult_field_exp_ksqr_dt(Array<cdouble,1> *V, Array<double,1> *k, double dt)
{
	// cout<<12<<*V<<endl;										
	for (int n=0; n<N; n++) 
	{
		int i = Get_i(n);
		(*V)(i) = (*V)(i) * exp(-nu*pow2((*k)(i))* dt);
        // (*V)(i) = (*V)(i) * exp( ( (-nu*pow2((*k)(i))) + (-drag_coeff*  (1/(pow2((*k)(i)))))  ) * dt); //With Drag coefficient //q = 1
        //(*V)(i) = (*V)(i) * exp( ( (-diss_coeff*pow2((*k)(i))) + (-drag_coeff*  (1/(pow4((*k)(i)))))  ) * dt); //With Drag coefficient //q = 1
	}
	// cout<<14<<*V<<endl;										
	
}

//
//

void Mult_nlin_exp_ksqr_dt(Array<cdouble,1> *nlin, Array<double,1> *k, double dt)
{
	for (int n=0; n<N; n++) 
	{
		int i = Get_i(n); 
		//or int i = n+2;
		(*nlin)(i) = (*nlin)(i) * exp(-nu*pow2((*k)(i))* dt);
        // (*nlin)(i) = (*nlin)(i) * exp( ( (-diss_coeff*pow2((*k)(i))) + (-drag_coeff*  (1/(pow2((*k)(i)))))  ) * dt); //With Drag coefficient //q = 1
       //(*nlin)(i) = (*nlin)(i) * exp( ( (-diss_coeff*pow2((*k)(i))) + (-drag_coeff*  (1/(pow4((*k)(i)))))  ) * dt); //With Drag coefficient //q = 1

	}
}


// void compute_E_flux(Array<cdouble,1> *u)
