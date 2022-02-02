#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <blitz/array.h>
#include <sys/stat.h>
#include <omp.h>
#include "mpi.h"
using namespace blitz;

typedef complex<double> cdouble;
extern int attempt;
extern string input_path;
extern string output_path;
extern string input_field_file;

extern const double k_0;								
extern const double q;
extern int N;
extern double nu;						// Viscosity		
extern double dt;				// time step
extern int T;				// # of steps at an interval of dt
extern cdouble f_amp;		// Forcing amplitude
extern double epsilon_u;				// Energy input in u
extern double epsilon_v;				// Energy input in v

extern int n_f;				// Forced shell number
extern int print_steps;
extern int IC_scheme;
extern int nic;
extern double Omega;
extern int perturbation_nod;	// Generate perturbed field: 0 for NO; 1 for YES
extern double perturbation_amp;	// Maximum amplitude of perturbation

// =================Input functions declaratation===============================
// void parameter_reading(ifstream&, int&, double&, double&, int&, cdouble&, int&, int&, int&);
void parameter_reading(ifstream&);
void set_environment(ifstream&);
void shell_wavenumbers(Array<double,1>*);
void u_initial(Array<cdouble,1>*, Array<double,1>*, int);
void initial_field(Array<cdouble,1>*, fstream&);
void Time_evolution(Array<cdouble,1>*, Array<cdouble,1>*, Array<cdouble,1>*, Array<double,1>*, double);
void u_perturbed(Array<cdouble,1>* , Array<cdouble,1>* , double , int);

// =================Output functions declaratation===============================
 void writeRealData (Array<double,1>*, ofstream&);



