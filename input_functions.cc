#include "shell.h"
using namespace blitz;

typedef complex<double> cdouble;
const long double PI = 3.141592653589793238L;


/*void parameter_reading(ifstream& para_file, int &no_of_shells, double &nu, double &dt, int &T,\
 cdouble &f_amp, int &n_f, int &print_steps, int &IC_scheme)*/
void parameter_reading(ifstream& para_file)
{
	string useless_line;
	if (para_file.is_open())
	{
//		getline(para_file, useless_line);
//		para_file>>k_0; getline(para_file, useless_line);
//		getline(para_file, useless_line);
//		para_file>>q; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>N; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>nu; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>dt; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>T; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>epsilon_u; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>epsilon_v; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>n_f; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>print_steps; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>IC_scheme; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>nic; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>Omega; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>perturbation_nod; getline(para_file, useless_line);
		getline(para_file, useless_line);
		para_file>>perturbation_amp;
		// para_file>>no_of_shells>>nu>>dt>>T>>f_amp>>n_f>>print_steps>>IC_scheme;
		// cout<<no_of_shells<<"		"<<nu<<"		"<<dt<<"		"<<T<<"		"<<endl;
	}
	else
	{
		cout<<"Parameter file is not open."<<endl;
	}

}

void set_environment(ifstream &env_file)
{
	string useless_line;
	if (env_file.is_open())
	{
		getline(env_file, useless_line);
		env_file>>attempt; getline(env_file, useless_line);
		getline(env_file, useless_line);
		env_file>>input_field_file;  getline(env_file, useless_line);
		getline(env_file, useless_line);
		env_file>>input_path; getline(env_file, useless_line);
		getline(env_file, useless_line);
		env_file>>output_path; getline(env_file, useless_line);
	}
	else
	{
		cout<<"Environment file is not open."<<endl;
	}

}

void shell_wavenumbers(Array<double,1>* k_n)
{
	for (int i = 0; i < N+2; ++i)
	{
		(*k_n)(i)=k_0*pow(q,i-1);			//i=2 corressponds to 1st shell
	}
	// cout<<1<<(*k_n)<<endl;
	// return (*k_n);
}

void u_initial(Array<cdouble,1>* u0, Array<double,1>* k_n, int seed)
{
	// int seed;
	srand(time(0) +seed);
	for(int i=2;i<N+2;i++){
		// u(i)=1.0*rand()/(RAND_MAX +1.0);
		// (*u0)(i) = pow((*k_n)(i),2);
		double theta = 2*PI*rand()/(RAND_MAX + 1.0);
		if (i<5)
		{
			// (*u0)(i) = cdouble (cos (theta), sin (theta))*sqrt((*k_n)(i));
			(*u0)(i) = cdouble (cos (theta), sin (theta))*pow((*k_n)(i),-1.0/3);
		}
		else
		{
			(*u0)(i) = cdouble (cos (theta), sin (theta))*pow((*k_n)(i),-1.0/3);
			// (*u0)(i) = cdouble (cos (theta), sin (theta))*sqrt((*k_n)(i))*exp(-pow((*k_n)(i),2));
			// (*u0)(i) = cdouble (cos (theta), sin (theta))*pow((*k_n)(i),0.5)*exp(-pow((*k_n)(i),2));

		}

	}

	// return (*u0);	
}

void u_perturbed(Array<cdouble,1>* u0, Array<cdouble,1>* U0, double AMPLITUDE, int seed)
{
	// int seed;
	srand(time(0) +seed);
	for(int i=2;i<N+2;i++){
		// double theta = 2*PI*rand()/(RAND_MAX + 1.0);
		double theta = AMPLITUDE*rand()/(RAND_MAX + 1.0);
		if (i<26)
		{
			(*U0)(i) = (*u0)(i);
		}
		else
		{
			(*U0)(i) = cdouble (cos (theta), sin (theta))*(*u0)(i);
		}

	}

	// return (*u0);	
}

void initial_field(Array<cdouble,1>* u0, fstream& initial_field_file)
{
	initial_field_file >> (*u0);
	cout  << "Reading of field configurations ended successfully" << endl; 
}
