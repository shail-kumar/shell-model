#include "shell.h"

#define MTAG1 1
#define MTAG2 2

using namespace blitz;

// ------------------- Global constructs --------------------------------------
// typedef complex<double> cdouble;
int attempt;					// Run number
string input_path;
string output_path;
string input_field_file;

//++++++++++++++++++++++++ Control Parameters +++++++++++++++++++++++++++++++++
// const long double PI = 3.141592653589793238L;

const double k_0=1.0/16;								
const double q = 0.5*(1+sqrt(5));	// Shell spacing
int N;	// # of Shells
double nu;						// Viscosity		
double dt;						// time step
int T;							// # of steps at an interval of dt
double epsilon_u;				// Energy input in u
double epsilon_v;				// Energy input in v
int n_f;						// Forced shell
int print_steps;				// Write the data in file at an interval of these many steps
double Dt;						// Time step of printing field
int IC_scheme;					// Initial condition scheme: 0 for random; 1 from  a file 	
int nic;						// Total number of initial conditions
double Omega;					// Rotation Strength
int perturbation_nod;			// Generate perturbed field: 0 for NO; 1 for YES
double perturbation_amp;		// Maximum amplitude of perturbation


// =====================I/O Files =============================================
ifstream para_file; ifstream environment_file; fstream initial_field_file; 
ofstream perturbed_field_file; ofstream final_field_file; ofstream out_field; 
ofstream out_modU; ofstream out_modV; ofstream energy_file;

// =====================main program ==========================================
int main()
{
	time_t wtime = time(0); 

	int process_id, num_procs, islave;
	MPI_Status status;
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_id);

	// ----------Setting environment--------------------------------------
	environment_file.open("environment.d");
	set_environment(environment_file);
	environment_file.close();
	//--------Read parameters---------------------------------------------
	stringstream path;
	path<<input_path<<"parameters.d";
	para_file.open(path.str().c_str());
	parameter_reading(para_file);
	para_file.close();

	//-------Construction of shells----------------------------------------
	Array<double,1> *kn;
	kn=new Array<double,1>(N+2);
	shell_wavenumbers(kn);
	// cout<<"shell_wavenumbers"<<*kn<<endl;

	// ---------- Velocity Fields -----------------------------------------
	Array<cdouble,1> *u;
	u= new Array<cdouble,1>(N+4); 
	Array<cdouble,1> *v;
	v= new Array<cdouble,1>(N+4);
	// ---------- Perturbed Velocity Fields -------------------------------
	Array<cdouble,1> *U;
	U= new Array<cdouble,1>(N+4); 
	Array<cdouble,1> *V;
	V= new Array<cdouble,1>(N+4);

	//---------------Force-------------------------------------------------
	Array<cdouble,1> *force;
	force=new Array<cdouble,1>(N+4);

	// -------------For output---------------------------------------------
	Array<double,1> *mod_u;
	mod_u= new Array<double,1>(N+2);
	Array<double,1> *mod_v;
	mod_v= new Array<double,1>(N+2);

	stringstream warning;

	Dt = dt*print_steps; 
	int data_points = 2*(N+4);			// Number of data points to be sent and receive in MPI routine, 1 complex has 2 data points
	// --------------------------------------------------------------------

	if (process_id==0) // ===== COMPUTATION ENVIRONMENT IS SETUP BY THE MASTER PROCESS HERE. =======
	{
		cout<<"# The program begins on "<<ctime(&wtime)<<endl;

		cout<<"INPUT PATH: "<<input_path<<endl;
		cout<<"OUTPUT PATH: "<<output_path<<endl;

			// ------------Create directories and file------------
		int dir;
		stringstream dir_field;
		dir_field<<output_path<<"velocity_field_"<<attempt;
		dir = mkdir(dir_field.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1==dir)
			{
				cout<<dir_field.str().c_str()<<" exists beforehand."<<endl; 
				cout<<"Overwriting avoided. All processes exit..\n"; MPI_Abort(MPI_COMM_WORLD,1);
			}
		else {cout<<dir_field.str().c_str()<<" created."<<endl;}

/*		stringstream dir_modU;
		dir_modU<<output_path<<"velocity_"<<attempt;
		dir = mkdir(dir_modU.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1==dir)
			{
				cout<<dir_modU.str().c_str()<<" exists beforehand."<<endl;
				cout<<"Overwriting avoided. All processes exit..\n"; MPI_Abort(MPI_COMM_WORLD,1);
			}
		else {cout<<dir_modU.str().c_str()<<" created."<<endl;}*/


/*		stringstream dir_Energy;
		dir_Energy<<output_path<<"energy_"<<attempt;
		dir = mkdir(dir_Energy.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1==dir)
			{
				cout<<dir_Energy.str().c_str()<<" exists beforehand."<<endl;
				cout<<"Overwriting avoided. All processes exit..\n"; MPI_Abort(MPI_COMM_WORLD,1);
			}
		else {cout<<dir_Energy.str().c_str()<<" created."<<endl;}*/

		if (perturbation_nod==1)
		{
			stringstream tmp;
			tmp<<input_path<<"p_field_"<<attempt<<".d"; 
			perturbed_field_file.open(tmp.str().c_str());
			cout<<"Perturbed field is generated and is put in file "<<tmp.str().c_str()<<endl;

		}
		else
		{
			cout<<"Perturbed field is not generated."<<endl;
		}
				//-------READ INITIAL CONDITIONS------------
		switch(IC_scheme)
		{
			case 0 :
				for (islave = 1; islave < num_procs; ++islave)
				{
					u_initial(u, kn, islave); //By generating randonm numbers
					u_initial(v, kn, islave + 200); // Assuming that the number of processes will not exceed 200.
					MPI_Send(u->data(), data_points, MPI_DOUBLE, islave, MTAG1, MPI_COMM_WORLD);				
					MPI_Send(v->data(), data_points, MPI_DOUBLE, islave, MTAG2, MPI_COMM_WORLD);				
				}
				break;
			case 1 :
				path.str("");
				path<<input_path<<input_field_file;
				initial_field_file.open(path.str().c_str());
				// cout<<path.str().c_str();
				// initial_field(u, initial_field_file);			// Reading from a file
				for (islave = 1; islave < num_procs; ++islave)
				{
					initial_field_file >> (*u);
					initial_field_file >> (*v);
					MPI_Send(u->data(), data_points, MPI_DOUBLE, islave, MTAG1, MPI_COMM_WORLD);				
					MPI_Send(v->data(), data_points, MPI_DOUBLE, islave, MTAG2, MPI_COMM_WORLD);

									// -------- PERTURBED INITIAL VELOCITY FIELD-------------
					if (perturbation_nod==1)
					{
						u_perturbed(u, U, perturbation_amp, islave); 
						u_perturbed(v, V, perturbation_amp, islave+200); 
				
						perturbed_field_file << *U << *V << endl;
					}
									
				}
				initial_field_file.close();
				if (perturbation_nod==1){perturbed_field_file.close();}
				break;
			default :
				cout<<"Choose proper IC_scheme in para_file. Remove the created directories before next attempt.\n"; MPI_Abort(MPI_COMM_WORLD,1);

		}
		// cout<<"# intial field u \t"<<*u<<endl;
		// cout<<"# intial field v \t"<<*v<<endl;

		// cout<<"Program terminated on purpose..."; MPI_Abort(MPI_COMM_WORLD,1);

			// --------WRITING FINAL OUTPUT STATUS ------------------------
		cout<<"NOTE: The number of initial conditions is controlled by the number of processes."<<endl;
		nic = num_procs -1;
		cout<<"# Shell spacing: "<<q<<", # Number of Shells: "<<N<<", Viscosity: "<<nu<<", time step-size: "<<dt<<\
		", Total # of steps: "<<T<<", energy input in u: "<<epsilon_u<<", energy input in v: "<<epsilon_v<<", forced shellnumber: "<<n_f\
		<<", printing steps: "<<print_steps<<",	IC Scheme: "<<IC_scheme<<", # of ICs: "<<nic<<\
		", Rot. Strength: "<<Omega<<endl;

		stringstream final_u;
		final_u<<output_path<<"final_field_"<<attempt<<".d"; 
		final_field_file.open(final_u.str().c_str());

		final_field_file<<"# The program begins on "<<ctime(&wtime)<<endl;
		final_field_file<<"# Shell spacing: "<<q<<", # Number of Shells: "<<N<<", Viscosity: "<<nu<<", time step-size: "<<dt<<\
		", Total # of steps: "<<T<<", energy input in u: "<<epsilon_u<<", energy input in v: "<<epsilon_v<<", forced shellnumber: "<<n_f\
		<<", printing steps: "<<print_steps<<",	IC Scheme: "<<IC_scheme<<", # of ICs: "<<nic<<\
		", Rot. Strength: "<<Omega<<endl;

		final_field_file<<"# shell wavenumbers: \n"<<(*kn)<<endl;

		for (islave = 1; islave < num_procs; ++islave)
		{
			int l;
			MPI_Probe(islave, 0, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_CHAR, &l);
			char *buf = new char[l];
			MPI_Recv(buf, l, MPI_CHAR, islave, 0, MPI_COMM_WORLD, &status);
			string proc_status(buf, l);
			delete [] buf;
			MPI_Recv(u->data(), data_points, MPI_DOUBLE, islave, MTAG1, MPI_COMM_WORLD, &status);
			MPI_Recv(v->data(), data_points, MPI_DOUBLE, islave, MTAG2, MPI_COMM_WORLD, &status);
	
			final_field_file<<islave<<"\t\t"<<*u<<endl;
			final_field_file<<islave<<"\t\t"<<*v<<endl;
			final_field_file<<proc_status<<endl;

		}
		time_t wtime = time(0);
		final_field_file<<"\n# Program successfully exits on "<<ctime(&wtime)<<endl;
		final_field_file.close(); 
		cout<<"\n# Program gracefully exits on "<<ctime(&wtime)<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	else // ============ COMPUTATION IS DONE BY OTHER PROCESSES HERE. ===================
	{
		MPI_Recv(u->data(), data_points, MPI_DOUBLE, 0, MTAG1, MPI_COMM_WORLD, &status);
		MPI_Recv(v->data(), data_points, MPI_DOUBLE, 0, MTAG2, MPI_COMM_WORLD, &status);
		// cout<<process_id<<"\t"<<*u<<"\t"<<*v<<endl;
		// out_field<<process_id<<endl;
		stringstream field_evolution;
		field_evolution<<output_path<<"velocity_field_"<<attempt<<"/u_vs_t_"<<process_id<<".d"; 
		out_field.open(field_evolution.str().c_str());
		out_field<<0<<"\t"<<*u<<"\n"; // Writing field at t = 0.

		stringstream modU;
		modU<<output_path<<"velocity_field_"<<attempt<<"/real_u_"<<process_id<<".d"; 
		out_modU.open(modU.str().c_str());
		stringstream modV;
		modV<<output_path<<"velocity_field_"<<attempt<<"/imag_u_"<<process_id<<".d"; 
		out_modV.open(modV.str().c_str());
/*		// Write wavenumbers for a given IC
		out_modU<<setw(10)<<0<<"\t\t"; writeRealData(kn,out_modU);
		out_modV<<setw(10)<<0<<"\t\t"; writeRealData(kn,out_modV);*/
		*mod_u = 0.0;
		*mod_v = 0.0;
		for (int j = 2; j < N+2; ++j)
		{
			(*mod_u)(j)=real((*u)(j));
			(*mod_v)(j)=imag((*u)(j));
		}
		out_modU<<setw(10)<<0<<"\t\t"; writeRealData(mod_u,out_modU);
		out_modV<<setw(10)<<0<<"\t\t"; writeRealData(mod_v,out_modV);


			  //================= MAIN COMPUTATION: TIME EVOLUTION OF SHELL MODEL =============

		double begin = MPI_Wtime();
		for(int Dt_index = 1; Dt_index <= T; Dt_index++)
		{
			for (int dt_index = 1; dt_index <= print_steps; ++dt_index)
			{
				Time_evolution(u,v,force, kn,epsilon_u); // Evolution of u
			}
			out_field<<Dt_index*Dt<<"\t"<<*u<<"\n";		// Writing blitz array to file
			if (process_id==1)
			{
				if(Dt_index%10==0){cout<<Dt_index<<"\n";}
			}
			for (int j = 2; j < N+2; ++j)
			{
				(*mod_u)(j)=real((*u)(j));
				(*mod_v)(j)=imag((*u)(j));
			}
			out_modU<<setw(10)<<Dt_index*Dt<<"\t\t"; writeRealData(mod_u,out_modU);
			out_modV<<setw(10)<<Dt_index*Dt<<"\t\t"; writeRealData(mod_v,out_modV);

			if (isnan(abs((*u)(2))) || isnan(abs((*v)(2))))
			{
				cout<<"# Termination step: "<<Dt_index<<"\t\t in process "<<process_id<<"\n # NUMERICAL OVERFLOW! Process exits with attempt: "<<attempt<<endl;
				double end = MPI_Wtime();
				warning<<"# Termination step: "<<Dt_index<<"\t\t in process "<<process_id<<"\n # NUMERICAL OVERFLOW! Process exits with attempt: "<<attempt<<" . \n # Program terminated in time = "<<(end - begin)/60<<" minutes \n ";
				MPI_Send(warning.str().c_str(), warning.str().size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
				MPI_Send(u->data(), data_points, MPI_DOUBLE, 0, MTAG1, MPI_COMM_WORLD);				
				MPI_Send(v->data(), data_points, MPI_DOUBLE, 0, MTAG2, MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				exit(0);
			}
		}
		double end = MPI_Wtime();
		warning<<"# Process "<<process_id<<" exited successfully. Computation time: "<<(end - begin)/60<<" minutes\n";
		MPI_Send(warning.str().c_str(), warning.str().size(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);

		out_modU.close();
		out_modV.close();
		out_field.close();
		energy_file.close();
		MPI_Send(u->data(), data_points, MPI_DOUBLE, 0, MTAG1, MPI_COMM_WORLD);				
		MPI_Send(v->data(), data_points, MPI_DOUBLE, 0, MTAG2, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}
/*===============================Main Program ends here====================*/
