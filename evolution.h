	// Computaion of RHS terms 
	void Compute_nlin(Array<cdouble,1>*, Array<cdouble,1>*, Array<double,1>*);
	void Compute_nlin_helical(Array<cdouble,1>*, Array<cdouble,1>*, Array<cdouble,1>*, Array<double,1>*); 
	void Compute_force(Array<cdouble,1>*, Array<cdouble,1>*, double);	    
	void Compute_force_u(Array<cdouble,1>*, Array<cdouble,1>*);	    
	void Compute_force_v(Array<cdouble,1>*, Array<cdouble,1>*);	    
	void Compute_rhs(Array<cdouble,1>*, Array<cdouble,1>*);
	// Component Schemes of RK4
	void Single_time_step_EULER(Array<cdouble,1>*, Array<cdouble,1>*, Array<double,1>*, double);
	void Single_time_step_Semi_implicit(Array<cdouble,1>*, Array<cdouble,1>*, Array<double,1>*, double);
	void Single_time_step_RK2(Array<cdouble,1>*, Array<cdouble,1>*, Array<double,1>*, double);
	// Functions for exponential trick
	void Mult_field_exp_ksqr_dt(Array<cdouble,1>*, Array<double,1>*, double);
	void Mult_nlin_exp_ksqr_dt(Array<cdouble,1>*, Array<double,1>*, double);
	// Functions for corrected shells and indices 
	inline int Get_i(int n); 
	inline int Get_n(int i);