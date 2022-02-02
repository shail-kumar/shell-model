#~ all: shell
INCLUDE = -I/home/shailen/local/include/
OBJ = shell.o input_functions.o evolution.o processing_functions.o output_functions.o
#~ OBJ = src/*.o
COMPILER = mpicxx
DEBUG = -g
PROFILING = -p
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)
#~ ATTEMPT = 0
#~ num_procs = 3
ATTEMPT ?= $(shell bash -c 'read -p "Target name: " pwd; echo $$pwd')
num_procs ?= $(shell bash -c 'read -p "No of procceses: " pwd; echo $$pwd')
OMP = -fopenmp

target: $(OBJ)
	$(COMPILER) $(DEBUG) $(OBJ) -o ../binaries/shell_$(ATTEMPT) $(INCLUDE)
	
shell.o: shell.h shell.cc
	$(COMPILER) $(CFLAGS) shell.cc $(INCLUDE) 
	
input_functions.o: shell.h input_functions.cc
	$(COMPILER) $(CFLAGS) input_functions.cc $(INCLUDE)
	
evolution.o: evolution.h evolution.cc shell.h
	$(COMPILER) $(CFLAGS) evolution.cc $(INCLUDE)
	
processing_functions.o: shell.h evolution.h processing_functions.cc
	$(COMPILER) $(CFLAGS) processing_functions.cc $(INCLUDE)
	
output_functions.o: shell.h output_functions.cc
	$(COMPILER) $(CFLAGS) output_functions.cc $(INCLUDE)
clean:
	rm ./*.o
	
tar:
	tar cvf shell_$(shell date +"%F_T%H-%M").tar shell.h shell.cc input_functions.cc \
	evolution.h evolution.cc processing_functions.cc output_functions.cc makefile \
	parameters.d initial_field.d environment.d Makefile_binary README.txt license.html \
	setup
	
run:
	@mpirun -np $(num_procs) ./shell_$(ATTEMPT) | tee std_out_3.d
	
flux: flux.o
	$(COMPILER) $(DEBUG) flux.o -o flux $(INCLUDE)
	      
flux.o: flux.cc shell.h
	$(COMPILER) $(CFLAGS) flux.cc $(INCLUDE)

hflux: hflux.o
	$(COMPILER) $(DEBUG) hflux.o -o hflux $(INCLUDE)
	      
hflux.o: hflux.cc shell.h
	$(COMPILER) $(CFLAGS) hflux.cc $(INCLUDE)
	
rf: real_field.o
	$(COMPILER) $(DEBUG) real_field.o -o rf $(INCLUDE)
	
real_field.o: real_field.cc shell.h
	$(COMPILER) $(CFLAGS) real_field.cc $(INCLUDE)

