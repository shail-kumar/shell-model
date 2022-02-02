This program is meant for computing the evolution of velocity fields from a dynamical systems perspective. It employs the shell model approach and RK4 integration scheme with exponential integrators. The program is parallelized with MPI.

Author: Shailendra Kumar Rathor 
Email: skrathor@iitk.ac.in, shailkumar22@gmail.com

This work is licensed under Creative Commons Attribution-NonCommercial 4.0 International License. If you find this work useful, please use it.



RELEVANT INFO
Run the executable named setup. This creates three directories namely binaries, input, and output in the parent directory and populates the directories with relevant files. 

#1 The source code is put in the src/ directory. Command <make> will generate the executable in the binaries/ directory.

#2 Parameters and initial conditions are kept in the input/ directory.

#3 The velocity field, moduli of velocity fields and 'energy vs time' is stored in output/ directory.

#4 The program can be run from the binaries/ directory by using the command <make run> or <make>.

#5 IMPORTANT: Set the environment (in binaries/ directory) before running the program.
