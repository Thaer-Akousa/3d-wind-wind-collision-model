# 3d-wind-wind-collision-model
my finals project

lib.c: functions for calculating the shock and post-shock structure.\\
numerical.c: functions and wrappers for numerical methods.\\
vars.c: data reading functions and value assignment of the initial parameters.\\
main.c: run() function for the calculation of spectra and luminosity and main() funciton.\\
run_bs.sh: script file for the phase cycle.\\

to compile and run you need:\\
data of opacity values.\\
data of the cooling function.\\
spectral grid from a planar-shock model.\\

commands:\\
make bs\\
make run\\

FUTURE NOTES:\\
If you want to include centrifugal forces then rewrite in C++ with objects or atleast creat structures.\\
Learn how to run multiple none-dependent function in parallel. \\
Assign initial parameters (terminal vel, mass loss, numb of grid points etc...) from a simple text file instead of delving into the code.\\
