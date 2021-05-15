# 3d-wind-wind-collision-model
my finals project

lib.c: functions for calculating the shock and post-shock structure.
numerical.c: functions and wrappers for numerical methods.
vars.c: data reading functions and value assignment of the initial parameters.
main.c: run() function for the calculation of spectra and luminosity and main() funciton.
run_bs.sh: script file for the phase cycle.

to compile and run you need:
data of opacity values.
data of the cooling function
spectral grid from a planar-shock model

make bs
make run
