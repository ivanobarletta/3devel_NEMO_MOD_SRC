This folder contains modified source code
of NEMO model, to include PETSc package. 

this version of modified source code solves
the elliptic equation for Filtered free surface
case on all domain ( the matrix size includes 
land point as well ) and works with jperio 0,1,4
( should work with jperio = 3 too, but has not been
tested yet )

to use PETSc package, nn_solv must be set to 3
