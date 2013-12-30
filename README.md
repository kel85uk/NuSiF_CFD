NuSiF_CFD v1.0
==============

C++ code for solving Incompressible Navier-Stokes (2D)

Written by kklloh (12/30/2013)

Method similar to the code in the book "Numerical Simulation in Fluid Dynamics" by Martin Griebel et al.

Added features:

->  User-defined enabled custom boundary conditions, initial conditions, and mesh
->  CG solver for the pressure poisson correction (continuity) equation
->  Arbitrary obstacle geometry from a png graphics file

TODO:
->  Precondition CG
->  Algebraic Multigrid
->  Extension to 3D
->  Solve transport equations (temperature, chemical species)
->  Parallelization via PETSc
->  Uncertainty quantification modules via MC methods
->  Data assimilation via EnKF
