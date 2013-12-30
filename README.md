NuSiF_CFD v1.0
==============

C++ code for solving Incompressible Navier-Stokes (2D)

Written by kklloh (12/30/2013)

Method similar to the code in the book "Numerical Simulation in Fluid Dynamics" by Martin Griebel et al.

Added features:

->  User-defined enabled custom boundary conditions, initial conditions, and mesh
->  CG solver for the pressure poisson correction (continuity) equation
->  Arbitrary obstacle geometry from a png graphics file

First example:
1. Create a build directory in the NuSiF_CFD root directory, e.g. "/${NuSiF_root}/build"
2. cd to the build directory, e.g. "cd /${NuSiF_root}/build"
3. In the build directory, run cmake with reference to the NuSiF_CFD root directory, e.g. "cmake /${NuSiF_root}"
4. Copy the files customMesh.png, and EnhancedFluidSimTest.par from Results_from_tests directory into the same directory as the compiled executable created inside the build directory, e.g. "cp /${NuSiF_root}/Results_from_tests/readFromPNG/EnhancedFluidSimTest.par /${NuSiF_root}/build/tests"
5. Run the code "./EnhancedFluidSimTest"
6. Post-process the vtk files with Paraview or Visit.

TODO:
->  Precondition CG
->  Algebraic Multigrid
->  Extension to 3D
->  Solve transport equations (temperature, chemical species)
->  Parallelization via PETSc
->  Uncertainty quantification modules via MC methods
->  Data assimilation via EnKF
