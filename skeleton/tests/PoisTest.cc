#include <iostream>
#include <cmath>
#include "SORSolver.hh"
#include "CGSolver.hh"

// 2D Setup 1:
void initGridSetup1( StaggeredGrid & grid )
{
   Array<real> & rhs = grid.rhs();
   Array<real> & p = grid.p();
   for (int i=0; i<rhs.getSize(0);i++)
       for (int j=0; j<rhs.getSize(1);j++) {
   //    - grid.p   : init with random values
           p(i,j) =  rand() % 10000-20000;
   //    - grid.rhs : init with zero
           rhs(i,j) = 0; }
}

// 2D Setup 2:
void initGridSetup2( StaggeredGrid & grid )
{
   Array<real> & rhs = grid.rhs();
   Array<real> & p = grid.p();
	for (int i = 0; i < grid.rhs().getSize(0); ++i)
		for (int j = 0; j < grid.rhs().getSize(1); ++j){
			p(i,j) =  rand() % 10000-20000;
			rhs(i,j) = sin(2.0*grid.XC()(i,j)*M_PI);
		}

}

void ReadFile(FileReader &configu)
{
 	 configu.registerRealParameter("xlength");
	 configu.registerRealParameter("ylength");
//	 configu.registerIntParameter("nGhost");
	 configu.registerIntParameter("imax");
	 configu.registerIntParameter("jmax");
//	 configu.registerStringParameter("type");
	 configu.registerIntParameter("itermax");
	 configu.registerIntParameter("checkfrequency");
	 configu.registerRealParameter("eps");
	 configu.registerRealParameter("omg");
	 configu.registerStringParameter("solver");	 
	 bool configread = configu.readFile("poisson.par");
	 CHECK_MSG(configread, "Could not open file 'poisson.par' which has to be in the current directory.");
}

int main()
{
   bool converged;
   std::cout<<"\nReading file ";
   FileReader configuration;
   ReadFile(configuration);
   std::cout<<"\nMaking Staggered Grid\n";
   StaggeredGrid grid(configuration);

   std::cout<<"\nCreating Solvers\n";
	CGSolver solverCG(configuration);
	SORSolver solverSOR(configuration);

   std::cout<<"\nInitializing Grid Setup 1\n";
   initGridSetup1( grid );
   std::cout<<"\nRunning Solver for Grid Setup 1\n";
   std::clock_t startcg = 0;
   double durationcg;
   converged = solverCG.solve(grid);
   durationcg = (double)( std::clock() - startcg ) / (double) CLOCKS_PER_SEC;
   std::cout<<"printf: "<< durationcg <<'\n';
  
   CHECK_MSG(converged , "Grid Setup 1 not converged");
   std::cout<<"\nInitializing Grid Setup 2\n";
   initGridSetup2( grid );
   std::cout<<"\nRunning Solver for Grid Setup 2\n";
   converged = solverCG.solve(grid);
   CHECK_MSG(converged, "Grid Setup 2 not converged");
   std::cout<<"\nCG Test 2D passed\n";
   std::cout<<"\nInitializing Grid Setup 1\n";
   initGridSetup1( grid );
   std::cout<<"\nRunning Solver for Grid Setup 1\n";
   std::clock_t startsor = 0;
   double durationsor;
   converged = solverSOR.solve(grid);
   durationsor = (double)( std::clock() - startsor ) / (double) CLOCKS_PER_SEC;
   std::cout<<"printf: "<< durationsor <<'\n';
  
   CHECK_MSG(converged , "Grid Setup 1 not converged");
   std::cout<<"\nInitializing Grid Setup 2\n";
   initGridSetup2( grid );
   std::cout<<"\nRunning Solver for Grid Setup 2\n";
   converged = solverSOR.solve(grid);
   CHECK_MSG(converged, "Grid Setup 2 not converged");
   std::cout<<"\nSOR Test 2D passed\n";   
   return 0;
}
