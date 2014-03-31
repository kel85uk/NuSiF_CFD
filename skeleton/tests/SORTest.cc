#include "SORSolver.hh"
#include "StaggeredGrid.hh"
#include "Debug.hh"

#include <cmath>
#include <iostream>
#define PI (4.0*atan(1.0))



void initGridSetup1( StaggeredGrid & grid )
{
	grid.p().randfill(-10,10);
	grid.rhs().fill(0);
}

void initGridSetup2( StaggeredGrid & grid )
{
	grid.p().randfill(-10,10);
	for (int i = 0; i < grid.rhs().getSize(0); ++i)
		for (int j = 0; j < grid.rhs().getSize(1); ++j)
			grid.rhs()(i,j) = sin(2.0*grid.XC()(i,j)*PI);
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
	 configu.registerRealParameter("eps");
	 configu.registerRealParameter("omg");
//	 configu.registerStringParameter("NormType");	 
	 bool configread = configu.readFile("poisson.par");
	 CHECK_MSG(configread, "Could not open file 'poisson.par' which has to be in the current directory.");
}

int main()
{
 	 FileReader configu;
 	 ReadFile(configu);
	 StaggeredGrid grid(configu);
	 SORSolver solver(configu);	 

	 initGridSetup1(grid);
	 bool res = solver.solve(grid);
	 CHECK(res);	 	 
	 CHECK(solver.residualnorm() < configu.getRealParameter("eps"));
	 initGridSetup2( grid );
	 res = solver.solve(grid);
	 CHECK(res);
	 CHECK(solver.residualnorm() < configu.getRealParameter("eps"));

   return 0;
}
