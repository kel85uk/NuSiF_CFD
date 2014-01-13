/*    Copyright (C) 2013  kklloh

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "PSolver.hh"
#include <MatrixCOO.hh>
#include "StaggeredGrid.hh"
#include "Geometry2D.hh"
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
	 configu.registerIntParameter("imax");
	 configu.registerIntParameter("jmax");
	 configu.registerIntParameter("itermax");
	 configu.registerRealParameter("eps");
	 configu.registerRealParameter("omg");
	 configu.registerIntParameter("Solver");
	 bool configread = configu.readFile("PSolver.par");
	 CHECK_MSG(configread, "Could not open file 'PSolver.par' which has to be in the current directory.");
}

int main()
{
 	 FileReader configu;
 	 ReadFile(configu);
	 StaggeredGrid grid(configu);
	 PSolver solver(configu);	 
	 Geometry2D mesh(grid);
	 mesh.geom_init();
	 initGridSetup1(grid);
	 bool res = solver.solve(grid,mesh,0);
	 CHECK(res);	 	 
//	 CHECK(solver.residualnorm() < configu.getRealParameter("eps"));
//*/
	 
/*	 initGridSetup2( grid );
	 res = solver.solve(grid);
	 CHECK(res);
	 CHECK(solver.residualnorm() < configu.getRealParameter("eps"));
*/
   return 0;
}
