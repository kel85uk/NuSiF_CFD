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

#include "SORSolver.hh"
#include "FluidSimulator.hh"
//#include "StaggeredGrid.hh"
#include "Debug.hh"

#include <cmath>
#include <iostream>
#define PI (4.0*atan(1.0))

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
	 configu.registerRealParameter("dt");
	 configu.registerRealParameter("tau");
	 configu.registerRealParameter("t_end");
	 configu.registerIntParameter("timesteps");
	 configu.registerIntParameter("checkfrequency");
	 configu.registerIntParameter("normalizationfrequency");
	 configu.registerIntParameter("writefrequency");
	 configu.registerRealParameter("gx");
	 configu.registerRealParameter("gy");
	 configu.registerRealParameter("gamma");
	 configu.registerRealParameter("Re");
	 configu.registerRealParameter("U_init");
	 configu.registerRealParameter("V_init");
   configu.registerRealParameter("P_init");
	 configu.registerIntParameter("boundary_condition_N");
   configu.registerIntParameter("boundary_condition_E");
   configu.registerIntParameter("boundary_condition_S");
   configu.registerIntParameter("boundary_condition_W");
   configu.registerRealParameter("boundary_velocityU_N");
   configu.registerRealParameter("boundary_velocityV_N");
   configu.registerRealParameter("boundary_velocityU_E");
   configu.registerRealParameter("boundary_velocityV_E");
   configu.registerRealParameter("boundary_velocityU_S");
   configu.registerRealParameter("boundary_velocityV_S");
   configu.registerRealParameter("boundary_velocityU_W");
   configu.registerRealParameter("boundary_velocityV_W");
   bool configread = configu.readFile("dcavity.par");
   CHECK_MSG(configread, "Could not open file 'dcavity.par' which has to be in the current directory.");
}

int main( )
{
 	 FileReader configu;
 	 ReadFile(configu);
     configu.printParameters();
 	 FluidSimulator FS(configu);
//	 StaggeredGrid grid(configu);
	 PROGRESS("Initialize first grid");
     FS.simulate(configu.getRealParameter("t_end"));
	 FS.grid().p().print();
     FS.grid().U().print();
     FS.grid().V().print();
	 FS.grid().rhs().print();
	 std::cout << FS.grid().rhs().sumI(FS.grid().NGhost()) << std::endl;
     configu.printParameters();
   return 0;
}
