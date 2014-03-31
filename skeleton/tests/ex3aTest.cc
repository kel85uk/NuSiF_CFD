#include "SORSolver.hh"
#include "FluidSimulator.hh"
#include "Debug.hh"

#include <cmath>
#include <iostream>
#define PI (4.0*atan(1.0))

void gridinit(FluidSimulator& FS)
{
    FS.initFields();
    FS.computeFG();
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
	 configu.registerRealParameter("dt");
	 configu.registerRealParameter("tau");
	 configu.registerRealParameter("t_end");
	 configu.registerIntParameter("timesteps");
	 configu.registerRealParameter("gx");
	 configu.registerRealParameter("gy");
	 configu.registerRealParameter("gamma");
	 configu.registerRealParameter("Re");
	 configu.registerRealParameter("U_init");
	 configu.registerRealParameter("V_init");
	 configu.registerRealParameter("P_init");	 
     bool configread = configu.readFile("ex3aTest.par");
     CHECK_MSG(configread, "Could not open file 'ex3aTest.par' which has to be in the current directory.");
}

int main( )
{
 	 FileReader configu;
 	 ReadFile(configu);
 	 FluidSimulator FS(configu);
//	 StaggeredGrid grid(configu);
	 PROGRESS("Initialize first grid");
     gridinit(FS);
	 FS.grid().p().print();
     FS.grid().U().print();
     FS.grid().V().print();
	 FS.grid().rhs().print();
     FS.grid().F().print();
     FS.grid().G().print();
	 std::cout << FS.grid().rhs().sumI(FS.grid().NGhost()) << std::endl;
     CHECK(FS.grid().rhs().sumI(FS.grid().NGhost()) < 1e-8);
     CHECK(FS.grid().F()(1,1) = FS.grid().U()(1,1));
     CHECK(FS.grid().G()(1,1) = FS.grid().V()(1,1));
//	 FS.grid().YFU().print();
//	 FS.grid().YFV().print();
//	 CHECK(res);
   return 0;
}
