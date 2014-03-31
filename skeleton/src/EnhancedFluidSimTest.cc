#include "SORSolver.hh"
#include "EnhancedFluidSimulator.hh"
#include "Debug.hh"
/** Include GrayScaleImage.hh as a custom user defined mesh which can be called from CustomMesh.hh in EnhancedFluidSimulator.cc */
#include "GrayScaleImage.hh"

#include <cmath>
#include <iostream>
#define PI (4.0*atan(1.0))

void ReadFile(FileReader &configu)
{
	configu.registerRealParameter("xlength");
	configu.registerRealParameter("ylength");
	configu.registerIntParameter("imax");
	configu.registerIntParameter("jmax");
	configu.registerStringParameter("SpecProb");
	configu.registerStringParameter("name");
	configu.registerIntParameter("customBC");
	configu.registerIntParameter("customInit");
	configu.registerIntParameter("customMesh");
	configu.registerRealParameter("RectangleX1");
	configu.registerRealParameter("RectangleY1");
	configu.registerRealParameter("RectangleX2");
	configu.registerRealParameter("RectangleY2");
	configu.registerRealParameter("CircleX");
	configu.registerRealParameter("CircleY");
	configu.registerRealParameter("CircleR");
	configu.registerIntParameter("itermax");
	configu.registerRealParameter("eps");
	configu.registerRealParameter("omg");
	configu.registerRealParameter("dt");
	configu.registerRealParameter("tau");
	configu.registerRealParameter("t_end");
	configu.registerIntParameter("steady");
	configu.registerIntParameter("customUDA");
	configu.registerRealParameter("UDA_init");
	configu.registerRealParameter("D_UDA");
	configu.registerIntParameter("boundary_conditionUDA_N");
	configu.registerIntParameter("boundary_conditionUDA_E");
	configu.registerIntParameter("boundary_conditionUDA_S");
	configu.registerIntParameter("boundary_conditionUDA_W");
	configu.registerIntParameter("boundary_conditionUDA_O");
	configu.registerRealParameter("boundary_UDA_N");
	configu.registerRealParameter("boundary_UDA_E");
	configu.registerRealParameter("boundary_UDA_S");
	configu.registerRealParameter("boundary_UDA_W");
	configu.registerRealParameter("boundary_UDA_O");	
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
	bool configread = configu.readFile("EnhancedFluidSimTest.par");
	CHECK_MSG(configread, "Could not open file 'EnhancedFluidSimTest.par' which has to be in the current directory.");
}

int main( )
{
	FileReader configu;
	ReadFile(configu);
	EnhancedFluidSimulator FS(configu);
	PROGRESS("Initialize first grid");
	FS.simulate(configu.getRealParameter("t_end"));
	configu.printParameters();
	return 0;
}
