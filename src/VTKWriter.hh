#ifndef VTKFILEWRITER_HH
#define VTKFILEWRITER_HH


#include "StaggeredGrid.hh"
#include "Geometry2D.hh"



//*******************************************************************************************************************
/*! Writes StaggeredGrid as vtk file

  - vtk files can for example be opened with Paraview ( http://www.paraview.org/ )
  - writes out pressure and/or velocity

  - Usage:
   \code
       VTKWriter vtkWriter ( myGrid, "lidDrivenCavity", true, true );
       // for each timestep:
      vtkWriter.write();
    \endcode
    This creates on file per timestep: "lidDrivenCavity_0001.vtk", ""lidDrivenCavity_0002.vtk" ...

*/
//*******************************************************************************************************************
class VTKWriter
{

public:
   VTKWriter(  const StaggeredGrid & grid, Geometry2D& mesh, const std::string & basename,
               bool writePressure = true, bool writeVelocity = true, bool writeObs = true, bool writeMeshVTK = false );
   VTKWriter(  const StaggeredGrid & grid, const std::string & basename,
               bool writePressure = true, bool writeVelocity = true);               
	 /** For backward compatibility reasons used in FluidSimulator.cc */
   void write();
   /** For writing custom mesh used in EnhancedFluidSimulator.cc */
   void write(Geometry2D& mesh, real time);

private:
   const StaggeredGrid & grid_;
   std::string baseName_;

   bool writeVelocity_;
   bool writePressure_;
   bool writeObs_;
   bool writeMeshVTK_;

   int counter_;
   std::string header_;

};



#endif




