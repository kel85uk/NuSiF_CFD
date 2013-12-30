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




