#include "VTKWriter.hh"
#include "Debug.hh"

#include <fstream>
#include <sstream>
#include <iomanip>


template<typename T> struct RealTypeToString         {};
template<>           struct RealTypeToString<float>  { static const char * str; };
template<>           struct RealTypeToString<double> { static const char * str; };

const char * RealTypeToString<float>::str  = "float";
const char * RealTypeToString<double>::str = "double";



VTKWriter::VTKWriter(  const StaggeredGrid & grid, Geometry2D& mesh, const std::string & basename, bool writePressure, bool writeVelocity, bool writeObs, bool writeMeshVTK )
      : grid_(grid),baseName_( basename ),
        writeVelocity_(writeVelocity), writePressure_(writePressure), writeObs_(writeObs), writeMeshVTK_(writeMeshVTK), counter_ (0 )
{
   ASSERT_MSG( writePressure_ || (writeVelocity_ || writeObs_) , "VTK Writer has to write at least velocity, pressure, or obstacle_mesh" );

   std::stringstream sstream;
   sstream << "# vtk DataFile Version 4.0\n";

/*   sstream << "DIMENSIONS " << grid_.Nx() << " " << grid_.Ny() << " 1\n";
   sstream << "ORIGIN 0 0 0 \n";
   sstream << "SPACING " << grid_.dx() << " " << grid_.dy() << " 1\n";
   sstream << "POINT_DATA " << grid_.Nx() * grid_.Ny() << " \n" << std::endl;  */ 
//   sstream << "ASCII\n";
//   sstream << "DATASET STRUCTURED_GRID\n";
//   sstream << "DIMENSIONS " << grid_.Nx() << " " << grid_.Ny() << " 1\n";
   if (writeMeshVTK_){
   	 sstream << "POINTS " << mesh.geom_nFluids() << " " << RealTypeToString<real>::str << "\n";
		 for ( int j = 0; j < grid_.Ny (); ++j )
		 	for ( int i = 0; i < grid_.Nx (); ++i )
		 		if(mesh.geom_isFluid(i+1,j+1))
			 		sstream << grid_.XC()(i+1,j+1) << " " << grid_.YC()(i+1,j+1) << " " << " 0\n";
//		 sstream << "POINT_DATA " << mesh.geom_nFluids() << " \n" << std::endl;
	 }
	 else{
/*	   sstream << "POINTS " << grid_.Nx()*grid_.Ny() << " " << RealTypeToString<real>::str << "\n";	 
		 for ( int j = 0; j < grid_.Ny (); ++j )
		 	for ( int i = 0; i < grid_.Nx (); ++i )
			 		sstream << grid_.XC()(i+1,j+1) << " " << grid_.YC()(i+1,j+1) << " " << " 0\n";
		 sstream << "POINT_DATA " << grid_.Nx()*grid_.Ny() << " \n" << std::endl; */
	 }
   header_ = sstream.str();
}

VTKWriter::VTKWriter(  const StaggeredGrid & grid, const std::string & basename, bool writePressure, bool writeVelocity )
      : grid_(grid), baseName_( basename ),
        writeVelocity_(writeVelocity), writePressure_(writePressure), counter_ (0 )
{
   ASSERT_MSG( writePressure_ || (writeVelocity_ || writeObs_) , "VTK Writer has to write at least velocity, pressure, or obstacle_mesh" );

   std::stringstream sstream;
   sstream << "# vtk DataFile Version 4.0\n";
   sstream << "Nusif VTK output\n";
   sstream << "ASCII\n";
   sstream << "DATASET STRUCTURED_POINTS\n";
   sstream << "DIMENSIONS " << grid_.Nx() << " " << grid_.Ny() << " 1\n";
   sstream << "ORIGIN 0 0 0 \n";
   sstream << "SPACING " << grid_.dx() << " " << grid_.dy() << " 1\n";
   sstream << "POINT_DATA " << grid_.Nx() * grid_.Ny() << " \n" << std::endl; 

   header_ = sstream.str();
}

void VTKWriter::write(Geometry2D& mesh, real time)
{
   std::stringstream fileName;
   fileName << baseName_ << "_" <<  std::setw(4) << std::setfill( '0') << counter_ << ".vtk";
   std::ofstream fileStream ( fileName.str().c_str() );
   fileStream << header_;
   fileStream << "Nusif VTK output\n";
   fileStream << "ASCII\n";
   fileStream << "DATASET STRUCTURED_POINTS\n";   
/*   fileStream << "FIELD FieldData 1\n";
   fileStream << "TIME 1 1 " << RealTypeToString<real>::str << "\n";
   fileStream << time << "\n\n"; */
   fileStream << "DIMENSIONS " << grid_.Nx() << " " << grid_.Ny() << " 1\n";
   fileStream << "ORIGIN 0 0 0 \n";
   fileStream << "SPACING " << grid_.dx() << " " << grid_.dy() << " 1\n";
   fileStream << "POINT_DATA " << grid_.Nx() * grid_.Ny() << " \n" << std::endl;   
	 if (!writeMeshVTK_){
		 if ( writeVelocity_ )
		 {
		    fileStream << "VECTORS velocity " << RealTypeToString<real>::str << "\n";

            for ( int j = 1; j <= grid_.Ny (); ++j )
               for ( int i = 1; i <= grid_.Nx (); ++i )
               {
                  const real u = 0.5 * ( grid_.U() ( i, j )  + grid_.U() ( i-1, j ) );
                  const real v = 0.5 * ( grid_.V() ( i, j )  + grid_.V() ( i, j-1 ) );
                  fileStream << u << " " << v << " " << " 0\n";
               }

		    fileStream << "\n";
		 }

		 if ( writePressure_ )
		 {
		    fileStream << "SCALARS pressure " << RealTypeToString<real>::str << " 1\n";
		    fileStream << "LOOKUP_TABLE default\n";

		    for ( int j = 0; j < grid_.Ny (); ++j )
		       for ( int i = 0; i < grid_.Nx (); ++i )
		          fileStream << grid_.p()( i+1, j+1 ) << "\n";
		 }
		 
		 if (writeObs_)
		 {
//		 		fileStream << "CELL_DATA " << 9801 << " \n" << std::endl;
		    fileStream << "SCALARS geom " << "int" << " 1\n";
		    fileStream << "LOOKUP_TABLE default\n";

		    for ( int j = 0; j < grid_.Ny (); ++j )
		       for ( int i = 0; i < grid_.Nx (); ++i )
		       	if(mesh.geom_isFluid(i+1,j+1))
		          fileStream << 1 << "\n";
		        else fileStream << 0 << "\n";
		 }
			fileStream << "SCALARS UDA " << RealTypeToString<real>::str << " 1\n";
			fileStream << "LOOKUP_TABLE default\n";
		    for ( int j = 0; j < grid_.Ny (); ++j )
		       for ( int i = 0; i < grid_.Nx (); ++i )
		          fileStream << grid_.UDA()( i+1, j+1 ) << "\n";
	 }
   ++counter_;
}

void VTKWriter::write(){
	std::stringstream fileName;
	fileName << baseName_ << "_" <<  std::setw(4) << std::setfill( '0') << counter_ << ".vtk";
	std::ofstream fileStream ( fileName.str().c_str() );

	fileStream << header_;
	 if ( writeVelocity_ )
	 {
		  fileStream << "VECTORS velocity " << RealTypeToString<real>::str << "\n";

          for ( int j = 1; j <= grid_.Ny (); ++j )
             for ( int i = 1; i <= grid_.Nx (); ++i )
		     {
                const real u = 0.5 * ( grid_.U() ( i, j )  + grid_.U() ( i-1, j ) );
                const real v = 0.5 * ( grid_.V() ( i, j )  + grid_.V() ( i, j-1 ) );
		        fileStream << u << " " << v << " " << " 0\n";
		     }

		  fileStream << "\n";
	 }

	 if ( writePressure_ )
	 {
		  fileStream << "SCALARS pressure " << RealTypeToString<real>::str << " 1\n";
		  fileStream << "LOOKUP_TABLE default\n";

		  for ( int j = 0; j < grid_.Ny (); ++j )
		     for ( int i = 0; i < grid_.Nx (); ++i )
		        fileStream << grid_.p()( i+1, j+1 ) << "\n";
	 }
   ++counter_;
}
