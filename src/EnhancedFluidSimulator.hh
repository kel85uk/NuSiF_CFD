#ifndef __EFLUID_SIMULATOR_H__
#define __EFLUID_SIMULATOR_H__

#include "StaggeredGrid.hh"
#include "VTKWriter.hh"
#include "SORSolver.hh"
#include "FileReader.hh"
#include "Geometry2D.hh"

class EnhancedFluidSimulator
{
  public:
  			EnhancedFluidSimulator();
      EnhancedFluidSimulator( const FileReader & conf );

      /// Simulates a given time-length
      void simulate             ( real duration              );
      void simulateTimeStepCount( unsigned int nrOfTimeSteps );


      // Getter functions for the internally stored StaggeredGrid
            StaggeredGrid & grid();
      const StaggeredGrid & grid() const;
      // Getter function for the geometry for calc
      				Geometry2D & MESH();
      	const Geometry2D & MESH() const;
      // Getter function for the geometry for visualization
      				Geometry2D & MESHVIS();
      	const Geometry2D & MESHVIS() const;      	
      
			void computeFG();
			void composeRHS2D();
			void setUVBC2D();
			void initFields();
			void adapUV2D();

  private:
			void setNBC2D(Array& bcVecU,Array& bcVecV, real Uw, int NX, int NY, int nGhost, int type);
			void setEBC2D(Array& bcVecU,Array& bcVecV, real Vw, int NX, int NY, int nGhost, int type);
			void setSBC2D(Array& bcVecU,Array& bcVecV, real Uw, int NX, int NY, int nGhost, int type);
			void setWBC2D(Array& bcVecU,Array& bcVecV, real Vw, int NX, int NY, int nGhost, int type);
			void setNBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type);
			void setEBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type);
			void setSBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type);
			void setWBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type);			
			void setDelta_t2D();
			void normalizePressure();
  			StaggeredGrid grid_;
  			FileReader conf_;
			SORSolver psolver_;
			Geometry2D mesh_,meshvis_;
			Array du2dx, duvdy, d2udx2, d2udy2, dpdx, duvdx, dv2dy, d2vdx2, d2vdy2, dpdy; //Kept in array form since post processing will be easier? can change to real if insufficient memory
			int Ntimes_,pnorm_,rescheckfreq_,writeVTKfreq_;
			std::string casename;
			bool customBC = false, customInit = false, customMesh = false;
			int specialProblems = 0;
			real gam_;
			real dt_, tau_ = 0.9, t_end;
			real RE_;
			real gx_;
			real gy_;
			real uinit_, vinit_, pinit_;
            real U_NBC_, V_NBC_, U_EBC_, V_EBC_, U_SBC_, V_SBC_, U_WBC_, V_WBC_ = 0.;
            int Nbctype_, Ebctype_, Sbctype_, Wbctype_ = 0;
};



#endif
