#ifndef __EFLUID_SIMULATOR_H__
#define __EFLUID_SIMULATOR_H__

#include "StaggeredGrid.hh"
#include "VTKWriter.hh"
#include "SORSolver.hh"
#include "CGSolver.hh"
#include "FileReader.hh"
#include "Geometry2D.hh"
#include "Types.hh"

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
      
			void computeFG();
			void composeRHS2D();
			void setUVBC2D();
			void setInteriorBC();
			void initFields();
			void adapUV2D();

  private:
			void setNBC2D(Array<real>& bcVecU,Array<real>& bcVecV, real Uw, int NX, int NY, int nGhost, int type);
			void setEBC2D(Array<real>& bcVecU,Array<real>& bcVecV, real Vw, int NX, int NY, int nGhost, int type);
			void setSBC2D(Array<real>& bcVecU,Array<real>& bcVecV, real Uw, int NX, int NY, int nGhost, int type);
			void setWBC2D(Array<real>& bcVecU,Array<real>& bcVecV, real Vw, int NX, int NY, int nGhost, int type);
			void setNBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type);
			void setEBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type);
			void setSBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type);
			void setWBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type);			
			void setDelta_t2D();
			void normalizePressure();
  			FileReader conf_;
//			SORSolver psolver_;
			CGSolver psolver_;
			Geometry2D mesh_;//,meshvis_;
			Array<real> du2dx, duvdy, d2udx2, d2udy2, dpdx, duvdx, dv2dy, d2vdx2, d2vdy2, dpdy; //Kept in array form since post processing will be easier? can change to real if insufficient memory
			      // Functions to compute derivatives
/*		   inline real dU2_dx(int i, int j);
		   inline real dV2_dy(int i, int j);
		   inline real dUV_dy(int i, int j);
		   inline real dVU_dx(int i, int j);
		   inline real d2U_dx2(int i, int j);
		   inline real d2U_dy2(int i, int j);
		   inline real d2V_dx2(int i, int j);
		   inline real d2V_dy2(int i, int j); */
			int Ntimes_,pnorm_,rescheckfreq_,writeVTKfreq_;
			int steadys_ = 0;
  			StaggeredGrid grid_;			
			std::string casename;
			bool customBC = false, customInit = false, customMesh = false, customUDA = false;
			int specialProblems = 0;
			real gam_;
			real dt_, tau_ = 0.9, t_end;
			real RE_, D_UDA = 0.;
			real gx_;
			real gy_;
			real uinit_, vinit_, pinit_;
            real U_NBC_ = 0., V_NBC_ = 0., U_EBC_ = 0., V_EBC_ = 0., U_SBC_ = 0., V_SBC_ = 0., U_WBC_ = 0., V_WBC_ = 0.;
            real UDA_NBC_ = 0., UDA_EBC_ = 0., UDA_SBC_ = 0., UDA_WBC_ = 0., UDA_OBC_ = 0.;
            int Nbctype_ = 0, Ebctype_ = 0, Sbctype_ = 0, Wbctype_ = 0;
            int NbctypeUDA_ = 0, EbctypeUDA_ = 0, SbctypeUDA_ = 0, WbctypeUDA_ = 0, ObctypeUDA_ = 0;
};



#endif
