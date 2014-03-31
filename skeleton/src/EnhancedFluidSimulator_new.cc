#include "EnhancedFluidSimulator.hh"
#include "Array_new.hh"
#include "GrayScaleImage.hh"
#include <cmath>

EnhancedFluidSimulator::EnhancedFluidSimulator(){
}

EnhancedFluidSimulator::EnhancedFluidSimulator( const FileReader & conf ):steadys_(conf.getIntParameter("steady")),grid_(StaggeredGrid(conf)),gamma(conf.getRealParameter("gamma")), Re_1(1.0/conf.getRealParameter("Re")), dx(grid_.dx()), dy(grid_.dy()), dx_2(1.0/(dx*dx)), dy_2(1.0/(dy*dy)),inv4dx(0.25/dx), inv4dy(0.25/dy){
	StaggeredGrid domain(conf);
	SORSolver pdomain(conf);
	Geometry2D mesh(domain);
//	grid_ = domain;
	conf_ = conf;
	psolver_ = pdomain;
	casename = conf.getStringParameter("name");
	std::string specprob1 = "backwardstep";
	std::string specprob2 = "cylinder2D";
	if(conf.getStringParameter("SpecProb") == specprob1)
		specialProblems = 1;
	else if(conf.getStringParameter("SpecProb") == specprob2)
		specialProblems = 2;	
	if(conf.getIntParameter("customMesh") != 0){ customMesh = true;}	
	mesh_ = mesh;
	mesh_.geom_init(conf,grid_);
	if (customMesh){
		std::cout << "Using custom geometry provided by the user\n\n";
		#include "CustomMesh.hh"
	}
	//	meshvis_ = mesh_;
	mesh_.geom_finalize();
	for (int i = 0; i < grid_.Nx() + 2; ++i)
		for (int j = 0; j < grid_.Ny() + 2; ++j){
			grid_.Geom()(i,j) = (mesh_.geom_isFluid(i,j)||(i != 0) || (i!=grid_.Nx() +1) || (j!=0) || (j!=grid_.Ny() +1))? 1 : 0;
		}
	dt_ = conf.getRealParameter("dt");
	tau_ = conf.getRealParameter("tau");
	t_end = conf.getRealParameter("t_end");
	Ntimes_ = conf.getIntParameter("timesteps");
	gx_ = conf.getRealParameter("gx");
	gy_ = conf.getRealParameter("gy");
	gam_ = conf.getRealParameter("gamma");
	RE_ = conf.getRealParameter("Re");
	uinit_ = conf.getRealParameter("U_init");
	vinit_ = conf.getRealParameter("V_init");
	pinit_ = conf.getRealParameter("P_init");
	if(conf.getIntParameter("customInit") != 0){ customInit = true;}
	if(conf.getIntParameter("customBC") != 0){ customBC = true; std::cout << "Using custom boundary conditions provided by the user\n\n";}
	pnorm_ = conf.getIntParameter("normalizationfrequency");
	writeVTKfreq_ = conf.getIntParameter("writefrequency");
	Nbctype_ = conf.getIntParameter("boundary_condition_N");
	Ebctype_ = conf.getIntParameter("boundary_condition_E");
	Sbctype_ = conf.getIntParameter("boundary_condition_S");
	Wbctype_ = conf.getIntParameter("boundary_condition_W");
	U_NBC_ = conf.getRealParameter("boundary_velocityU_N");
	V_NBC_ = conf.getRealParameter("boundary_velocityV_N");
	U_EBC_ = conf.getRealParameter("boundary_velocityU_E");
	V_EBC_ = conf.getRealParameter("boundary_velocityV_E");
	U_SBC_ = conf.getRealParameter("boundary_velocityU_S");
	V_SBC_ = conf.getRealParameter("boundary_velocityV_S");
	U_WBC_ = conf.getRealParameter("boundary_velocityU_W");
	V_WBC_ = conf.getRealParameter("boundary_velocityV_W");
}

/** Simulates a given duration with variable dt within the CFL limit for stability */
void EnhancedFluidSimulator::simulate             ( real duration              ){
	EnhancedFluidSimulator::initFields();
	int i = 0;
	real time = 0.;
	duration = t_end;
	VTKWriter vtkWriter ( grid_,mesh_, casename, true, true, true, false );
//	VTKWriter vtkWriteMesh (grid_,mesh_,"lidDrivenCavity_obs",true,true,true,true);    	
  while ((time < duration)&&(i<Ntimes_)){
		EnhancedFluidSimulator::setDelta_t2D();  
		EnhancedFluidSimulator::setUVBC2D();
		EnhancedFluidSimulator::computeFG();
		EnhancedFluidSimulator::composeRHS2D();
		psolver_.solve(grid_,mesh_);
		EnhancedFluidSimulator::adapUV2D();
		if(i%pnorm_ == 0){
			EnhancedFluidSimulator::normalizePressure();
		}
		if(i%writeVTKfreq_ == 0){
			vtkWriter.write(mesh_,time);
			std::cout << "Writing to file\n";
		}		
		time += dt_;
		++i;
		std::cout << "Timestep: " << i << " Time: " << time << " Delta_t: " << dt_ << " Pres. Iterations: " << psolver_.iterations() << " Pres. residual: " << psolver_.residualnorm() << " : " << grid_.rhs().sum() << std::endl;			
	}
	vtkWriter.write(mesh_,time);
	PROGRESS("Finished time simulations\n");

}
/** Simulates for a pre-determined number of time steps (Usually for testing new algorithms) */
void EnhancedFluidSimulator::simulateTimeStepCount( unsigned int nrOfTimeSteps ){
	EnhancedFluidSimulator::initFields();
	nrOfTimeSteps = Ntimes_;
	VTKWriter vtkWriter ( grid_,mesh_, casename, true, true, true, false );
	real time = 0.;
	for (unsigned int i = 0; i < nrOfTimeSteps; ++i){
		EnhancedFluidSimulator::setUVBC2D();
		EnhancedFluidSimulator::computeFG();
		EnhancedFluidSimulator::composeRHS2D();
		psolver_.solve(grid_,mesh_);
		EnhancedFluidSimulator::adapUV2D();
		if(i%pnorm_ == 0){
			EnhancedFluidSimulator::normalizePressure();
		}
		if(i%writeVTKfreq_ == 0){
			vtkWriter.write(mesh_,time);
			std::cout << "Writing to file\n";
		}
		time += dt_;
		std::cout << "Timestep: " << i << " Time: " << time << " Delta_t: " << dt_ << " Pres. Iterations: " << psolver_.iterations() << " Pres. residual: " << psolver_.residualnorm() << " : " << grid_.rhs().sum() << std::endl;			
	}
	vtkWriter.write(mesh_,time);
	PROGRESS("Finished time simulations\n");
}


/** Getter functions for the internally stored StaggeredGrid */
StaggeredGrid & EnhancedFluidSimulator::grid(){
	return grid_;
}
const StaggeredGrid & EnhancedFluidSimulator::grid() const{
	return grid_;
}

/** Getter function for the internally stored Geometry (for computation) with read-write */
Geometry2D & EnhancedFluidSimulator::MESH(){
	return mesh_;
}
/** Getter function for the internally stored Geometry (for computation) with read-only */
const Geometry2D & EnhancedFluidSimulator::MESH() const{
	return mesh_;
}

void EnhancedFluidSimulator::initFields(){
	real pinit = pinit_, uinit = uinit_, vinit = vinit_;
	grid_.p().fill(0);
	grid_.U().fill(0);
	grid_.V().fill(0);
	for (int i = 1; i <= grid_.Nx(); ++i)
		for (int j = 1; j <= grid_.Ny(); ++j){
			if(mesh_.geom_isFluid(i,j)){
				grid_.p()(i,j) = pinit;
				grid_.U()(i,j) = uinit;
				grid_.V()(i,j) = vinit;
			}
		}
	if(customInit){
		std::cout << "Using custom initial conditions provided by the user\n\n";
		#include "CustomInit.hh"
	}
}
/** Normalizes the pressure based on mean pressure, E[p] */
void EnhancedFluidSimulator::normalizePressure(){
	real p_mean = grid_.p().sum()/mesh_.geom_nFluids();//((grid_.Nx()+2)*(grid_.Ny()+2));
	for (int i = 0; i < grid_.p().getSize(); ++i){
		grid_.p()(i) = grid_.p()(i) - p_mean;
	}
	std::cout << "Normalizing pressure\n";
}

/** Calculates the initial guess of the velocities without pressure correction */
void EnhancedFluidSimulator::computeFG(){
//	Array<real> &p = grid_.p();
	Array<real> &U = grid_.U();
	Array<real> &V = grid_.V();
	Array<real> &F = grid_.F();
	Array<real> &G = grid_.G();
//	duvdy = du2dx = d2udx2 = d2udy2 = dpdx = duvdx = dv2dy = d2vdx2 = d2vdy2 = dpdy = p;	
	real dt = dt_;
//	real RE = RE_;
	real gx = gx_;
	real gy = gy_;
	int imax = grid_.Nx(),jmax = grid_.Ny();
   for (int i=1; i<imax; i++) 
       for (int j=1; j<=jmax; j++)
           F(i,j) = grid_.u(i,j,CENTER) + dt *  ( Re_1*(d2U_dx2(i,j)+d2U_dy2(i,j)) - dU2_dx(i,j) - dUV_dy(i,j) + gx);

   for (int i=1; i<=imax; i++) 
       for (int j=1; j<jmax; j++)
           G(i,j) = grid_.v(i,j,CENTER) + dt *  ( Re_1*(d2V_dx2(i,j)+d2V_dy2(i,j)) - dV2_dy(i,j) - dVU_dx(i,j) + gy);

   for (int i=0; i<=imax; i++){
       G(i,0) = V(i,0);
       G(i,jmax) = V(i,jmax); }        

   for (int j=0; j<=jmax; j++) {
       F(0,j) = U(0,j);
       F(imax,j) = U(imax,j); }  
}

void EnhancedFluidSimulator::composeRHS2D(){
	real dt = dt_;
	int imax = grid_.Nx(),jmax = grid_.Ny();
//	PROGRESS("Allocating RHS");
   for (int i=1; i<=imax; i++) 
       for (int j=1; j<=jmax; j++)
           if (mesh_.geom_isFluid(i,j))
              grid_.rhs()(i,j) = (1.0/dt) * ( (grid_.f(i,j,CENTER) - grid_.f(i,j,WEST))/dx + (grid_.g(i,j,CENTER) - grid_.g(i,j,SOUTH))/dy );
/*	real rhs_mean = grid_.rhs().sum()/mesh_.geom_nFluids();
	for (int i = 0; i < grid_.rhs().getSize(); ++i){
		grid_.rhs()(i) = grid_.rhs()(i) - rhs_mean;
	} */
}

void EnhancedFluidSimulator::adapUV2D(){
//	setDelta_t2D();
	real dt = dt_;
	for (int i = 1; i <= grid_.Nx()-1; ++i)
		for (int j = 1; j <= grid_.Ny(); ++j)
			if( mesh_.geom_isFluid(i,j) && mesh_.geom_isFluid(i+1,j) )
				grid_.U()(i,j) = grid_.F()(i,j) - dt/dx*(grid_.p()(i+1,j) - grid_.p()(i,j));
	for (int i = 1; i <= grid_.Nx(); ++i)
		for (int j = 1; j<= grid_.Ny()-1; ++j)
			if(mesh_.geom_isFluid(i,j) && mesh_.geom_isFluid(i,j+1) )	
				grid_.V()(i,j) = grid_.G()(i,j) - dt/dy*(grid_.p()(i,j+1) - grid_.p()(i,j));
}

void EnhancedFluidSimulator::setDelta_t2D(){
	real tau = tau_;
	if (tau > 0){
	dt_ = std::min(RE_/(2*(1/pow(grid_.dx(),2) + 1/pow(grid_.dy(),2))),grid_.dx()/grid_.U().maximum());
	dt_ = tau*std::min(dt_,grid_.dy()/grid_.V().maximum());
	}
}

void EnhancedFluidSimulator::setUVBC2D(){
  real U_NBC = U_NBC_,V_NBC = V_NBC_, U_EBC = U_EBC_,V_EBC = V_EBC_, U_SBC = U_SBC_,V_SBC = V_SBC_, U_WBC= U_WBC_,V_WBC = V_WBC_;
	int NGhost = grid_.NGhost();
	int nx = grid_.Nx();
	int ny = grid_.Ny();
	// Sets custom boundary conditions. User can add profiles (boundary layer?)
	if (customBC){
		Array<real> U_NBCv(nx + 2*NGhost);Array<real> V_NBCv(nx + 2*NGhost);
		U_NBCv.fill(U_NBC_); V_NBCv.fill(V_NBC_);
		Array<real> U_EBCv(ny + 2*NGhost);Array<real> V_EBCv(ny + 2*NGhost);
		U_EBCv.fill(U_EBC_); V_EBCv.fill(V_EBC_);
		Array<real> U_SBCv(nx + 2*NGhost);Array<real> V_SBCv(nx + 2*NGhost);
		U_SBCv.fill(U_SBC_); V_SBCv.fill(V_SBC_);
		Array<real> U_WBCv(ny + 2*NGhost);Array<real> V_WBCv(ny + 2*NGhost);
		U_WBCv.fill(U_WBC_); V_WBCv.fill(V_WBC_);
		#include "CustomBC.hh" // Input for user to modify boundary conditions	
		EnhancedFluidSimulator::setNBC2D(U_NBCv,V_NBCv,U_NBC,nx,ny,NGhost,Nbctype_); //Types 1 = Symmetry, 2 = outflow, 3 = inlet, 4 = periodic, Else = Wall with Uw
		EnhancedFluidSimulator::setEBC2D(U_EBCv,V_EBCv,V_EBC,nx,ny,NGhost,Ebctype_);
		EnhancedFluidSimulator::setSBC2D(U_SBCv,V_SBCv,U_SBC,nx,ny,NGhost,Sbctype_);
		EnhancedFluidSimulator::setWBC2D(U_WBCv,V_WBCv,V_WBC,nx,ny,NGhost,Wbctype_);
	}
	else{
		EnhancedFluidSimulator::setNBC2D(U_NBC,V_NBC,nx,ny,NGhost,Nbctype_); //Types 1 = Symmetry, 2 = outflow, 3 = inlet, 4 = periodic, Else = Wall with Uw
		EnhancedFluidSimulator::setEBC2D(U_EBC,V_EBC,nx,ny,NGhost,Ebctype_);
		EnhancedFluidSimulator::setSBC2D(U_SBC,V_SBC,nx,ny,NGhost,Sbctype_);
		EnhancedFluidSimulator::setWBC2D(U_WBC,V_WBC,nx,ny,NGhost,Wbctype_);
  }
  
}

void EnhancedFluidSimulator::setNBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type){
    if(type == 1){
        for (int i = nGhost; i < NX + nGhost; ++i){
            grid_.U()(i,NY+1) = grid_.U()(i,NY);
            grid_.V()(i,NY) = 0;
        }
//        PROGRESS("Set North as symmetry");
    }
    else if(type == 2){
        for (int i = nGhost; i < NX + nGhost; ++i){
            grid_.U()(i,NY+1) = grid_.U()(i,NY);
            grid_.V()(i,NY) = grid_.V()(NY-1,1);
        }
//        PROGRESS("Set North as outflow");
    }
    else if(type == 3){
        for (int i = nGhost; i < NX + nGhost; ++i){
            grid_.U()(i,NY+1) = -grid_.U()(i,NY) + 2*Uw;
            grid_.V()(i,NY) = Vw;
        }
//        PROGRESS("Set North as inlet");
    }
    else if(type == 4){
        for (int i = nGhost; i < NX + nGhost; ++i){
            grid_.U()(i,NY) = grid_.U()(i,1);
            grid_.V()(i,NY+1) = grid_.V()(i,2);
        }
//        PROGRESS("Set North as periodic");
    }    
    else{
        for (int i = nGhost; i < NX + nGhost; ++i){
            grid_.V()(i,NY) = 0.;
            grid_.U()(i,NY+1) = -grid_.U()(i,NY) + 2*Uw;
        }
//        PROGRESS("Set North as wall");
    }
}

void EnhancedFluidSimulator::setEBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type){
    if(type == 1){
        for (int j = nGhost; j < NY + nGhost; ++j){
            grid_.U()(NX,j) = 0;//-grid_.U()(NX-1,j);
            grid_.V()(NX+1,j) = grid_.V()(NX,j);
        }
//        PROGRESS("Set East as symmetry");
    }
    else if(type == 2){
        for (int j = nGhost; j < NY + nGhost; ++j){
            grid_.U()(NX,j) = grid_.U()(NX-1,j);
            grid_.V()(NX+1,j) = grid_.V()(NX,j);
        }
//        PROGRESS("Set East as outflow");
    }
    else if(type == 3){
        for (int j = nGhost; j < NY + nGhost; ++j){
            grid_.U()(NX,j) = Uw;
            grid_.V()(NX+1,j) = -grid_.V()(NX,j) + 2*Vw;
        }
//        PROGRESS("Set East as inlet");
    }
    else if(type == 4){
        for (int j = nGhost; j < NY + nGhost; ++j){
            grid_.U()(NX,j) = grid_.U()(1,j);
            grid_.V()(NX+1,j) = grid_.V()(2,j);
        }
//        PROGRESS("Set East as periodic");
    }    
    else{
        for (int j = nGhost; j < NY + nGhost; ++j){
            grid_.V()(NX+1,j) = -grid_.V()(NX,j) + 2*Vw;
            grid_.U()(NX,j) = 0.;
        }
//        PROGRESS("Set East as wall");
    }
}

void EnhancedFluidSimulator::setSBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type){

    if(type == 1){
        for (int i = nGhost; i < NX + nGhost; ++i){
            grid_.U()(i,0) = grid_.U()(i,1);
            grid_.V()(i,0) = 0; //-grid_.V()(i,1);
        }
//        PROGRESS("Set South as symmetry");
    }
    else if(type == 2){
        for (int i = nGhost; i < NX + nGhost; ++i){
            grid_.U()(i,0) = grid_.U()(i,1);
            grid_.V()(i,0) = grid_.V()(i,1);
        }
//        PROGRESS("Set South as outflow");
    }
    else if(type == 3){
        for (int i = nGhost; i < NX + nGhost; ++i){
            grid_.U()(i,0) = -grid_.U()(i,1) + 2*Uw;
            grid_.V()(i,0) = Vw;
        }
//        PROGRESS("Set South as inlet");
    }
    else if(type == 4){
        for (int i = nGhost; i < NX + nGhost; ++i){
            grid_.V()(i,0) = grid_.V()(i,NY-1);
            grid_.U()(i,0) = grid_.U()(i,NY-1);
            grid_.U()(i,1) = grid_.U()(i,NY);
        }
//        PROGRESS("Set South as periodic");
    }     
    else{
        for (int i = nGhost; i < NX + nGhost; ++i){
            grid_.V()(i,0) = 0.;
            grid_.U()(i,0) = -grid_.U()(i,1) + 2*Uw;
        }
//        PROGRESS("Set South as wall");
    }
}

void EnhancedFluidSimulator::setWBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type){
    if(type == 1){
        for (int j = nGhost; j < NY + nGhost; ++j){
            grid_.U()(0,j) = 0; //-grid_.U()(1,j);
            grid_.V()(0,j) = grid_.V()(1,j);
        }
//        PROGRESS("Set West as symmetry");
    }
    else if(type == 2){
        for (int j = nGhost; j < NY + nGhost; ++j){
            grid_.U()(0,j) = grid_.U()(1,j);
            grid_.V()(0,j) = grid_.V()(1,j);
        }
//        PROGRESS("Set West as outflow");
    }
    else if(type == 3){
        for (int j = nGhost; j < NY + nGhost; ++j){
            grid_.U()(0,j) = Uw;
            grid_.V()(0,j) = -grid_.V()(1,j) + 2*Vw;
        }
//        PROGRESS("Set West as inlet");
    }
    else if(type == 4){
        for (int j = nGhost; j < NY + nGhost; ++j){
            grid_.U()(0,j) = grid_.U()(NX-1,j);
            grid_.V()(0,j) = grid_.V()(NX-1,j);
            grid_.V()(1,j) = grid_.V()(NX,j);
        }
//        PROGRESS("Set West as periodic");
    }    
    else{
        for (int j = nGhost; j < NY + nGhost; ++j){
            grid_.V()(0,j) = -grid_.V()(1,j) + 2*Vw;
            grid_.U()(0,j) = 0.;
        }
//        PROGRESS("Set West as wall");
    }
}

void EnhancedFluidSimulator::setNBC2D(Array<real>& bcVecU,Array<real>& bcVecV, real Uw, int NX, int NY, int nGhost, int type){
	CHECK_MSG((bcVecU.getSize() == NX + 2*nGhost)&&(bcVecV.getSize() == NX + 2*nGhost),"BC Array not set correctly (NORTH)");
  if(type == 1){
      for (int i = nGhost; i < NX + nGhost; ++i){
          grid_.U()(i,NY+1) = grid_.U()(i,NY);
          grid_.V()(i,NY) = 0;
      }
//        PROGRESS("Set North as symmetry");
  }
  else if(type == 2){
      for (int i = nGhost; i < NX + nGhost; ++i){
          grid_.U()(i,NY+1) = grid_.U()(i,NY);
          grid_.V()(i,NY) = grid_.V()(NY-1,1);
      }
//        PROGRESS("Set North as outflow");
  }
  else if(type == 3){
      for (int i = nGhost; i < NX + nGhost; ++i){
          grid_.U()(i,NY+1) = -grid_.U()(i,NY) + 2*bcVecU(i);
          grid_.V()(i,NY) = bcVecV(i);
      }
//        PROGRESS("Set North as inlet");
  }
  else if(type == 4){
      for (int i = nGhost; i < NX + nGhost; ++i){
          grid_.U()(i,NY) = grid_.U()(i,1);
          grid_.V()(i,NY+1) = grid_.V()(i,2);
      }
//        PROGRESS("Set North as periodic");
  }    
  else{
      for (int i = nGhost; i < NX + nGhost; ++i){
          grid_.V()(i,NY) = 0.;
          grid_.U()(i,NY+1) = -grid_.U()(i,NY) + 2*Uw;
      }
//        PROGRESS("Set North as wall");
  }
}

void EnhancedFluidSimulator::setEBC2D(Array<real>& bcVecU,Array<real>& bcVecV, real Vw, int NX, int NY, int nGhost, int type){
	CHECK_MSG((bcVecU.getSize() == NY + 2*nGhost)&&(bcVecV.getSize() == NY + 2*nGhost),"BC Array not set correctly (EAST)");
  if(type == 1){
      for (int j = nGhost; j < NY + nGhost; ++j){
          grid_.U()(NX,j) = 0;//-grid_.U()(NX-1,j);
          grid_.V()(NX+1,j) = grid_.V()(NX,j);
      }
//        PROGRESS("Set East as symmetry");
  }
  else if(type == 2){
      for (int j = nGhost; j < NY + nGhost; ++j){
          grid_.U()(NX,j) = grid_.U()(NX-1,j);
          grid_.V()(NX+1,j) = grid_.V()(NX,j);
      }
//        PROGRESS("Set East as outflow");
  }
  else if(type == 3){
      for (int j = nGhost; j < NY + nGhost; ++j){
          grid_.U()(NX,j) = bcVecU(j);
          grid_.V()(NX+1,j) = -grid_.V()(NX,j) + 2*bcVecV(j);
      }
//        PROGRESS("Set East as inlet");
  }
  else if(type == 4){
      for (int j = nGhost; j < NY + nGhost; ++j){
          grid_.U()(NX,j) = grid_.U()(1,j);
          grid_.V()(NX+1,j) = grid_.V()(2,j);
      }
//        PROGRESS("Set East as periodic");
  }    
  else{
      for (int j = nGhost; j < NY + nGhost; ++j){
          grid_.V()(NX+1,j) = -grid_.V()(NX,j) + 2*Vw;
          grid_.U()(NX,j) = 0.;
      }
//        PROGRESS("Set East as wall");
  }
}

void EnhancedFluidSimulator::setSBC2D(Array<real>& bcVecU,Array<real>& bcVecV, real Uw, int NX, int NY, int nGhost, int type){
	CHECK_MSG((bcVecU.getSize() == NX + 2*nGhost)&&(bcVecV.getSize() == NX + 2*nGhost),"BC Array not set correctly (SOUTH)");
  if(type == 1){
      for (int i = nGhost; i < NX + nGhost; ++i){
          grid_.U()(i,0) = grid_.U()(i,1);
          grid_.V()(i,0) = 0; //-grid_.V()(i,1);
      }
//        PROGRESS("Set South as symmetry");
  }
  else if(type == 2){
      for (int i = nGhost; i < NX + nGhost; ++i){
          grid_.U()(i,0) = grid_.U()(i,1);
          grid_.V()(i,0) = grid_.V()(i,1);
      }
//        PROGRESS("Set South as outflow");
  }
  else if(type == 3){
      for (int i = nGhost; i < NX + nGhost; ++i){
          grid_.U()(i,0) = -grid_.U()(i,1) + 2*bcVecU(i);
          grid_.V()(i,0) = bcVecV(i);
      }
//        PROGRESS("Set South as inlet");
  }
  else if(type == 4){
      for (int i = nGhost; i < NX + nGhost; ++i){
          grid_.V()(i,0) = grid_.V()(i,NY-1);
          grid_.U()(i,0) = grid_.U()(i,NY-1);
          grid_.U()(i,1) = grid_.U()(i,NY);
      }
//        PROGRESS("Set South as periodic");
  }     
  else{
      for (int i = nGhost; i < NX + nGhost; ++i){
          grid_.V()(i,0) = 0.;
          grid_.U()(i,0) = -grid_.U()(i,1) + 2*Uw;
      }
//        PROGRESS("Set South as wall");
  }
}

void EnhancedFluidSimulator::setWBC2D(Array<real>& bcVecU,Array<real>& bcVecV, real Vw, int NX, int NY, int nGhost, int type){
	CHECK_MSG((bcVecU.getSize() == NY + 2*nGhost)&&(bcVecV.getSize() == NY + 2*nGhost),"BC Array not set correctly (WEST)");
  if(type == 1){
      for (int j = nGhost; j < NY + nGhost; ++j){
          grid_.U()(0,j) = 0; //-grid_.U()(1,j);
          grid_.V()(0,j) = grid_.V()(1,j);
      }
//        PROGRESS("Set West as symmetry");
  }
  else if(type == 2){
      for (int j = nGhost; j < NY + nGhost; ++j){
          grid_.U()(0,j) = grid_.U()(1,j);
          grid_.V()(0,j) = grid_.V()(1,j);
      }
//        PROGRESS("Set West as outflow");
  }
  else if(type == 3){
      for (int j = nGhost; j < NY + nGhost; ++j){
          grid_.U()(0,j) = bcVecU(j);
          grid_.V()(0,j) = -grid_.V()(1,j) + 2*bcVecV(j);
      }
//        PROGRESS("Set West as inlet");
  }
  else if(type == 4){
      for (int j = nGhost; j < NY + nGhost; ++j){
          grid_.U()(0,j) = grid_.U()(NX-1,j);
          grid_.V()(0,j) = grid_.V()(NX-1,j);
          grid_.V()(1,j) = grid_.V()(NX,j);
      }
//        PROGRESS("Set West as periodic");
  }    
  else{
      for (int j = nGhost; j < NY + nGhost; ++j){
          grid_.V()(0,j) = -grid_.V()(1,j) + 2*Vw;
          grid_.U()(0,j) = 0.;
      }
//        PROGRESS("Set West as wall");
  }
}

inline real EnhancedFluidSimulator::dU2_dx(int i, int j)
{
   return (inv4dx * (  pow(grid_.u(i,j,CENTER)+grid_.u(i,j,EAST),2) - pow(grid_.u(i,j,WEST)+grid_.u(i,j,CENTER),2) 
                             + gamma * (  std::abs(grid_.u(i,j,CENTER)+grid_.u(i,j,EAST))
                                          * (grid_.u(i,j,CENTER)-grid_.u(i,j,EAST))
                                        - std::abs(grid_.u(i,j,WEST)+grid_.u(i,j,CENTER))
                                          * (grid_.u(i,j,WEST)-grid_.u(i,j,CENTER)) ) ) );
}
inline real EnhancedFluidSimulator::dV2_dy(int i, int j)
{
   return (inv4dy * (  pow(grid_.v(i,j,CENTER)+grid_.v(i,j,NORTH),2) - pow(grid_.v(i,j,SOUTH)+grid_.v(i,j,CENTER),2) 
                             + gamma * (  std::abs(grid_.v(i,j,CENTER)+grid_.v(i,j,NORTH))
                                          * (grid_.v(i,j,CENTER)-grid_.v(i,j,NORTH)) 
                                        - std::abs(grid_.v(i,j,SOUTH)+grid_.v(i,j,CENTER))
                                          * (grid_.v(i,j,SOUTH)-grid_.v(i,j,CENTER)) ) ) );
}
inline real EnhancedFluidSimulator::dUV_dy(int i, int j)
{
   return (inv4dy * (  (grid_.v(i,j,CENTER)+grid_.v(i,j,EAST))*(grid_.u(i,j,CENTER)+grid_.u(i,j,NORTH))
                      - (grid_.v(i,j,SOUTH)+grid_.v(i+1,j,SOUTH))*(grid_.u(i,j,SOUTH)+grid_.u(i,j,CENTER))
                      + gamma * (  std::abs(grid_.v(i,j,CENTER)+grid_.v(i,j,EAST)) * (grid_.u(i,j,CENTER)-grid_.u(i,j,NORTH))  
                                 - std::abs(grid_.v(i,j,SOUTH)+grid_.v(i+1,j,SOUTH)) * (grid_.u(i,j,SOUTH)-grid_.u(i,j,CENTER)) ) ) );
}
inline real EnhancedFluidSimulator::dVU_dx(int i, int j)
{
   return (inv4dx * (  (grid_.u(i,j,CENTER)+grid_.u(i,j,NORTH))*(grid_.v(i,j,CENTER)+grid_.v(i,j,EAST))
                      - (grid_.u(i,j,WEST)+grid_.u(i,j+1,WEST))*(grid_.v(i,j,WEST)+grid_.v(i,j,CENTER)) 
                      + gamma * (  std::abs(grid_.u(i,j,CENTER)+grid_.u(i,j,NORTH)) * (grid_.v(i,j,CENTER)-grid_.v(i,j,EAST)) 
                                 - std::abs(grid_.u(i,j,WEST)+grid_.u(i,j+1,WEST)) * (grid_.v(i,j,WEST)-grid_.v(i,j,CENTER)) ) ) );
}
inline real EnhancedFluidSimulator::d2U_dx2(int i, int j)
{
   return (dx_2 * ( grid_.u(i,j,EAST)-2.0*grid_.u(i,j,CENTER)+grid_.u(i,j,WEST) ) );
}
inline real EnhancedFluidSimulator::d2U_dy2(int i, int j)
{
   return (dy_2 * ( grid_.u(i,j,NORTH)-2.0*grid_.u(i,j,CENTER)+grid_.u(i,j,SOUTH) ) );
}
inline real EnhancedFluidSimulator::d2V_dx2(int i, int j)
{
   return (dx_2 * ( grid_.v(i,j,EAST)-2.0*grid_.v(i,j,CENTER)+grid_.v(i,j,WEST) ) );
}
inline real EnhancedFluidSimulator::d2V_dy2(int i, int j)
{
   return (dy_2 * ( grid_.v(i,j,NORTH)-2.0*grid_.v(i,j,CENTER)+grid_.v(i,j,SOUTH) ) );
}
