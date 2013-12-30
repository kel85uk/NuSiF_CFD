#include "FluidSimulator.hh"
#include "Array.hh"

FluidSimulator::FluidSimulator(){
}

FluidSimulator::FluidSimulator( const FileReader & conf ){
	StaggeredGrid domain(conf);
	SORSolver pdomain(conf);
	grid_ = domain;
	psolver_ = pdomain;
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

/// Simulates a given time-length
void FluidSimulator::simulate             ( real duration              ){
	FluidSimulator::initFields();
//	bool resp;
	int i = 0;
	real time = 0.;
	duration = t_end;
	VTKWriter vtkWriter ( grid_, "lidDrivenCavity", true, true);	
  while ((time < duration)&&(i<Ntimes_)){
		FluidSimulator::setDelta_t2D();  
		FluidSimulator::setUVBC2D();
		FluidSimulator::computeFG();
		FluidSimulator::composeRHS2D();
		psolver_.solve(grid_);
//        ASSERT(resp);
		FluidSimulator::adapUV2D();
		if(i%pnorm_ == 0){
			FluidSimulator::normalizePressure();
		}
		if(i%writeVTKfreq_ == 0){
			vtkWriter.write();
			std::cout << "Writing to file\n";
		}		
		time += dt_;
		++i;
		vtkWriter.write();
		std::cout << "Timestep: " << i << "\tDelta_t: " << dt_ << "\tPres. Iterations: " << psolver_.iterations() << std::endl;		
	}
	PROGRESS("Finished time simulations\n");

}
void FluidSimulator::simulateTimeStepCount( unsigned int nrOfTimeSteps ){
	FluidSimulator::initFields();
//	bool resp;
    nrOfTimeSteps = Ntimes_;
	VTKWriter vtkWriter ( grid_, "lidDrivenCavity", true, true );    
    for (unsigned int i = 0; i < nrOfTimeSteps; ++i){
		FluidSimulator::setUVBC2D();
		FluidSimulator::computeFG();
		FluidSimulator::composeRHS2D();
		psolver_.solve(grid_);
//        ASSERT(resp);
		FluidSimulator::adapUV2D();
		if(i%pnorm_ == 0){
			FluidSimulator::normalizePressure();
		}
		if(i%writeVTKfreq_ == 0){
			vtkWriter.write();
			std::cout << "Writing to file\n";
		}
		std::cout << "Timestep: " << i << "\tDelta_t: " << dt_ << "\tPres. Iterations: " << psolver_.iterations() << std::endl;			
	}
	PROGRESS("Finished time simulations\n");
}


// Getter functions for the internally stored StaggeredGrid
StaggeredGrid & FluidSimulator::grid(){
	return grid_;
}
const StaggeredGrid & FluidSimulator::grid() const{
	return grid_;
}

void FluidSimulator::initFields(){
	real pinit = pinit_, uinit = uinit_, vinit = vinit_;
	grid_.p().fill(pinit);
	grid_.U().fill(uinit);
	grid_.V().fill(vinit);
}

void FluidSimulator::normalizePressure(){
	real p_mean = grid_.p().sumI(0)/((grid_.Nx()+1)*(grid_.Ny()+1));
	for (int i = 0; i < grid_.p().getSize(); ++i){
		grid_.p()(i) = grid_.p()(i) - p_mean;
	}
	std::cout << "Normalizing pressure\n";
}

void FluidSimulator::computeFG(){
	Array &p = grid_.p();
	Array &U = grid_.U();
	Array &V = grid_.V();
	duvdy = du2dx = d2udx2 = d2udy2 = dpdx = duvdx = dv2dy = d2vdx2 = d2vdy2 = dpdy = p;	

	real dx = grid_.dx();
  real dx2 = dx*dx;
	real dy = grid_.dy();
  real dy2 = dy*dy;
	real gam = gam_;
	real dt = dt_;
	real RE = RE_;
	real gx = gx_;
	real gy = gy_;
	int Nx = grid_.Nx();
	int Ny = grid_.Ny();
	for (int i = 1; i <= Nx; ++i)
		for (int j = 1; j <= Ny; ++j){
			if(i != Nx){
				du2dx(i,j) = 1/dx*(pow((U(i,j)+U(i+1,j))/2,2) - pow((U(i-1,j)+U(i,j))/2,2)) + gam/dx*(std::abs(U(i,j) + U(i+1,j))*(U(i,j) - U(i+1,j))/4 - std::abs(U(i-1,j) + U(i,j))*(U(i-1,j) - U(i,j))/4);
				duvdy(i,j) = 1/dy*((V(i,j) + V(i+1,j))*(U(i,j) + U(i,j+1))/4 - (V(i,j-1) + V(i+1,j-1))*(U(i,j-1) + U(i,j))/4) + gam/dy*(std::abs(V(i,j) + V(i+1,j))*(U(i,j) - U(i,j+1))/4 - std::abs(V(i,j-1) + V(i+1,j-1))*(U(i,j-1) - U(i,j))/4);
				d2udx2(i,j) = (U(i+1,j) - 2*U(i,j) + U(i-1,j))/dx2;
				d2udy2(i,j) = (U(i,j+1) - 2*U(i,j) + U(i,j-1))/dy2;
				dpdx(i,j) = (p(i+1,j) - p(i,j))/dx;
				grid_.F()(i,j) = U(i,j) + dt*((d2udx2(i,j) + d2udy2(i,j))/RE - du2dx(i,j) - duvdy(i,j) + gx);
			}
			if(j != Ny){
				duvdx(i,j) = 1/dx*((U(i,j) + U(i,j+1))*(V(i,j) + V(i+1,j))/4 - (U(i-1,j) + U(i-1,j+1))*(V(i-1,j) + V(i,j))/4) + gam/dx*(std::abs(U(i,j) + U(i,j+1))*(V(i,j) - V(i+1,j))/4 - std::abs(U(i-1,j) + U(i-1,j+1))*(V(i-1,j) - V(i,j))/4);
				dv2dy(i,j) = 1/dy*(pow((V(i,j)+V(i,j+1))/2,2) - pow((V(i,j-1)+V(i,j))/2,2)) + gam/dy*(std::abs(V(i,j) + V(i,j+1))*(V(i,j) - V(i,j+1))/4 - std::abs(V(i,j-1) + V(i,j))*(V(i,j-1) - V(i,j))/4);
				d2vdx2(i,j) = (V(i+1,j) - 2*V(i,j) + V(i-1,j))/dx2;
				d2vdy2(i,j) = (V(i,j+1) - 2*V(i,j) + V(i,j-1))/dy2;
				dpdy(i,j) = (p(i,j+1) - p(i,j))/dy;
				grid_.G()(i,j) = V(i,j) + dt*((d2vdx2(i,j) + d2vdy2(i,j))/RE - dv2dy(i,j) - duvdx(i,j) + gy);
			}
		}

	for (int j = 1; j <= Ny; ++j){
		grid_.F()(0,j) = grid_.U()(0,j);
		grid_.F()(Nx,j) = grid_.U()(Nx,j);
	}
	for (int i = 1; i <= Nx; ++i){
		grid_.G()(i,0) = grid_.V()(i,0);
		grid_.G()(i,Ny) = grid_.V()(i,Ny);
	}

}

void FluidSimulator::composeRHS2D(){
	real dt = dt_;
	real dx = grid_.dx();
	real dy = grid_.dy();
//	PROGRESS("Allocating RHS");
	for (int i = grid_.NGhost(); i < grid_.Nx() + grid_.NGhost(); ++i)
		for (int j = grid_.NGhost(); j < grid_.Ny() + grid_.NGhost(); ++j){
			grid_.rhs()(i,j) = 1/dt*((grid_.F()(i,j) - grid_.F()(i-1,j))/dx + (grid_.G()(i,j) - grid_.G()(i,j-1))/dy);
		}
}

void FluidSimulator::adapUV2D(){
//	setDelta_t2D();
	real dt = dt_;
	real dx = grid_.dx();
	real dy = grid_.dy();
	for (int i = 1; i <= grid_.Nx(); ++i)
		for (int j = 1; j <= grid_.Ny(); ++j){
			grid_.U()(i,j) = grid_.F()(i,j) - dt/dx*(grid_.p()(i+1,j) - grid_.p()(i,j));
			grid_.V()(i,j) = grid_.G()(i,j) - dt/dy*(grid_.p()(i,j+1) - grid_.p()(i,j));
		}
}

void FluidSimulator::setDelta_t2D(){
	real tau = tau_;
	if (tau > 0){
	dt_ = std::min(RE_/(2*(1/pow(grid_.dx(),2) + 1/pow(grid_.dy(),2))),grid_.dx()/grid_.U().amax());
	dt_ = tau*std::min(dt_,grid_.dy()/grid_.V().amax());
	}
}

void FluidSimulator::setUVBC2D(){
  real U_NBC = U_NBC_,V_NBC = V_NBC_, U_EBC = U_EBC_,V_EBC = V_EBC_, U_SBC = U_SBC_,V_SBC = V_SBC_, U_WBC= U_WBC_,V_WBC = V_WBC_;
	int NGhost = grid_.NGhost();
	int nx = grid_.Nx();
	int ny = grid_.Ny();
  FluidSimulator::setNBC2D(U_NBC,V_NBC,nx,ny,NGhost,Nbctype_); //Types 1 = Symmetry, 2 = outflow, 3 = inlet, 4 = periodic, Else = Wall with Uw
  FluidSimulator::setEBC2D(U_EBC,V_EBC,nx,ny,NGhost,Ebctype_);
  FluidSimulator::setSBC2D(U_SBC,V_SBC,nx,ny,NGhost,Sbctype_);
  FluidSimulator::setWBC2D(U_WBC,V_WBC,nx,ny,NGhost,Wbctype_);
}

void FluidSimulator::setNBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type){
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

void FluidSimulator::setEBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type){
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

void FluidSimulator::setSBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type){

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

void FluidSimulator::setWBC2D(real Uw,real Vw, int NX, int NY, int nGhost, int type){
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
