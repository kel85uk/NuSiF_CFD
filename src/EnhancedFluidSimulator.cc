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

#include "EnhancedFluidSimulator.hh"
#include "Array.hh"
#include "GrayScaleImage.hh"

EnhancedFluidSimulator::EnhancedFluidSimulator(){
}

EnhancedFluidSimulator::EnhancedFluidSimulator( const FileReader & conf ){
	StaggeredGrid domain(conf);
	SORSolver pdomain(conf);
	Geometry2D mesh(domain);
	grid_ = domain;
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
	meshvis_ = mesh_;
	mesh_.geom_finalize();	
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

  while ((time < duration)&&(i<Ntimes_)){
		EnhancedFluidSimulator::setDelta_t2D();  
		EnhancedFluidSimulator::setUVBC2D();
		EnhancedFluidSimulator::computeFG();
		EnhancedFluidSimulator::composeRHS2D();
		psolver_.solve_SOR1(grid_,mesh_);
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
		std::cout << "Timestep: " << i << "\tTime: " << time << "\tDelta_t: " << dt_ << "\tPres. Iterations: " << psolver_.iterations() << std::endl;			
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
		psolver_.solve_CG(grid_,mesh_);
		EnhancedFluidSimulator::adapUV2D();
		if(i%pnorm_ == 0){
			EnhancedFluidSimulator::normalizePressure();
		}
		if(i%writeVTKfreq_ == 0){
			vtkWriter.write(mesh_,time);
			std::cout << "Writing to file\n";
		}
		time += dt_;
		std::cout << "Timestep: " << i << "\tTime: " << time << "\tDelta_t: " << dt_ << "\tPres. Iterations: " << psolver_.iterations() << std::endl;			
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
/** Getter function for the internally stored Geometry (for visualization) with read-write */
Geometry2D & EnhancedFluidSimulator::MESHVIS(){
	return meshvis_;
}
/** Getter function for the internally stored Geometry (for visualization) with read-only */
const Geometry2D & EnhancedFluidSimulator::MESHVIS() const{
	return meshvis_;
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
	real p_mean = grid_.p().sumI(0)/((grid_.Nx()+1)*(grid_.Ny()+1));
	for (int i = 0; i < grid_.p().getSize(); ++i){
		grid_.p()(i) = grid_.p()(i) - p_mean;
	}
	std::cout << "Normalizing pressure\n";
}
/** Calculates the initial guess of the velocities without pressure correction */
void EnhancedFluidSimulator::computeFG(){
	Array &p = grid_.p();
	Array &U = grid_.U();
	Array &V = grid_.V();
	Array &F = grid_.F();
	Array &G = grid_.G();
	duvdy = du2dx = d2udx2 = d2udy2 = dpdx = duvdx = dv2dy = d2vdx2 = d2vdy2 = dpdy = p;	
	real dx = grid_.dx();
  real dx2 = dx*dx;
	real dy = grid_.dy();
  real dy2 = dy*dy;
	real gamma = gam_;
//	real dt = dt_;
	real RE = RE_;
	real gx = gx_;
	real gy = gy_;
	int Nx = grid_.Nx(),Ny = grid_.Ny(),i,j;
	for (i=1;i<=Nx-1;++i)
	  for (j=1;j<=Ny;++j)
      {
       /* only if both adjacent cells are fluid cells */
       if( ((mesh_(i,j) & C_F)) &&
           ((mesh_(i+1,j) & C_F)) )
         {
          du2dx(i,j) = ((U(i,j)+U(i+1,j))*(U(i,j)+U(i+1,j)) + gamma*std::abs(U(i,j)+U(i+1,j))*(U(i,j)-U(i+1,j)) - (U(i-1,j)+U(i,j))*(U(i-1,j)+U(i,j)) - gamma*std::abs(U(i-1,j)+U(i,j))*(U(i-1,j)-U(i,j)))/(4.0*dx);
          duvdy(i,j) = ((V(i,j)+V(i+1,j))*(U(i,j)+U(i,j+1)) + gamma*std::abs(V(i,j)+V(i+1,j))*(U(i,j)-U(i,j+1)) - (V(i,j-1)+V(i+1,j-1))*(U(i,j-1)+U(i,j)) - gamma*std::abs(V(i,j-1)+V(i+1,j-1))*(U(i,j-1)-U(i,j)))/(4.0*dy);
          d2udx2(i,j) = (U(i+1,j)-2.0*U(i,j)+U(i-1,j))/dx2;
          d2udy2(i,j) = (U(i,j+1)-2.0*U(i,j)+U(i,j-1))/dy2;
   
          F(i,j) = U(i,j)+dt_*((d2udx2(i,j)+d2udy2(i,j))/RE-du2dx(i,j)-duvdy(i,j)+gx);
         }
       else
          F(i,j) = U(i,j);
      }

 for (i=1;i<=Nx;++i)
    for (j=1;j<=Ny-1;++j)
      {
       /* only if both adjacent cells are fluid cells */
       if( ((mesh_(i,j) & C_F)) &&
           ((mesh_(i,j+1) & C_F)) )
         {
          duvdx(i,j) = ((U(i,j)+U(i,j+1))*(V(i,j)+V(i+1,j)) + gamma*std::abs(U(i,j)+U(i,j+1))*(V(i,j)-V(i+1,j)) - (U(i-1,j)+U(i-1,j+1))*(V(i-1,j)+V(i,j)) - gamma*std::abs(U(i-1,j)+U(i-1,j+1))*(V(i-1,j)-V(i,j)))/(4.0*dx);
          dv2dy(i,j) = ((V(i,j)+V(i,j+1))*(V(i,j)+V(i,j+1)) + gamma*std::abs(V(i,j)+V(i,j+1))*(V(i,j)-V(i,j+1)) - (V(i,j-1)+V(i,j))*(V(i,j-1)+V(i,j)) - gamma*std::abs(V(i,j-1)+V(i,j))*(V(i,j-1)-V(i,j)))/(4.0*dy);

          d2vdx2(i,j) = (V(i+1,j)-2.0*V(i,j)+V(i-1,j))/dx2;
          d2vdy2(i,j) = (V(i,j+1)-2.0*V(i,j)+V(i,j-1))/dy2;

          G(i,j) = V(i,j)+dt_*((d2vdx2(i,j)+d2vdy2(i,j))/RE-duvdx(i,j)-dv2dy(i,j)+gy);	      
         }
       else
          G(i,j) = V(i,j);
      }
 /* F and G at external boundary */
 /*------------------------------*/ 
 for (j=1;j<=Ny;++j)
   {
    F(0,j)    = U(0,j);
    F(Nx,j) = U(Nx,j);
   }
 for (i=1;i<=Nx;++i)
   {
    G(i,0)    = V(i,0);
    G(i,Ny) = V(i,Ny);
   }
}

void EnhancedFluidSimulator::composeRHS2D(){
	real dt = dt_;
	real dx = grid_.dx();
	real dy = grid_.dy();
//	PROGRESS("Allocating RHS");
	for (int i = 1; i < grid_.Nx() + 1; ++i)
		for (int j = 1; j < grid_.Ny() + 1; ++j){
			if ((mesh_(i,j) & C_F))
				grid_.rhs()(i,j) = ((grid_.F()(i,j) - grid_.F()(i-1,j))/dx + (grid_.G()(i,j) - grid_.G()(i,j-1))/dy)/dt;
		}
}

void EnhancedFluidSimulator::adapUV2D(){
//	setDelta_t2D();
	real dt = dt_;
	real dx = grid_.dx();
	real dy = grid_.dy();
	for (int i = 1; i <= grid_.Nx()-1; ++i)
		for (int j = 1; j <= grid_.Ny(); ++j)
			if( ((mesh_(i,j) & C_F)) && ((mesh_(i+1,j) & C_F)) )
				grid_.U()(i,j) = grid_.F()(i,j) - dt/dx*(grid_.p()(i+1,j) - grid_.p()(i,j));
	for (int i = 1; i <= grid_.Nx(); ++i)
		for (int j = 1; j<= grid_.Ny()-1; ++j)
			if( ((mesh_(i,j) & C_F)) && ((mesh_(i,j+1) & C_F)) )		
				grid_.V()(i,j) = grid_.G()(i,j) - dt/dy*(grid_.p()(i,j+1) - grid_.p()(i,j));
}

void EnhancedFluidSimulator::setDelta_t2D(){
	real tau = tau_;
	if (tau > 0){
	dt_ = std::min(RE_/(2*(1/pow(grid_.dx(),2) + 1/pow(grid_.dy(),2))),grid_.dx()/grid_.U().amax());
	dt_ = tau*std::min(dt_,grid_.dy()/grid_.V().amax());
	}
}

void EnhancedFluidSimulator::setUVBC2D(){
  real U_NBC = U_NBC_,V_NBC = V_NBC_, U_EBC = U_EBC_,V_EBC = V_EBC_, U_SBC = U_SBC_,V_SBC = V_SBC_, U_WBC= U_WBC_,V_WBC = V_WBC_;
	int NGhost = grid_.NGhost();
	int nx = grid_.Nx();
	int ny = grid_.Ny();

	if (customBC){
		Array U_NBCv(nx + 2*NGhost);Array V_NBCv(nx + 2*NGhost);
		U_NBCv.fill(U_NBC_); V_NBCv.fill(V_NBC_);
		Array U_EBCv(ny + 2*NGhost);Array V_EBCv(ny + 2*NGhost);
		U_EBCv.fill(U_EBC_); V_EBCv.fill(V_EBC_);
		Array U_SBCv(nx + 2*NGhost);Array V_SBCv(nx + 2*NGhost);
		U_SBCv.fill(U_SBC_); V_SBCv.fill(V_SBC_);
		Array U_WBCv(ny + 2*NGhost);Array V_WBCv(ny + 2*NGhost);
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
  /* setting the boundary values at inner obstacle cells */
  /*                  (only no-slip)                     */
  /*-----------------------------------------------------*/
	for(int i=1;i<=nx;++i)
		for(int j=1;j<=ny;++j)
			if(mesh_(i,j) & 0x000f)		/* The mask 0x000f filters the */
		   													/* obstacle cells adjacent to  */
		   													/* fluid cells                 */
				switch (mesh_(i,j))
				{
					case B_N:  { 
						grid_.V()(i,j)   = 0.0;
						grid_.U()(i,j)   = -grid_.U()(i,j+1);
						grid_.U()(i-1,j) = -grid_.U()(i-1,j+1);
						break;
					}
					case B_O:  { 
						grid_.U()(i,j)   = 0.0;
						grid_.V()(i,j)   = -grid_.V()(i+1,j);
						grid_.V()(i,j-1) = -grid_.V()(i+1,j-1);
						break;
					}
					case B_S:  { 
						grid_.V()(i,j-1) = 0.0;
						grid_.U()(i,j)   = -grid_.U()(i,j-1);
						grid_.U()(i-1,j) = -grid_.U()(i-1,j-1);
						break;
					}
					case B_W:  { 
						grid_.U()(i-1,j) = 0.0;
						grid_.V()(i,j)   = -grid_.V()(i-1,j);
						grid_.V()(i,j-1) = -grid_.V()(i-1,j-1);
						break;
					}
					case B_NO: { 
						grid_.V()(i,j)   = 0.0;
						grid_.U()(i,j)   = 0.0;
						grid_.V()(i,j-1) = -grid_.V()(i+1,j-1);
						grid_.U()(i-1,j) = -grid_.U()(i-1,j+1);
						break;
					}
					case B_SO: { 
						grid_.V()(i,j-1) = 0.0;
						grid_.U()(i,j)   = 0.0;
						grid_.V()(i,j)   = -grid_.V()(i+1,j);
						grid_.U()(i-1,j) = -grid_.U()(i-1,j-1);
						break;
					}
					case B_SW: { 
						grid_.V()(i,j-1) = 0.0;
						grid_.U()(i-1,j) = 0.0;
						grid_.V()(i,j)   = -grid_.V()(i-1,j);
						grid_.U()(i,j)   = -grid_.U()(i,j-1);
						break;
					}
					case B_NW: { 
						grid_.V()(i,j)   = 0.0;
						grid_.U()(i-1,j) = 0.0;
						grid_.V()(i,j-1) = -grid_.V()(i-1,j-1);
						grid_.U()(i,j)   = -grid_.U()(i,j+1);
						break;
					}
					default : break;
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

void EnhancedFluidSimulator::setNBC2D(Array& bcVecU,Array& bcVecV, real Uw, int NX, int NY, int nGhost, int type){
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

void EnhancedFluidSimulator::setEBC2D(Array& bcVecU,Array& bcVecV, real Vw, int NX, int NY, int nGhost, int type){
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

void EnhancedFluidSimulator::setSBC2D(Array& bcVecU,Array& bcVecV, real Uw, int NX, int NY, int nGhost, int type){
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

void EnhancedFluidSimulator::setWBC2D(Array& bcVecU,Array& bcVecV, real Vw, int NX, int NY, int nGhost, int type){
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
