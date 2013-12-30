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
	Ntimes_ = conf.getIntParameter("timesteps");
	gx_ = conf.getRealParameter("gx");
	gy_ = conf.getRealParameter("gy");
	gam_ = conf.getRealParameter("gamma");
	RE_ = conf.getRealParameter("Re");
	uinit_ = conf.getRealParameter("U_init");
	vinit_ = conf.getRealParameter("V_init");
	pinit_ = conf.getRealParameter("P_init");		
}

/// Simulates a given time-length
void FluidSimulator::simulate             ( real duration              ){
	FluidSimulator::initFields();
	bool resp;
	int i = 0;
	real time = 0.;
	FluidSimulator::setUVBC2D();	
	FluidSimulator::computeFG();
	FluidSimulator::composeRHS2D();
	while (time < duration){
		FluidSimulator::setUVBC2D();
		FluidSimulator::computeFG();
		FluidSimulator::composeRHS2D();
		resp = psolver_.solve(grid_);
//		ASSERT(resp);
		FluidSimulator::adapUV2D();
		time += dt_;
		++i;
		std::cout << "Timestep: " << i << "\tDelta_t: " << dt_ << "\tPres. Iterations: " << psolver_.iterations() << std::endl;		
	}
	PROGRESS("Finished time simulations\n");

}
void FluidSimulator::simulateTimeStepCount( unsigned int nrOfTimeSteps ){
	FluidSimulator::initFields();
	bool resp;
	FluidSimulator::setUVBC2D();	
	FluidSimulator::computeFG();
	FluidSimulator::composeRHS2D();	
	for (int i = 0; i < nrOfTimeSteps; ++i){
		FluidSimulator::setUVBC2D();
		FluidSimulator::computeFG();
		FluidSimulator::composeRHS2D();
		resp = psolver_.solve(grid_);
//		ASSERT(resp);
		FluidSimulator::adapUV2D();
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

void FluidSimulator::computeFG(){
	Array p = grid_.p();
	Array U = grid_.U();
	Array V = grid_.V();
	Array F = grid_.F();
	Array G = grid_.G();
	duvdy = du2dx = d2udx2 = d2udy2 = dpdx = duvdx = dv2dy = d2vdx2 = d2vdy2 = dpdy = p;	
	real dx = grid_.dx();
	real dx2 = pow(dx,2);
	real dy = grid_.dy();
	real dy2 = pow(dy,2);
	real gam = gam_;
	real dt = dt_;
	real RE = RE_;
	real gx = gx_;
	real gy = gy_;
	int NGhost = grid_.NGhost();
	int Nx = grid_.Nx();
	int Ny = grid_.Ny();
	for (int i = NGhost; i < Nx + NGhost; ++i)
		for (int j = NGhost; j < Ny + NGhost; ++j){
			du2dx(i,j) = 1/dx*(pow((U(i,j)+U(i+1,j))/2,2) - pow((U(i-1,j)+U(i,j))/2,2)) + gam/dx*(std::abs(U(i,j) + U(i+1,j))*(U(i,j) - U(i+1,j))/4 - std::abs(U(i-1,j) + U(i,j))*(U(i-1,j) - U(i,j))/4);
			duvdy(i,j) = 1/dy*((V(i,j) + V(i+1,j))*(U(i,j) + U(i,j+1))/4 - (V(i,j-1) + V(i+1,j-1))*(U(i,j-1) + U(i,j))/4) + gam/dy*(std::abs(V(i,j) + V(i+1,j))*(U(i,j) - U(i,j+1))/4 - std::abs(V(i,j-1) + V(i+1,j-1))*(U(i,j-1) - U(i,j))/4);
			d2udx2(i,j) = (U(i+1,j) - 2*U(i,j) + U(i-1,j))/dx2;
			d2udy2(i,j) = (U(i,j+1) - 2*U(i,j) + U(i,j-1))/dy2;
			dpdx(i,j) = (p(i+1,j) - p(i,j))/dx;
			duvdx(i,j) = 1/dx*((V(i,j) + V(i+1,j))*(U(i,j) + U(i,j+1))/4 - (V(i-1,j) + V(i,j))*(U(i-1,j) + U(i-1,j+1))/4) + gam/dx*(std::abs(U(i,j) + U(i,j+1))*(V(i,j) - V(i+1,j))/4 - std::abs(U(i-1,j) + U(i-1,j+1))*(V(i-1,j) - V(i,j))/4);
			dv2dy(i,j) = 1/dy*(pow((V(i,j)+V(i,j+1))/2,2) - pow((V(i,j-1)+V(i,j))/2,2)) + gam/dy*(std::abs(V(i,j) + V(i,j+1))*(V(i,j) - V(i,j+1))/4 - std::abs(V(i,j-1) + V(i,j))*(V(i,j-1) - V(i,j))/4);
			d2vdx2(i,j) = (V(i+1,j) - 2*V(i,j) + V(i-1,j))/dx2;
			d2vdy2(i,j) = (V(i,j+1) - 2*V(i,j) + V(i,j-1))/dy2;
			dpdy(i,j) = (p(i,j+1) - p(i,j))/dy;
			grid_.F()(i,j) = U(i,j) + dt*((d2udx2(i,j) + d2udy2(i,j))/RE - du2dx(i,j) - duvdy(i,j) + gx);
			grid_.G()(i,j) = V(i,j) + dt*((d2vdx2(i,j) + d2vdy2(i,j))/RE - dv2dy(i,j) - duvdx(i,j) + gy);
		}
//	PROGRESS("Setting BC for F,G");
/*	Array vecx(Nx);
	Array vecy(Ny);
	vecx.fill(0);
	Array vecx1(vecx);	
	vecy.fill(0);
	Array vecy1(vecy);	
	vecy = U.row(0);
	vecy1 = U.row(Ny);
	grid_.F().col(0,vecy);
	grid_.F().col(Ny,vecy1);
	vecx = V.row(0);
	vecx1 = V.row(Nx);
	grid_.G().row(0,vecx);
	grid_.G().row(Nx,vecx1); */
	for (int j = 1; j < Ny+1; ++j){
		grid_.F()(0,j) = grid_.U()(0,j);
		grid_.F()(Nx,j) = grid_.U()(Nx,j);
	}
	for (int i = 1; i < Nx+1; ++i){
		grid_.G()(i,0) = grid_.V()(i,0);
		grid_.G()(i,Ny) = grid_.V()(i,Ny);
	}
//	PROGRESS("Done setting BC for F,G");	

//	PROGRESS("All F and G allocated");
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
	setDelta_t2D();
	real dt = dt_;
	real dx = grid_.dx();
	real dy = grid_.dy();
	for (int i = grid_.NGhost(); i < grid_.Nx() + grid_.NGhost(); ++i)
		for (int j = grid_.NGhost(); j < grid_.Ny() + grid_.NGhost(); ++j){
			grid_.U()(i,j) = grid_.F()(i,j) - dt/dx*(grid_.p()(i+1,j) - grid_.p()(i,j));
			grid_.V()(i,j) = grid_.G()(i,j) - dt/dy*(grid_.p()(i,j+1) - grid_.p()(i,j));
		}
}

void FluidSimulator::setDelta_t2D(){
	real tau = 0.9;//tau_; Changed to hardcode as 0.9 for now
	dt_ = std::min(RE_/(2*(1/pow(grid_.dx(),2) + 1/pow(grid_.dy(),2))),grid_.dx()/grid_.U().amax());
	dt_ = tau*std::min(dt_,grid_.dy()/grid_.V().amax());
}

void FluidSimulator::setUVBC2D(){
    real U_NBC = 0.1,V_NBC = 0.0, U_EBC = 0.0,V_EBC = 0.0, U_SBC = 0.0,V_SBC = 0.0, U_WBC= 0.0,V_WBC = 0.0;
	int NGhost = grid_.NGhost();
	int nx = grid_.Nx();
	int ny = grid_.Ny();
	Array tempUx(nx+2*NGhost), tempVx(nx+2*NGhost);
	Array tempUy(ny+2*NGhost), tempVy(ny+2*NGhost);	
	tempUx.fill(U_NBC);
	tempVx.fill(V_NBC);
	FluidSimulator::setNBC2D(tempUx,tempVx,nx,ny,NGhost,1);
	tempUy.fill(U_EBC);
	tempVy.fill(V_EBC);
	FluidSimulator::setEBC2D(tempUy,tempVy,nx,ny,NGhost,1);
	tempUx.fill(U_SBC);
	tempVx.fill(V_SBC);
	FluidSimulator::setSBC2D(tempUx,tempVx,nx,ny,NGhost,1);
	tempUy.fill(U_WBC);
	tempVy.fill(V_WBC);
	FluidSimulator::setWBC2D(tempUy,tempVy,nx,ny,NGhost,1);			
}

void FluidSimulator::setNBC2D(Array& bcVecU,Array& bcVecV, int NX, int NY, int nGhost, int type){
	Array tempU;
	tempU = bcVecU + bcVecU - grid_.U().col(NY);
	grid_.U().col(NY+1,tempU); //Only no-slip for now, use type to change
    Array tempV(NY);
    tempV.fill(0);
    grid_.V().col(NY,tempV);
}

void FluidSimulator::setEBC2D(Array& bcVecU,Array& bcVecV, int NX, int NY, int nGhost, int type){
	Array tempV;
	tempV = bcVecV + bcVecV - grid_.V().row(NX);
	grid_.V().row(NX+1,tempV);
    Array tempU(NX);
    tempU.fill(0);
    grid_.U().row(NX,tempU);
}

void FluidSimulator::setSBC2D(Array& bcVecU,Array& bcVecV, int NX, int NY, int nGhost, int type){
	Array tempU;
	tempU = bcVecU + bcVecU - grid_.U().col(1);
    Array tempV(NY);
    tempV.fill(0);
    grid_.V().col(0,tempV);
	grid_.U().col(0,tempU);
}

void FluidSimulator::setWBC2D(Array& bcVecU,Array& bcVecV, int NX, int NY, int nGhost, int type){
	Array tempV;
	tempV = bcVecV + bcVecV - grid_.V().row(1);
    Array tempU(NX);
    tempU.fill(0);
	grid_.V().row(0,tempV);
    grid_.U().row(0,tempU);
}
