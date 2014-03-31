#include <SORSolver.hh>

SORSolver::SORSolver(){
}
// Constructor to manually create SORSolver
SORSolver::SORSolver ( int iterMax, real TOL, real W, std::string nrmtype ){
	iterMax_ = iterMax;
	TOL_ = TOL;
	W_ = W;
	nrmtype_ = nrmtype;
}

// Constructor to create a SORSolver from a parsed configuration file
SORSolver::SORSolver ( const FileReader & configuration ){
	iterMax_ = configuration.getIntParameter("itermax");
	TOL_ = configuration.getRealParameter("eps");
	W_ = configuration.getRealParameter("omg");
	rescheckfreq_ = configuration.getIntParameter("checkfrequency");
	nrmtype_ = "2"; // Hard coded to comply with SORTest.c
  Hbctype_ = configuration.getIntParameter("boundary_condition_E");
  Vbctype_ = configuration.getIntParameter("boundary_condition_S");
//	nrmtype_ = configuration.getStringParameter("NormType");
}

int SORSolver::iterations(){
	return iter_;
}


// solve the pressure equation on the staggered grid (Older version: need to change)
bool SORSolver::solve( StaggeredGrid & grid ){
	real tol=1;
	int Nx = grid.Nx();
	int Ny = grid.Ny();
	int nfluid = Nx*Ny;
	int NGhost = grid.NGhost();
	real unew, uold, ubar, add;//, resnrm;

	real Wm1 = 1.0 - W_;
	real dx = grid.dx();
	real dy = grid.dy();
	real odx2 = 1/pow(dx,2);
	real ody2 = 1/pow(dy,2);
	real oddiag = 2*odx2 + 2*ody2;
	Array<real> &u_sol = grid.p();
	Array<real> &rhs_sol = grid.rhs();
//	Array<real> res(Nx,Ny);
//	Array<real> res_in(res.getSize());	
	
	iter_ = 0;
	while((tol > TOL_)&&(iter_ < iterMax_)){
		++iter_;
//		u_old = u_sol;		
		setBC(u_sol,NGhost);
//		PROGRESS("Starting solution sweep");
		for (int i = 0; i < Nx; ++i){
			for (int j = 0; j < Ny; ++j)
			{
				uold = u_sol(i+NGhost,j+1+NGhost)*ody2 + u_sol(i+1+NGhost,j+NGhost)*odx2;
				unew = u_sol(i+NGhost,j-1+NGhost)*ody2 + u_sol(i-1+NGhost,j+NGhost)*odx2;
				ubar = (uold + unew - rhs_sol(i+NGhost,j+NGhost))/(oddiag);
				u_sol(i+NGhost,j+NGhost) = W_*ubar + Wm1*u_sol(i+NGhost,j+NGhost);
			}
		}
//		PROGRESS("Done sweeping interior points");
		//Update BC
		setBC(u_sol,NGhost);
		if (iter_%rescheckfreq_ == 0){
			tol = 0.;
			for (int i = 0; i < Nx; ++i)
				for (int j = 0; j < Ny; ++j){
					add = -rhs_sol(i+NGhost,j+NGhost) - (oddiag)*(u_sol(i+NGhost,j+NGhost)) + ody2*(u_sol(i+NGhost,j+1+NGhost) + u_sol(i+NGhost,j-1+NGhost)) + odx2*(u_sol(i+1+NGhost,j+NGhost) + u_sol(i-1+NGhost,j+NGhost));
					tol += add*add;
				} 
			tol = sqrt((tol)/nfluid);
//			std::cout << "Iteration: " << iter_ << "\t Residual: " << tol << std::endl;
		}
      if (iter_ %100 == 0 || iter_ == 2)
		std::cout<<"Iteration no = "<< iter_ << "\tResidual = "<< tol <<"\n";
	}
	if(tol < TOL_){
		residual = tol;
		return true;
	}
	else{
		residual = tol;
		return false;
	}
}

bool SORSolver::solve(StaggeredGrid& grid, Geometry2D& mesh){
	int i,j,nfluid;
	real rdx2,rdy2;
	real add,beta_2,tol;
	int imax = grid.Nx();
	int jmax = grid.Ny();
	real delx = grid.dx();
	real dely = grid.dy();
	real omg = W_;
	Array<real> &P = grid.p();
	Array<real> &RHS = grid.rhs();
	rdx2 = 1./delx/delx;
	rdy2 = 1./dely/dely;
	beta_2 = -omg/(2.0*(rdx2+rdy2));
//	maxP = P.maximum();
	nfluid = mesh.geom_nFluids();
	tol = 1.0;

	for (iter_=0;(tol > TOL_)&&(iter_ < iterMax_);++iter_){

		for (i=0;i<=imax+1;++i)
			for (j=0;j<=jmax+1;++j)
				if (mesh.geom_isFluid(i,j))
					P(i,j) = (1.-omg)*P(i,j) - beta_2*( (P(i+1,j)+P(i-1,j)) *rdx2
                                  +(P(i,j+1)+P(i,j-1)) *rdy2
                                  - RHS(i,j) );

		setBC(P,mesh);

		if (iter_%rescheckfreq_ == 0){
			tol = 0.0;
			setBC(P,mesh);	
			for (i=0;i<=imax+1;i++)
				for (j=0;j<=jmax+1;j++)
					if (mesh.geom_isFluid(i,j))   
					{
						add =  (P(i+1,j)-2*P(i,j)+P(i-1,j))*rdx2 + (P(i,j+1)-2*P(i,j)+P(i,j-1))*rdy2 - RHS(i,j);
						tol += add*add;
//						tol = (tol < std::abs(add))? std::abs(add): tol;
					}
				tol = sqrt((tol)/nfluid);
		}
      if (iter_ %100 == 0 || iter_ == 2)
		std::cout<<"Iteration no = "<< iter_ << "\tResidual = "<< tol <<"\n";
	}
	if(tol < TOL_){
		residual = tol;
		return true;
	}
	else{
		residual = tol;
		return false;
	}
}

real SORSolver::residualnorm(){ return residual; }

void SORSolver::setBC(Array<real>& Q,Geometry2D& mesh){
	setBC(Q,1);
	int nx = Q.getSize(0);
	int ny = Q.getSize(1);	
	for (int i=0;i < nx;++i)
		for (int j=0;j < ny;++j)
		  if (mesh(i,j) >= B_N && mesh(i,j) <= B_SE) 
			 switch (mesh(i,j))
			{
				case B_N:{  Q(i,j) = Q(i,j+1);                 break;}
				case B_E:{  Q(i,j) = Q(i+1,j);                 break;}
				case B_S:{  Q(i,j) = Q(i,j-1);                 break;} 
				case B_W:{  Q(i,j) = Q(i-1,j);                 break;}
				case B_NE:{ Q(i,j) = 0.5*(Q(i,j+1)+Q(i+1,j)); break;}
				case B_SE:{ Q(i,j) = 0.5*(Q(i,j-1)+Q(i+1,j)); break;}
				case B_SW:{ Q(i,j) = 0.5*(Q(i,j-1)+Q(i-1,j)); break;}
				case B_NW:{ Q(i,j) = 0.5*(Q(i,j+1)+Q(i-1,j)); break;}
				default:                                         break;
			}
}

void SORSolver::setBC(Array<real> &Q,const int & nGhost){
	int nx = Q.getSize(0);
	int ny = Q.getSize(1);
	for (int i = 1; i < nx-1; ++i){
		if(Vbctype_ == 4){ //Periodic BC treatment for pressure equation (N-S)
			Q(i,1) = Q(i,ny-2);
			Q(i,0) = Q(i,1);
			Q(i,ny-1) = Q(i,ny-2);
		}
		else{
			Q(i,0) = Q(i,1);
			Q(i,ny-1) = Q(i,ny-2);
		}
	}
	for (int j = 1; j < ny-1; ++j){
		if(Hbctype_ == 4){ //Periodic BC treatment for pressure equation (E-W)
			Q(0,j) = Q(1,j);
			Q(nx-1,j) = Q(nx-2,j);		
			Q(1,j) = Q(nx-2,j);
		}
		else{
			Q(0,j) = Q(1,j);
			Q(nx-1,j) = Q(nx-2,j);
		}
	}
}
