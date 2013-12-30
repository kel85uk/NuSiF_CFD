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
	int NGhost = grid.NGhost();
  real unew, uold, ubar;//, resnrm;

	real Wm1 = 1.0 - W_;
	real dx = grid.dx();
	real dy = grid.dy();
	real odx2 = 1/pow(dx,2);
	real ody2 = 1/pow(dy,2);
	real oddiag = 2*odx2 + 2*ody2;
	Array &u_sol(grid.p());
	Array &rhs_sol(grid.rhs());
	Array res(Nx,Ny);
	Array res_in(res.getSize());	
	
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
			for (int i = 0; i < Nx; ++i)
				for (int j = 0; j < Ny; ++j){
					res(i,j) = -rhs_sol(i+NGhost,j+NGhost) - (oddiag)*(u_sol(i+NGhost,j+NGhost)) + ody2*(u_sol(i+NGhost,j+1+NGhost) + u_sol(i+NGhost,j-1+NGhost)) + odx2*(u_sol(i+1+NGhost,j+NGhost) + u_sol(i-1+NGhost,j+NGhost));
					res_in(i + j*Nx) = res(i,j);
				} 
			tol = res_in.norm(nrmtype_)/sqrt(res_in.getSize());
//			std::cout << "Iteration: " << iter_ << "\t Residual: " << tol << std::endl;
		}
	}
	if(tol < TOL_){
		grid.p() = u_sol;
		residual = tol;
		return true;
	}
	else{
		return false;
	}
}

bool SORSolver::solve_CG(StaggeredGrid& grid, Geometry2D& mesh){
	int i,j,nfluid;
	int epsE, epsN, epsW, epsS;
	real rdx2,rdy2;
	real add,beta_2,tol, tol0;
	int imax = grid.Nx();
	int jmax = grid.Ny();
	real delx = grid.dx();
	real dely = grid.dy();
	Array &U = grid.p();
	Array &RHS = grid.rhs();
	rdx2 = 1./delx/delx;
	rdy2 = 1./dely/dely;
	nfluid = mesh.geom_nFluids();
	tol = 1.0;
  real u0 = 0.;
	for (i=1;i<=imax;i++)
		for (j=1;j<=jmax;j++)
 			if ((mesh(i,j) & C_F))
				u0 += U(i,j)*U(i,j);
	u0 = sqrt(u0/nfluid);
	if (u0 < 0.0001)
		u0 = 1.0;	
	/** Set variables for CG */
	Array Ap(U);
	Array Res(U);
	Res.fill(0);
	Ap.fill(0);
	Array p(U);
	p.fill(0);
	real resdot, resdotold, resdot0, pkdot, betak, alphak = 0.;
	
	/** Initialize CG (Calculate residual of initial condition) */
	setBC(U,mesh,1);
	for (i=1;i<=imax;i++)
    for (j=1;j<=jmax;j++)
			if ((mesh(i,j) & C_F))   
		    {
					Res(i,j) = (U(i+1,j)-2*U(i,j)+U(i-1,j))*rdx2 + (U(i,j+1)-2*U(i,j)+U(i,j-1))*rdy2 - RHS(i,j); // res = b - Ax
					p(i,j) = Res(i,j); // p = res
					resdot += Res(i,j)*Res(i,j); // resdot_0
			  }
	resdot0 = resdot;
//	setBC(p,mesh,1);
                    /* CG-iteration */
                    /*---------------*/
	for (iter_=0;(resdot > TOL_*TOL_*resdot0 )&&(iter_ < iterMax_);++iter_){
//		std::cout << iter_ << ":" << tol << std::endl;
		/** Start CG iterations */
		pkdot = 0.;
		for (i = 1; i <= imax; ++i) /* p_k Update search direction vector */
			for (j = 1; j <= jmax; ++j)
				if ((mesh(i,j) & C_F)){
					Ap(i,j) = -(p(i+1,j)-2*p(i,j)+p(i-1,j))*rdx2 - (p(i,j+1)-2*p(i,j)+p(i,j-1))*rdy2;
					pkdot += p(i,j)*(Ap(i,j));
				}
//		setBC(p,mesh,1);
		alphak = resdot/pkdot;
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
				if ((mesh(i,j) & C_F)){
					U(i,j) = U(i,j) + alphak*p(i,j); // Update iterate
				}
		setBC(U,mesh,1);
		if (iter_%50 == 0){
			for (i = 1; i <= imax; ++i)
				for (j = 1; j <= jmax; ++j)
					if ((mesh(i,j) & C_F)){
						Res(i,j) = (U(i+1,j)-2*U(i,j)+U(i-1,j))*rdx2 + (U(i,j+1)-2*U(i,j)+U(i,j-1))*rdy2 - RHS(i,j);
					}
		} 
		else{
			for (i = 1; i <= imax; ++i)
				for (j = 1; j <= jmax; ++j)
					if ((mesh(i,j) & C_F))
						Res(i,j) = Res(i,j) - alphak*Ap(i,j); // Update residual		
		}
		resdotold = resdot; // resdot_k-2
		resdot = 0.;
		for (i = 1; i <= imax; ++i) /* resdot_k-1 */
			for (j = 1; j <= jmax; ++j)
				if ((mesh(i,j) & C_F))
					resdot += Res(i,j)*Res(i,j);
					
		betak = resdot/resdotold; // betak = resdot_k/resdot_k-1
		
		for (i = 1; i <= imax; ++i) /* p_k Update search direction vector */
			for (j = 1; j <= jmax; ++j)
				if ((mesh(i,j) & C_F))
					p(i,j) = Res(i,j) + betak*p(i,j);
//		setBC(p,mesh,1);
		tol = sqrt(resdot/nfluid);
	}
	/** Boundary treatment */
	if((tol < TOL_)||(tol < tol0)){
		residual = tol;
		return true;
	}
	else{
		return false;
	}
}

bool SORSolver::solve_PCG(StaggeredGrid& grid, Geometry2D& mesh){
	int i,j,nfluid;
	int epsE, epsN, epsW, epsS;
	real rdx2,rdy2;
	real add,beta_2,tol, tol0;
	int imax = grid.Nx();
	int jmax = grid.Ny();
	real delx = grid.dx();
	real dely = grid.dy();
	Array &U = grid.p();
	Array &RHS = grid.rhs();
	rdx2 = 1./delx/delx;
	rdy2 = 1./dely/dely;
	nfluid = mesh.geom_nFluids();
	tol = 1.0;
  real u0 = 0.;
	for (i=1;i<=imax;i++)
		for (j=1;j<=jmax;j++)
 			if ((mesh(i,j) & C_F))
				u0 += U(i,j)*U(i,j);
	u0 = sqrt(u0/nfluid);
	if (u0 < 0.0001)
		u0 = 1.0;	
	/** Set variables for CG */
	Array Ap(U);
	Array Res(U);
	Res.fill(0);
	Ap.fill(0);
	Array p(U);
	p.fill(0);
	real resdot, resdotold, resdot0, pkdot, betak, alphak = 0.;
	real diagMi = 1/(2*rdx2 + 2*rdy2);
	real Miplus1 = 0;
	real MiplusNx = 0;
	
	/** Initialize CG (Calculate residual of initial condition) */
	setBC(U,mesh,1);
	for (i=1;i<=imax;i++)
    for (j=1;j<=jmax;j++)
			if ((mesh(i,j) & C_F))   
		    {
					Res(i,j) = (U(i+1,j)-2*U(i,j)+U(i-1,j))*rdx2 + (U(i,j+1)-2*U(i,j)+U(i,j-1))*rdy2 - RHS(i,j); // res = b - Ax
			  }
  for (i = 1; i <= imax; ++i)
  		for (j = 1; j <= jmax; ++j)
  			if((mesh(i,j) & C_F)){
				p(i,j) = Miplus1*(Res(i+1,j) + Res(i-1,j)) + MiplusNx*(Res(i,j+1) + Res(i,j-1)) + diagMi*Res(i,j); // p = res
				resdot += Res(i,j)*p(i,j); // resdot_0
			}
	resdot0 = resdot;
//	setBC(p,mesh,1);
                    /* CG-iteration */
                    /*---------------*/
	for (iter_=0;(resdot > TOL_*TOL_*resdot0 )&&(iter_ < iterMax_);++iter_){
//		std::cout << iter_ << ":" << tol << std::endl;
		/** Start CG iterations */
		pkdot = 0.;
		for (i = 1; i <= imax; ++i) /* p_k Update search direction vector */
			for (j = 1; j <= jmax; ++j)
				if ((mesh(i,j) & C_F)){
					Ap(i,j) = -(p(i+1,j)-2*p(i,j)+p(i-1,j))*rdx2 - (p(i,j+1)-2*p(i,j)+p(i,j-1))*rdy2;
					pkdot += p(i,j)*(Ap(i,j));
				}
//		setBC(p,mesh,1);
		alphak = resdot/pkdot;
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
				if ((mesh(i,j) & C_F)){
					U(i,j) = U(i,j) + alphak*p(i,j); // Update iterate
				}
		setBC(U,mesh,1);
		if (iter_%50 == 0){
			for (i = 1; i <= imax; ++i)
				for (j = 1; j <= jmax; ++j)
					if ((mesh(i,j) & C_F)){
						Res(i,j) = (U(i+1,j)-2*U(i,j)+U(i-1,j))*rdx2 + (U(i,j+1)-2*U(i,j)+U(i,j-1))*rdy2 - RHS(i,j);
					}
		} 
		else{
			for (i = 1; i <= imax; ++i)
				for (j = 1; j <= jmax; ++j)
					if ((mesh(i,j) & C_F))
						Res(i,j) = Res(i,j) - alphak*Ap(i,j); // Update residual		
		}
		resdotold = resdot; // resdot_k-2
		resdot = 0.;
		for (i = 1; i <= imax; ++i) /* resdot_k-1 */
			for (j = 1; j <= jmax; ++j)
				if ((mesh(i,j) & C_F))
					resdot += Res(i,j)*(Miplus1*(Res(i+1,j) + Res(i-1,j)) + MiplusNx*(Res(i,j+1) + Res(i,j-1)) + diagMi*Res(i,j));
					
		betak = resdot/resdotold; // betak = resdot_k/resdot_k-1
		
		for (i = 1; i <= imax; ++i) /* p_k Update search direction vector */
			for (j = 1; j <= jmax; ++j)
				if ((mesh(i,j) & C_F))
					p(i,j) = Miplus1*(Res(i+1,j) + Res(i-1,j)) + MiplusNx*(Res(i,j+1) + Res(i,j-1)) + diagMi*Res(i,j) + betak*p(i,j);
//		setBC(p,mesh,1);
		tol = sqrt(resdot/nfluid);
	}
	/** Boundary treatment */
	if((tol < TOL_)||(tol < tol0)){
		residual = tol;
		return true;
	}
	else{
		return false;
	}
}

bool SORSolver::solve_SOR1(StaggeredGrid& grid, Geometry2D& mesh){
	int i,j,nfluid;
	int epsE, epsW, epsN, epsS = 0;
	real rdx2,rdy2;
	real add,beta_2,beta_mod,tol;
	int imax = grid.Nx();
	int jmax = grid.Ny();
	real delx = grid.dx();
	real dely = grid.dy();
	real omg = W_;
	Array &P = grid.p();
	Array &RHS = grid.rhs();
	rdx2 = 1./delx/delx;
	rdy2 = 1./dely/dely;
	beta_2 = -omg/(2.0*(rdx2+rdy2));
	nfluid = mesh.geom_nFluids();
	tol = 1.0;
  real p0 = 0.;
	for (i=1;i<=imax;i++)
		for (j=1;j<=jmax;j++)
 			if ((mesh(i,j) & C_F))
				p0 += P(i,j)*P(i,j);
	p0 = sqrt(p0/nfluid);
	if (p0 < 0.0001)
	 p0 = 1.0;  
                                            /* SOR-iteration */
                                            /*---------------*/
	for (iter_=0;(tol > TOL_*p0)&&(iter_ < iterMax_);++iter_){
		                       /* copy values at external boundary */
		                       /*----------------------------------*/
		setBC(P,mesh,1);
	
		/* relaxation for fluid cells */
		/*----------------------------*/
		for (i=1;i<=imax;i+=1)
		        for (j=1;j<=jmax;j+=1)
			  if ((mesh(i,j) & C_F))
			    P(i,j) = (1.-omg)*P(i,j) - beta_2*((P(i+1,j)+P(i-1,j))*rdx2 +	(P(i,j+1)+P(i,j-1))*rdy2 - RHS(i,j));
	
		/* computation of residual */
		/*-------------------------*/
		if (iter_%rescheckfreq_ == 0){
			tol = 0.0;
			for (i=1;i<=imax;i++)
	      for (j=1;j<=jmax;j++)
					if ((mesh(i,j) & C_F))   
				  /* only fluid cells */
				  /*------------------*/
			    {
						add =  (P(i+1,j)-2*P(i,j)+P(i-1,j))*rdx2+(P(i,j+1)-2*P(i,j)+P(i,j-1))*rdy2-RHS(i,j);
						tol += add*add;
				  }
	
				    tol = sqrt((tol)/nfluid);
			/* convergence? */
			/*--------------*/
		}
	}
	if(tol < TOL_){
		residual = tol;
		return true;
	}
	else{
		return false;
	}
}

real SORSolver::residualnorm(){ return residual; }

void SORSolver::setBC(Array &Q,const int & nGhost){
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

void SORSolver::setBC(Array &Q,Geometry2D& mesh,const int & nGhost){
	int nx = Q.getSize(0);
	int ny = Q.getSize(1);
	int imax = nx - 2;
	int jmax = ny - 2;
	
	for (int i = 1; i <= imax; ++i){
		if(Vbctype_ == 4){ //Periodic BC treatment for pressure equation (N-S)
			Q(i,1) = Q(i,jmax);
			Q(i,0) = Q(i,1);
			Q(i,jmax+1) = Q(i,jmax);
		}
		else{
			Q(i,0) = Q(i,1);
			Q(i,jmax+1) = Q(i,jmax);
		}
	}
	for (int j = 1; j <= jmax; ++j){
		if(Hbctype_ == 4){ //Periodic BC treatment for pressure equation (E-W)
			Q(0,j) = Q(1,j);
			Q(imax+1,j) = Q(imax,j);		
			Q(1,j) = Q(imax,j);
		}
		else{
			Q(0,j) = Q(1,j);
			Q(imax+1,j) = Q(imax,j);
		}
	}
	for (int i=1;i<=imax;++i)
		for (int j=1;j<=jmax;++j)
			if (mesh(i,j) >=B_N && mesh(i,j) <=B_SO) 
				switch (mesh(i,j))
				{
					case B_N:{  Q(i,j) = Q(i,j+1);                 break;}
					case B_O:{  Q(i,j) = Q(i+1,j);                 break;}
					case B_S:{  Q(i,j) = Q(i,j-1);                 break;} 
					case B_W:{  Q(i,j) = Q(i-1,j);                 break;}
					case B_NO:{ Q(i,j) = 0.5*(Q(i,j+1)+Q(i+1,j)); break;}
					case B_SO:{ Q(i,j) = 0.5*(Q(i,j-1)+Q(i+1,j)); break;}
					case B_SW:{ Q(i,j) = 0.5*(Q(i,j-1)+Q(i-1,j)); break;}
					case B_NW:{ Q(i,j) = 0.5*(Q(i,j+1)+Q(i-1,j)); break;}
					default:                                         break;
				}
}
