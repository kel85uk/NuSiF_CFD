#include "CGSolver.hh"

CGSolver::CGSolver ( const FileReader & configuration ):
                     imax(configuration.getIntParameter("imax")), jmax(configuration.getIntParameter("jmax")),
                     itermax(configuration.getIntParameter("itermax")), Hbctype_(0), Vbctype_(0),
                     checkfrequency(configuration.getIntParameter("checkfrequency")),
                     eps(configuration.getRealParameter("eps")),
                     r(imax+2,jmax+2), Ark(imax+2,jmax+2), Residual(imax+2,jmax+2)
{
	Hbctype_ = configuration.getIntParameter("boundary_condition_E");
	Vbctype_ = configuration.getIntParameter("boundary_condition_S");
}

bool CGSolver::solve( StaggeredGrid & grid )
{
	unsigned int i,j,numFluid_ = grid.Nx()*grid.Ny();
	real dx_2,dy_2;
	real delx = grid.dx();
	real dely = grid.dy();	
	Array<real> &p_ = grid.p();
	Array<real> &rhs_ = grid.rhs();
	dx_2 = 1./(delx*delx);
	dy_2 = 1./(dely*dely);
	real resid = 1e100;
//	real pMax = p_.maximum();
	unsigned int iterno = 0;

	/** Set variables for CG */
	r.fill(0);
	Ark.fill(0);
	Residual.fill(0);
	real resdot = 0., resdotold, rkdot, betak, alphak = 0.;

	/** Initialize CG (Calculate residual of initial condition) */
	setBC(p_,1);
	for (i=1;i<=imax;i++)
		for (j=1;j<=jmax;j++)
		    {
                        // res = b - Ax
				Residual(i,j) =   (p_(i+1,j) + p_(i-1,j) - 2.0*p_(i,j)) *dx_2
                                        + (p_(i,j+1) + p_(i,j-1) - 2.0*p_(i,j)) *dy_2
                                        - rhs_(i,j) ;                                       
				r(i,j) = Residual(i,j); 		   // p = res
				resdot += Residual(i,j)*Residual(i,j); // resdot_0
			}
	//resdot0 = sqrt(resdot/numFluid_);
// CG iteration
   do {
            iterno++;
		// Start CG iterations
		rkdot = 0.0;

		setBC(p_,1);
		setBC(r,1);
//		setBC(Residual,mesh);
		// r_k Update search direction vector
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j){
					Ark(i,j) = -(r(i+1,j)-2*r(i,j)+r(i-1,j))*dx_2 - (r(i,j+1)-2*r(i,j)+r(i,j-1))*dy_2;
					rkdot += r(i,j)*(Ark(i,j));
				}
		alphak = resdot/rkdot;
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j){
					p_(i,j) = p_(i,j) + alphak*r(i,j);    // Update pressure
				}
		setBC(p_,1);
		setBC(r,1);
//		setBC(Residual,mesh);		
		if (iterno%checkfrequency == 0){
			for (i = 1; i <= imax; ++i)
				for (j = 1; j <= jmax; ++j)
						Residual(i,j) =   (p_(i+1,j) + p_(i-1,j) - 2.0*p_(i,j)) *dx_2
                                        + (p_(i,j+1) + p_(i,j-1) - 2.0*p_(i,j)) *dy_2
                                        - rhs_(i,j) ;
		} 
		else{
			for (i = 1; i <= imax; ++i)
				for (j = 1; j <= jmax; ++j)
						Residual(i,j) = Residual(i,j) - alphak*Ark(i,j); // Update residual		
		}
      if (iterno %100 == 0 || iterno == 2)
		std::cout<<"Iteration no = "<<iterno<< "\tResidual = "<< resid<<"\n";
		resdotold = resdot;
		resdot = 0.;

		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
					resdot += Residual(i,j)*Residual(i,j);
					
		betak = resdot/resdotold;
		
		// r_k Update search direction vector
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
					r(i,j) = Residual(i,j) + betak*r(i,j);
					
//		resid = Residual.amax();
		resid = sqrt(resdot/numFluid_);
	}
	while (resid>=eps && iterno<itermax);   
	residual = resid;
	iter_ = iterno;
	if(resid <= eps){
		return true;
	}
	else{
		return false;
	}
}

bool CGSolver::solve( StaggeredGrid & grid, Geometry2D& mesh )
{
	unsigned int i,j,numFluid_ = mesh.geom_nFluids();
	real dx_2,dy_2;
	real delx = grid.dx();
	real dely = grid.dy();	
	Array<real> &p_ = grid.p();
	Array<real> &rhs_ = grid.rhs();
	dx_2 = 1./(delx*delx);
	dy_2 = 1./(dely*dely);
	real resid = 1e100;
//	real pMax = p_.maximum();
	unsigned int iterno = 0;

	/** Set variables for CG */
	r.fill(0);
	Ark.fill(0);
	Residual.fill(0);
	real resdot = 0., resdotold, rkdot, betak, alphak = 0.;

	/** Initialize CG (Calculate residual of initial condition) */
	setBC(p_,mesh);
	for (i=1;i<=imax;i++)
		for (j=1;j<=jmax;j++)
			if (mesh.geom_isFluid(i,j))   
		    {
                        // res = b - Ax
				Residual(i,j) =   (p_(i+1,j) + p_(i-1,j) - 2.0*p_(i,j)) *dx_2
                                        + (p_(i,j+1) + p_(i,j-1) - 2.0*p_(i,j)) *dy_2
                                        - rhs_(i,j) ;                                       
				r(i,j) = Residual(i,j); 		   // p = res
				resdot += Residual(i,j)*Residual(i,j); // resdot_0
			}
	//resdot0 = sqrt(resdot/numFluid_);
// CG iteration
   do {
            iterno++;
		// Start CG iterations
		rkdot = 0.0;

		setBC(p_,mesh);
		setBC(r,mesh);
//		setBC(Residual,mesh);
		// r_k Update search direction vector
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
				if (mesh.geom_isFluid(i,j)){
					Ark(i,j) = -(r(i+1,j)-2*r(i,j)+r(i-1,j))*dx_2 - (r(i,j+1)-2*r(i,j)+r(i,j-1))*dy_2;
					rkdot += r(i,j)*(Ark(i,j));
				}
		alphak = resdot/rkdot;
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
				if (mesh.geom_isFluid(i,j)){
					p_(i,j) = p_(i,j) + alphak*r(i,j);    // Update pressure
				}
		setBC(p_,mesh);
		setBC(r,mesh);
//		setBC(Residual,mesh);		
		if (iterno%checkfrequency == 0){
			for (i = 1; i <= imax; ++i)
				for (j = 1; j <= jmax; ++j)
					if (mesh.geom_isFluid(i,j))
						Residual(i,j) =   (p_(i+1,j) + p_(i-1,j) - 2.0*p_(i,j)) *dx_2
                                        + (p_(i,j+1) + p_(i,j-1) - 2.0*p_(i,j)) *dy_2
                                        - rhs_(i,j) ;
		} 
		else{
			for (i = 1; i <= imax; ++i)
				for (j = 1; j <= jmax; ++j)
					if (mesh.geom_isFluid(i,j))
						Residual(i,j) = Residual(i,j) - alphak*Ark(i,j); // Update residual		
		}
      if (iterno %100 == 0 || iterno == 2)
		std::cout<<"Iteration no = "<<iterno<< "\tResidual = "<< resid<<"\n";
		resdotold = resdot;
		resdot = 0.;

		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
				if (mesh.geom_isFluid(i,j))
					resdot += Residual(i,j)*Residual(i,j);
					
		betak = resdot/resdotold;
		
		// r_k Update search direction vector
		for (i = 1; i <= imax; ++i)
			for (j = 1; j <= jmax; ++j)
				if (mesh.geom_isFluid(i,j))
					r(i,j) = Residual(i,j) + betak*r(i,j);
					
//		resid = Residual.amax();
		resid = sqrt(resdot/numFluid_);
	}
	while (resid>=eps && iterno<itermax);   
	residual = resid;
	iter_ = iterno;
	if(resid <= eps){
		return true;
	}
	else{
		return false;
	}
}

void CGSolver::setBC(Array<real>& Q,Geometry2D& mesh){
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

void CGSolver::setBC(Array<real> &Q,const int & nGhost){
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

real CGSolver::residualnorm(){ return residual; }
int CGSolver::iterations(){ return iter_; }
