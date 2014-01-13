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

#include <PSolver.hh>

PSolver::PSolver(){
}
// Constructor to manually create PSolver
PSolver::PSolver ( int solver, int iterMax, real TOL, real W, std::string nrmtype ){
	iterMax_ = iterMax;
	TOL_ = TOL;
	W_ = W;
	nrmtype_ = nrmtype;
	solver_ = solver;
}

// Constructor to create a PSolver from a parsed configuration file
PSolver::PSolver ( const FileReader & configuration ){
	iterMax_ = configuration.getIntParameter("itermax");
	TOL_ = configuration.getRealParameter("eps");
	rescheckfreq_ = configuration.getIntParameter("checkfrequency");
	nrmtype_ = "2"; // Hard coded to comply with SORTest.c
  solver_ = configuration.getIntParameter("Solver");
  if(solver_ == 0){
		W_ = configuration.getRealParameter("omg");
	}
}

int PSolver::iterations(){
	return iter_;
}

MatrixCOO PSolver::mat_assemble( StaggeredGrid& grid, Geometry2D& mesh ){
	int NxG = grid.Nx() + 2;
	int NyG = grid.Ny() + 2;
	real dx = grid.dx();
	real dy = grid.dy();
	real idx2 = 1./(dx*dx);
	real idy2 = 1./(dy*dy);
	MatrixCOO Amat;
	int index = 0;
	for (int i = 0; i<NxG; ++i)
		for (int j = 0; j<NyG; ++j){
			index = i + NxG*j;
			if((i!=0)&&(i!=NxG-1)&&(j!=0)&&(j!=NyG-1)){
				if(mesh(i,j) & C_F){
					Amat.mat_set(index,index,-2*idx2-2*idy2);
					Amat.mat_set(index,index+1,idx2);
					Amat.mat_set(index,index-1,idx2);
					Amat.mat_set(index,index+NxG,idy2);
					Amat.mat_set(index,index-NxG,idy2);
				}
				if (mesh(i,j) >=B_N && mesh(i,j) <=B_SO) 
					switch (mesh(i,j))
					{
						case B_N:{ Amat.mat_set(index,index,1); Amat.mat_set(index,index+NxG,-1); break; }
						case B_O:{ Amat.mat_set(index,index,1); Amat.mat_set(index,index+1,-1); break; }
						case B_S:{ Amat.mat_set(index,index,1); Amat.mat_set(index,index-NxG,-1); break; } 
						case B_W:{ Amat.mat_set(index,index,1); Amat.mat_set(index,index-1,-1); break; }
						case B_NO:{ Amat.mat_set(index,index,1); Amat.mat_set(index,index+NxG,-0.5); Amat.mat_set(index,index+1,-0.5); break; }
						case B_SO:{ Amat.mat_set(index,index,1); Amat.mat_set(index,index-NxG,-0.5); Amat.mat_set(index,index+1,-0.5); break;}
						case B_SW:{ Amat.mat_set(index,index,1); Amat.mat_set(index,index-NxG,-0.5); Amat.mat_set(index,index-1,-0.5); break; }
						case B_NW:{ Amat.mat_set(index,index,1); Amat.mat_set(index,index+NxG,-0.5); Amat.mat_set(index,index-1,-0.5); break; }
						default: break;
					}
			}
			if((i==0)&&(j!=0)&&(j!=NyG-1)){
				Amat.mat_set(index,index,1);
				Amat.mat_set(index,index+1,-1);
			}
			if((i==NxG-1)&&(j!=0)&&(j!=NyG-1)){
				Amat.mat_set(index,index,1);
				Amat.mat_set(index,index-1,-1);
			}
			if((j==0)&&(i!=0)&&(i!=NxG-1)){
				Amat.mat_set(index,index,1);
				Amat.mat_set(index,index+NxG,-1);
			}
			if((j==NyG-1)&&(i!=0)&&(i!=NxG-1)){
				Amat.mat_set(index,index,1);
				Amat.mat_set(index,index-NxG,-1);
			}
			if( ((i==0)&&(j==0))||((i==0)&&(j==NyG-1))||((i==NxG-1)&&(j==0))||((i==NxG-1)&&(j==NyG-1)) ){
				Amat.mat_set(index,index,1);
			}
		}
	Amat_ = Amat;
	return Amat;
}

bool PSolver::solve( StaggeredGrid& grid, Geometry2D& mesh, int solver ){
	Amat_ = mat_assemble( grid,mesh );
	nfluid = mesh.geom_nFluids();
	switch(solver){
		case 0:{return solve_PCG(grid);}
		case 1:{return solve_SORRB(grid);}
		default: return false;
	}
}

bool PSolver::solve_PCG(StaggeredGrid& grid){
	int NxG = grid.p().getSize(0);
	int NyG = grid.p().getSize(1);
	Array SOL, b;
	SOL = SOL.vectorize(grid.p());
	b = b.vectorize(grid.rhs());
	real tol = 1.0;
/** Set variables for CG */
	Array q(SOL);
	Array Res(SOL);
	Res.fill(0);
	q.fill(0);
	Array d(SOL);
	d.fill(0);
	real p0 = sqrt(SOL.dotNC(SOL));
	real dnew, dold, d0, betak, alphak = 0.;
	/** Initialize the CG loop **/
	Res = -b - Amat_.mvmult(SOL); // res = b-Ax
	d = Res;
	dnew = Res.dotNC(Res);
	d0 = dnew;
	std::cout << Res.getSize(0) << std::endl;
	std::cout << SOL.getSize(0) << std::endl;
	std::cout << iter_ << ":" << dnew << std::endl;
	while(iter_ < iterMax_ && dnew > TOL_*d0){
		q = Amat_.mvmult(d);
		alphak = dnew/(d.dotNC(q));
		SOL = SOL + d*alphak;
		if (iter_%50 == 0){
			Res = -b - Amat_.mvmult(SOL);
		}
		else{
			Res = Res - q*alphak;
		}
		dold = dnew;
		dnew = Res.dotNC(Res);
		betak = dnew/dold;
		d = Res + d*betak;
		++iter_;
		std::cout << iter_ << ":" << dnew << std::endl;	
	}
	grid.p() = grid.p().reshape(SOL,NxG,NyG);
	std::cout << iter_ << ":" << dnew << std::endl;	
	return true;
}

bool PSolver::solve_SORRB(StaggeredGrid& grid){
	return true;
}

real PSolver::residualnorm(){ return residual; }
