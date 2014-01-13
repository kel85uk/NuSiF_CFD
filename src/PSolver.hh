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

#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"
#include "Geometry2D.hh"
#include <MatrixCOO.hh>
#include <omp.h>


class PSolver
{
public:
   // Constructor to manually create PSolver
   PSolver();
   PSolver ( int solver, int iterMax, real TOL, real W, std::string nrmtype );

   // Constructor to create a PSolver from a parsed configuration file
   PSolver ( const FileReader & configuration );
   // Solve the pressure equation on the staggered grid
   bool solve( StaggeredGrid& grid, Geometry2D& mesh, int solver );
   MatrixCOO mat_assemble( StaggeredGrid& grid, Geometry2D& mesh );
   // solve the pressure equation on the staggered grid with obstacles
   bool solve_PCG( StaggeredGrid & grid);
   bool solve_SORRB( StaggeredGrid & grid);
   int iterations();
   real residualnorm();
   inline int eps_E(int i, int j, Geometry2D&mesh);
   inline int eps_N(int i, int j, Geometry2D&mesh);
   inline int eps_W(int i, int j, Geometry2D&mesh);
   inline int eps_S(int i, int j, Geometry2D&mesh);

private:
   int iterMax_,Hbctype_,Vbctype_;
   int iter_=0, rescheckfreq_=1, solver_=0;
   int nfluid = 1;
   real W_ = 1.0;
   real TOL_ = 1e-2;
   real residual;
   MatrixCOO Amat_;
   std::string nrmtype_;
};


inline int PSolver::eps_E(int i, int j, Geometry2D&mesh){
	int epsE = (mesh.geom_isFluid(i+1,j))? 1 : 0;
	return epsE;
}
inline int PSolver::eps_S(int i, int j, Geometry2D&mesh){
	int epsS = (mesh.geom_isFluid(i,j-1))? 1 : 0;
	return epsS;
}
inline int PSolver::eps_W(int i, int j, Geometry2D&mesh){
	int epsW = (mesh.geom_isFluid(i-1,j))? 1 : 0;
	return epsW;
}
inline int PSolver::eps_N(int i, int j, Geometry2D&mesh){
	int epsN = (mesh.geom_isFluid(i,j+1))? 1 : 0;
	return epsN;
}


#endif //SOR_SOLVER_HH




