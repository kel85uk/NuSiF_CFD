#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"
#include "Geometry2D.hh"


class SORSolver
{
public:
   // Constructor to manually create SORSolver
   SORSolver();
   SORSolver ( int iterMax, real TOL, real W, std::string nrmtype );

   // Constructor to create a SORSolver from a parsed configuration file
   SORSolver ( const FileReader & configuration );
   // solve the pressure equation on the staggered grid
   bool solve( StaggeredGrid & grid );
   // solve the pressure equation on the staggered grid with obstacles
   bool solve_CG( StaggeredGrid & grid, Geometry2D& mesh);
   bool solve_PCG( StaggeredGrid & grid, Geometry2D& mesh);
   bool solve_SOR1( StaggeredGrid & grid, Geometry2D& mesh);
   int iterations();
   real residualnorm();
   inline int eps_E(int i, int j, Geometry2D&mesh);
   inline int eps_N(int i, int j, Geometry2D&mesh);
   inline int eps_W(int i, int j, Geometry2D&mesh);
   inline int eps_S(int i, int j, Geometry2D&mesh);

private:
   int iterMax_,Hbctype_,Vbctype_;
   int iter_, rescheckfreq_;
   real W_;
   real TOL_;
   real residual;
   Array sol_;
   std::string nrmtype_;
	 void setBC(Array &Q,const int & nGhost);
	 void setBC(Array &Q,Geometry2D &mesh,const int & nGhost);
};


inline int SORSolver::eps_E(int i, int j, Geometry2D&mesh){
	int epsE = (mesh.geom_isFluid(i+1,j))? 1 : 0;
	return epsE;
}
inline int SORSolver::eps_S(int i, int j, Geometry2D&mesh){
	int epsS = (mesh.geom_isFluid(i,j-1))? 1 : 0;
	return epsS;
}
inline int SORSolver::eps_W(int i, int j, Geometry2D&mesh){
	int epsW = (mesh.geom_isFluid(i-1,j))? 1 : 0;
	return epsW;
}
inline int SORSolver::eps_N(int i, int j, Geometry2D&mesh){
	int epsN = (mesh.geom_isFluid(i,j+1))? 1 : 0;
	return epsN;
}


#endif //SOR_SOLVER_HH




