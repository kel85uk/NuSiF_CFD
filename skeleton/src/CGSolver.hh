#ifndef CG_SOLVER_HH
#define CG_SOLVER_HH

#include "StaggeredGrid.hh"
#include "Geometry2D.hh"
//#include "Solver.hh"

class CGSolver//: public Solver
{
public:
	CGSolver(){;}
   // Constructor to create a CGSolver from a parsed configuration file
   CGSolver ( const FileReader & configuration );

   // Solve the pressure equation on the staggered grid. Returns true if converged
   bool solve( StaggeredGrid & grid );
   
   // solve the pressure equation on the staggered grid with obstacles
   bool solve( StaggeredGrid & grid, Geometry2D& mesh);

   int iterations();
   real residualnorm();    

private:
   unsigned int imax,jmax, iter_, itermax, Hbctype_,Vbctype_,checkfrequency;
   real eps, residual;
   Array<real> r;
   Array<real> Ark;
   Array<real> Residual;
   // copy inner points to boundary
	void setBC(Array<real>& Q,const int & nGhost);
	void setBC(Array<real>& Q,Geometry2D& mesh);
};
#endif //CG_SOLVER_HH
