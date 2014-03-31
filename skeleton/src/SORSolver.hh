#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"
#include "Geometry2D.hh"
#include "Solver.hh"


class SORSolver//: public Solver
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
   bool solve( StaggeredGrid & grid, Geometry2D& mesh);
   int iterations();
   real residualnorm();

private:
   int iterMax_,Hbctype_,Vbctype_;
   int iter_, rescheckfreq_;
   real W_;
   real TOL_;
   real residual;
   Array<real> sol_;
   std::string nrmtype_;
	void setBC(Array<real>& Q,const int & nGhost);
	void setBC(Array<real>& Q,Geometry2D& mesh);
};






#endif //SOR_SOLVER_HH




