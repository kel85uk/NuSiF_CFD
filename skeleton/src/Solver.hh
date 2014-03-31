#ifndef SOLVER_HH
#define SOLVER_HH

#include "StaggeredGrid.hh"
#include "Geometry2D.hh"

class Solver
{
public:

   // Constructor to create a Solver
   Solver () {}

   // Solve the pressure equation on the staggered grid. Returns true if converged
   virtual bool solve( StaggeredGrid & grid ) { return false;}
   virtual bool solve( StaggeredGrid & grid, Geometry2D& mesh ) { return false;}
   virtual int iterations();
   virtual real residualnorm();
};
#endif //SOLVER_HH
