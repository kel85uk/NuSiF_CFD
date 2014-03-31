#ifndef STAGGERED_GRID_HH
#define STAGGERED_GRID_HH


#include <Types.hh>
#include <Array_new.hh>
#include <FileReader.hh>

//*******************************************************************************************************************
/*! Class for storing all arrays required for the simulation
*
* For now it only contains an array for the pressure and another array for the
* right hand side of the pressure equation.
* In following assignments this will be extended and also contain
* arrays for x/y velocity and some arrays for intermediate values.
*
* Feel free to add member functions or variables to this class, but please don't
* delete or rename any existing functions, since future skeletons rely on these functions.
*
*/
//*******************************************************************************************************************
class StaggeredGrid
{
public:
	// Default constructor
	StaggeredGrid (); // TODO implement!
   // Constructors to manually create staggered grid
   StaggeredGrid ( int xSize, int ySize, int nGhost, real dx1, real dy1 ); // TODO implement!

   // Constructor to create a staggered grid from a parsed configuration file
   StaggeredGrid ( const FileReader & configuration );       // TODO implement!

   StaggeredGrid& operator= (const StaggeredGrid &rihs);

	void meshgrid(Array<real> x, Array<real> y, Array<real> &X, Array<real> &Y);
	void setCoord(Array<real> &x, real dx1, int nGhost, bool stagger);
   // Getters / Setters for member variables
   Array<real> & p();//    { return p_;    }
   Array<real> & rhs();//  { return rhs_;  }
   Array<real> & U();
   Array<real> & V();
   Array<real> & F();
   Array<real> & G();
   Array<real> & UDA();
//	Array<unsigned int> & Geom();

   const Array<real> & p()   const;// { return p_;   }
   const Array<real> & rhs() const;// { return rhs_; }
   const Array<real> & U()   const;
   const Array<real> & V()   const;
   const Array<real> & F()   const;
   const Array<real> & G()   const;
   const Array<real> & UDA() const;          
	const Array<real> & XC() const;
	const Array<real> & YC() const;	
	const Array<real> & XFU() const;
	const Array<real> & YFU() const;	
	const Array<real> & XFV() const;
	const Array<real> & YFV() const;
//	const Array<unsigned int> & Geom() const;

   real dx() const;// { return dx; }
   real dy() const;// { return dy; }
   
   int Nx() const;
   int Ny() const;
   int NGhost() const;
   /*
   inline real& u( int ii, int jj, Direction dir);
   inline real& v( int ii, int jj, Direction dir);
   inline real& f( int ii, int jj, Direction dir);
   inline real& g( int ii, int jj, Direction dir);
   inline real& p( int ii, int jj, Direction dir);   */

protected:
	int xSize_;
	int nGhost_ = 1;
	int ySize_;
	real dx_;   //< distance between two grid points in x direction
	real dy_;   //< distance between two grid points in y direction	
	Array<real> p_;   //< pressure field
	Array<real> rhs_; //< right hand side of the pressure equation	
	Array<real> U_;
	Array<real> V_;
	Array<real> F_;
	Array<real> G_;
	Array<real> user_array_;
	Array<real> XGridc_; // x-coordinates of cell
	Array<real> YGridc_; // y-coordinates of cell 	
	Array<real> XGridfU_; // x-coordinates of faces U/F
	Array<real> YGridfU_; // y-coordinates of faces U/F
	Array<real> XGridfV_; // x-coordinates of faces V/G
	Array<real> YGridfV_; // y-coordinates of faces V/G
//	Array<unsigned int> Geom_; // Internal geometry array for staggered grid   
};
/*
inline real& StaggeredGrid::u( int x, int y, Direction dir)
{
   switch(dir)
   {
      case NORTH: if (!(Geom_(x,y+1) || Geom_(x+1,y+1)))
                      U_(x,y+1) = -U_(x,y);
                  else if (Geom_(x,y+1) && Geom_(x+1,y+1))
                      return U_(x,y+1);
                  else
                      U_(x,y+1) = 0;
                  return U_(x,y+1);
      case SOUTH:  if (!(Geom_(x,y-1) || Geom_(x+1,y-1)))
                      U_(x,y-1) = -U_(x,y);
                  else if (Geom_(x,y-1) && Geom_(x+1,y-1))
                      return U_(x,y-1);
                  else
                      U_(x,y-1) = 0;
                  return U_(x,y+1);
      case EAST:  if (Geom_(x+1,y))
                     return u(x+1,y,CENTER);
                  else
                     U_(x+1,y) = 0.0;
                  return U_(x+1,y);
      case WEST:  if (!Geom_(x-1,y))
                     U_(x-1,y) = 0.0;
                  return U_(x-1,y);
      default:    if (Geom_(x+1,y))
                     return U_(x,y);
                  else 
                     U_(x,y) = 0.0;
                  return U_(x,y);
   }
}

inline real& StaggeredGrid::v(int x, int y, Direction dir)
{
   switch(dir)
   {
      case NORTH: if (Geom_(x,y+1))
                     return v(x,y+1,CENTER);
                  else
                     V_(x,y+1)=0.0;
                     return V_(x,y+1);
      case SOUTH: if (!Geom_(x,y-1))
                     V_(x,y-1)=0.0;
                  return V_(x,y-1);
      case EAST: if (!(Geom_(x+1,y) || Geom_(x+1,y+1)))
                      V_(x+1,y) = -V_(x,y);
                  else if (Geom_(x+1,y) && Geom_(x+1,y+1))
                      return V_(x+1,y);
                  else
                      V_(x+1,y) = 0;
                  return V_(x+1,y);
      case WEST: if (!(Geom_(x-1,y) || Geom_(x-1,y+1)))
                      V_(x-1,y) = -V_(x,y);
                  else if (Geom_(x-1,y) && Geom_(x-1,y+1))
                      return V_(x-1,y);
                  else
                      V_(x-1,y) = 0;
                  return V_(x-1,y);
      default:    if (!Geom_(x,y+1))
                     V_(x,y) = 0.0;
                  return V_(x,y);
   }
}

inline real& StaggeredGrid::f( int x, int y, Direction dir)
{
   switch(dir)
   {
      case CENTER: if (Geom_(x+1,y))
                     return F_(x,y);
                  else
                     return U_(x,y);
      case WEST:  if (Geom_(x-1,y))
                     return F_(x-1,y);
                  else
                     return U_(x-1,y);
      default:    return F_(x,y);
   }
}
inline real& StaggeredGrid::g(int x, int y, Direction dir)
{
   switch(dir)
   {
      case CENTER: if (Geom_(x,y+1))
                     return G_(x,y);
                  else
                     return V_(x,y);
      case SOUTH: if (Geom_(x,y-1))
                     return G_(x,y-1);
                  else
                     return V_(x,y-1);
      default:    return G_(x,y);
   }
}

inline real& StaggeredGrid::p(int x, int y, Direction dir)
{
   switch(dir)
   {
      case NORTH: if (Geom_(x,y+1))
                     return p_(x,y+1);
                  else
                     return p_(x,y);
      case SOUTH: if (Geom_(x,y-1))
                     return p_(x,y-1);
                  else
                     return p_(x,y);
      case EAST:  if (Geom_(x+1,y))
                     return p_(x+1,y);
                  else
                     return p_(x,y);
      case WEST:  if (Geom_(x-1,y))
                     return p_(x-1,y);
                  else
                     return p_(x,y);
      default:    return p_(x,y);
   }
}
*/
#endif //STAGGERED_GRID_HH

