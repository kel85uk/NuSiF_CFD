#ifndef STAGGERED_GRID_HH
#define STAGGERED_GRID_HH


#include <Types.hh>
#include <Array.hh>
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

   StaggeredGrid& operator= (StaggeredGrid &rihs);

	 void meshgrid(Array x, Array y, Array &X, Array &Y);
     void setCoord(Array &x, real dx1, int nGhost, bool stagger);
   // Getters / Setters for member variables
   Array & p();//    { return p_;    }
   Array & rhs();//  { return rhs_;  }
   Array & U();
   Array & V();
   Array & F();
   Array & G();
   Array & UDA();

   const Array & p()   const;// { return p_;   }
   const Array & rhs() const;// { return rhs_; }
   const Array & U()   const;
   const Array & V()   const;
   const Array & F()   const;
   const Array & G()   const;
   const Array & UDA() const;          
	 const Array & XC() const;
	 const Array & YC() const;	
	 Array & XFU();
	 Array & YFU();	
	 Array & XFV();
	 Array & YFV();		 

   real dx() const;// { return dx; }
   real dy() const;// { return dy; }
   
   int Nx() const;
   int Ny() const;
   int NGhost() const;

protected:
	 int xSize_;
	 int nGhost_;
	 int ySize_;
   Array p_;   //< pressure field
   Array rhs_; //< right hand side of the pressure equation	
   Array U_;
   Array V_;
   Array F_;
   Array G_;
   Array user_array_;
   Array XGridfU_; // x-coordinates of faces U/F
   Array YGridfU_; // y-coordinates of faces U/F
   Array XGridfV_; // x-coordinates of faces V/G
   Array YGridfV_; // y-coordinates of faces V/G
   Array XGridc_; // x-coordinates of cell
   Array YGridc_; // y-coordinates of cell   
   real dx_;   //< distance between two grid points in x direction
   real dy_;   //< distance between two grid points in y direction
   
friend class FluidSimulator;   
};



#endif //STAGGERED_GRID_HH

