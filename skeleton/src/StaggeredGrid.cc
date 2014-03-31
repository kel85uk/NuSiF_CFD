#include <StaggeredGrid.hh>
#include "Debug.hh"

StaggeredGrid::StaggeredGrid(){
}

// Constructors to manually create staggered grid
StaggeredGrid::StaggeredGrid ( int xSize, int ySize, int nGhost, real dx1, real dy1 ):
					xSize_(xSize), ySize_(ySize), dx_(dx1), dy_(dy1),
					p_(xSize_ + 2,ySize_ + 2), rhs_(xSize_ + 2,ySize_ + 2),
					U_(xSize_ + 2,ySize_ + 2), V_(xSize_ + 2,ySize_ + 2),
					F_(xSize_ + 2,ySize_ + 2), G_(xSize_ + 2,ySize_ + 2), user_array_(xSize_ + 2,ySize_ + 2)
//					Geom_(xSize_ + 2,ySize_ + 2)
{
	Array<real> temp_xGridc_ (xSize_ + 2*nGhost_);
	Array<real> temp_yGridc_ (ySize_ + 2*nGhost_);
	setCoord(temp_xGridc_,dx_,nGhost_,false);
	setCoord(temp_yGridc_,dy_,nGhost_,false);
	meshgrid(temp_xGridc_,temp_yGridc_,XGridc_,YGridc_);
	setCoord(temp_xGridc_,dx_,nGhost_,true);
	setCoord(temp_yGridc_,dy_,nGhost_,false);
	meshgrid(temp_xGridc_,temp_yGridc_,XGridfU_,YGridfU_);
	setCoord(temp_xGridc_,dx_,nGhost_,false);
	setCoord(temp_yGridc_,dy_,nGhost_,true);
	meshgrid(temp_xGridc_,temp_yGridc_,XGridfV_,YGridfV_);			
	ASSERT(XGridc_.getSize(0) == xSize_ + 2*nGhost_);	
	ASSERT(XGridc_.getSize(1) == ySize_ + 2*nGhost_);
	ASSERT(p_.getSize(0) == xSize_ + 2*nGhost_);
	ASSERT(p_.getSize(1) == ySize_ + 2*nGhost_);	
	ASSERT(rhs_.getSize(0) == xSize_ + 2*nGhost_);
	ASSERT(rhs_.getSize(1) == ySize_ + 2*nGhost_);	
}

// Constructor to create a staggered grid from a parsed configuration file
StaggeredGrid::StaggeredGrid ( const FileReader & config ):
					xSize_(config.getIntParameter("imax")), ySize_(config.getIntParameter("jmax")),
					p_(xSize_ + 2,ySize_ + 2), rhs_(xSize_ + 2,ySize_ + 2),
					U_(xSize_ + 2,ySize_ + 2), V_(xSize_ + 2,ySize_ + 2),
					F_(xSize_ + 2,ySize_ + 2), G_(xSize_ + 2,ySize_ + 2), user_array_(xSize_ + 2,ySize_ + 2),
					XGridc_(xSize_ + 2,ySize_ + 2), YGridc_(xSize_ + 2,ySize_ + 2),
					XGridfU_(xSize_ + 2,ySize_ + 2), YGridfU_(xSize_ + 2,ySize_ + 2),
					XGridfV_(xSize_ + 2,ySize_ + 2), YGridfV_(xSize_ + 2,ySize_ + 2)
//					Geom_(xSize_ + 2,ySize_ + 2)
{
	real xlength = config.getRealParameter("xlength");
	real ylength = config.getRealParameter("ylength");
	dx_ = xlength/(xSize_);
	dy_ = ylength/(ySize_);
	Array<real> temp_xGridc_ (xSize_ + 2*nGhost_);
	Array<real> temp_yGridc_ (ySize_ + 2*nGhost_);
	setCoord(temp_xGridc_,dx_,nGhost_,false);
	setCoord(temp_yGridc_,dy_,nGhost_,false);
	meshgrid(temp_xGridc_,temp_yGridc_,XGridc_,YGridc_);
	setCoord(temp_xGridc_,dx_,nGhost_,true);
	setCoord(temp_yGridc_,dy_,nGhost_,false);
	meshgrid(temp_xGridc_,temp_yGridc_,XGridfU_,YGridfU_);
	setCoord(temp_xGridc_,dx_,nGhost_,false);
	setCoord(temp_yGridc_,dy_,nGhost_,true);
	meshgrid(temp_xGridc_,temp_yGridc_,XGridfV_,YGridfV_);			
	ASSERT(XGridc_.getSize(0) == xSize_ + 2*nGhost_);	
	ASSERT(XGridc_.getSize(1) == ySize_ + 2*nGhost_);
	ASSERT(p_.getSize(0) == xSize_ + 2*nGhost_);
	ASSERT(p_.getSize(1) == ySize_ + 2*nGhost_);	
	ASSERT(rhs_.getSize(0) == xSize_ + 2*nGhost_);
	ASSERT(rhs_.getSize(1) == ySize_ + 2*nGhost_);		
}

StaggeredGrid& StaggeredGrid::operator = (const StaggeredGrid &rihs)
{
	if(this != &rihs)
	{
		xSize_ = rihs.Nx();
		nGhost_ = rihs.NGhost();
		ySize_ = rihs.Ny();
		p_ = rihs.p();   //< pressure field
		rhs_ = rihs.rhs(); //< right hand side of the pressure equation	
		U_ = rihs.U();
		V_ = rihs.V();
		F_ = rihs.F();
		G_ = rihs.G();
		user_array_ = rihs.UDA();
//		Geom_ = rihs.Geom();
		XGridfU_ = rihs.XFU(); // x-coordinates of faces U/F
		YGridfU_ = rihs.YFU(); // y-coordinates of faces U/F
		XGridfV_ = rihs.XFV(); // x-coordinates of faces V/G
		YGridfV_ = rihs.YFV(); // y-coordinates of faces V/G
		XGridc_ = rihs.XC(); // x-coordinates of cell
		YGridc_ = rihs.YC(); // y-coordinates of cell   
		dx_ = rihs.dx();   //< distance between two grid points in x direction
		dy_ = rihs.dy();   //< distance between two grid points in y direction		
	}
	return *this;
}

void StaggeredGrid::meshgrid(Array<real> x, Array<real> y, Array<real> &X, Array<real> &Y) // Similar to matlab's implementation but X(i,j),Y(i,j) gives the coordinates of point i,j
{
	ASSERT((x.getDimen() == 0)&&(y.getDimen() == 0));
	int nx = x.getSize();
	int ny = y.getSize();
	Array<real> X_temp(nx,ny);
	Array<real> Y_temp(nx,ny);
	for (int ii = 0; ii < ny; ++ii)
		for (int jj = 0; jj < nx; ++jj)
		 X_temp(jj,ii) = x(jj);
	for (int ii = 0; ii < nx; ++ii)
		for (int jj = 0; jj < ny; ++jj)
		 Y_temp(ii,jj) = y(jj);
	X = X_temp;
	Y = Y_temp;	    
}

void StaggeredGrid::setCoord(Array<real> &x, real dx1, int nGhost,bool stagger) //TODO = Generalized to include U,V staggering
{
	ASSERT(x.getDimen()==0);
	if (!stagger){
		for (int i = 0; i < x.getSize(); ++i)
            x(i) = dx1*(i - nGhost) + dx1/2;
	}
	if (stagger){
		for (int i = 0; i < x.getSize(); ++i)
            x(i) = dx1*i;
	}
}

// Getters / Setters for member variables
Array<real> & StaggeredGrid::p()    { return p_;    }
Array<real> & StaggeredGrid::rhs()  { return rhs_;  }
Array<real> & StaggeredGrid::U()    { return U_;    }
Array<real> & StaggeredGrid::V()  { return V_;  }
Array<real> & StaggeredGrid::F()    { return F_;    }
Array<real> & StaggeredGrid::G()  { return G_;  }
Array<real> & StaggeredGrid::UDA() {return user_array_;}
//Array<unsigned int> & StaggeredGrid::Geom() {return Geom_;}

const Array<real> & StaggeredGrid::p()   const { return p_;   }
const Array<real> & StaggeredGrid::rhs() const { return rhs_; }
const Array<real> & StaggeredGrid::U()   const { return U_;   }
const Array<real> & StaggeredGrid::V() const { return V_; }
const Array<real> & StaggeredGrid::F()   const { return F_;   }
const Array<real> & StaggeredGrid::G() const { return G_; }
const Array<real> & StaggeredGrid::UDA() const {return user_array_;}

//const Array<unsigned int> & StaggeredGrid::Geom() const {return Geom_;}
const Array<real> & StaggeredGrid::XC() const { return XGridc_; }
const Array<real> & StaggeredGrid::YC() const { return YGridc_; }
const Array<real> & StaggeredGrid::XFU() const { return XGridfU_; }
const Array<real> & StaggeredGrid::YFU() const { return YGridfU_; }
const Array<real> & StaggeredGrid::XFV() const { return XGridfV_; }
const Array<real> & StaggeredGrid::YFV() const { return YGridfV_; }
real StaggeredGrid::dx() const { return dx_; }
real StaggeredGrid::dy() const { return dy_; }
int StaggeredGrid::Nx() const { return xSize_; }
int StaggeredGrid::Ny() const { return ySize_; }
int StaggeredGrid::NGhost() const { return nGhost_; }
