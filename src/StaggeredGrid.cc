#include <StaggeredGrid.hh>
#include "Debug.hh"

StaggeredGrid::StaggeredGrid(){
}

// Constructors to manually create staggered grid
StaggeredGrid::StaggeredGrid ( int xSize, int ySize, int nGhost, real dx1, real dy1 ){
	xSize_ = xSize;
	ySize_ = ySize;
  dx_ = dx1;
  dy_ = dy1;
  nGhost_ = 1;
	Array temp_p_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
	Array temp_rhs_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
  Array temp_U_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
  Array temp_V_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
  Array temp_F_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
  Array temp_G_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
  Array temp_user_array_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
	p_ = temp_p_;	
	rhs_ = temp_rhs_;
	U_ = temp_U_;
	V_ = temp_V_;
	F_ = temp_F_;
	G_ = temp_G_;
	user_array_ = temp_user_array_;
	Array temp_XGridfH_ (xSize_ + 2*nGhost_, ySize_ + 2*nGhost_); //TODO! Need to correct these values at a later date
	Array temp_YGridfH_ (xSize_ + 2*nGhost_, ySize_ + 2*nGhost_);
	Array temp_XGridfV_ (xSize_ + 2*nGhost_, ySize_ + 2*nGhost_);
	Array temp_YGridfV_ (xSize_ + 2*nGhost_, ySize_ + 2*nGhost_);
	Array temp_xGridc_ (xSize_ + 2*nGhost_);
	Array temp_yGridc_ (ySize_ + 2*nGhost_);
	setCoord(temp_xGridc_,dx_,nGhost_,false);
	setCoord(temp_yGridc_,dy_,nGhost_,false);
	Array temp_XGridc_ (xSize_ + nGhost_,ySize_ + nGhost_);
	Array temp_YGridc_ (xSize_ + nGhost_,ySize_ + nGhost_);
	meshgrid(temp_xGridc_,temp_yGridc_,temp_XGridc_,temp_YGridc_);
	XGridc_ = temp_XGridc_;
	YGridc_ = temp_YGridc_;
	setCoord(temp_xGridc_,dx_,nGhost_,true);
	setCoord(temp_yGridc_,dy_,nGhost_,false);
	meshgrid(temp_xGridc_,temp_yGridc_,temp_XGridfH_,temp_YGridfH_);
	XGridfU_ = temp_XGridfH_;
	YGridfU_ = temp_YGridfH_;
	setCoord(temp_xGridc_,dx_,nGhost_,false);
	setCoord(temp_yGridc_,dy_,nGhost_,true);
	meshgrid(temp_xGridc_,temp_yGridc_,temp_XGridfV_,temp_YGridfV_);
	XGridfV_ = temp_XGridfV_;
	YGridfV_ = temp_YGridfV_;			
	ASSERT(XGridc_.getSize(0) == xSize_ + 2*nGhost_);	
	ASSERT(XGridc_.getSize(1) == ySize_ + 2*nGhost_);
	ASSERT(p_.getSize(0) == xSize_ + 2*nGhost_);
	ASSERT(p_.getSize(1) == ySize_ + 2*nGhost_);	
	ASSERT(rhs_.getSize(0) == xSize_ + 2*nGhost_);
	ASSERT(rhs_.getSize(1) == ySize_ + 2*nGhost_);	
}

// Constructor to create a staggered grid from a parsed configuration file
StaggeredGrid::StaggeredGrid ( const FileReader & configuration ){
	FileReader config_param = configuration;
	real xlength = config_param.getRealParameter("xlength");
	real ylength = config_param.getRealParameter("ylength");	
	xSize_ = config_param.getIntParameter("imax");
	ySize_ = config_param.getIntParameter("jmax");
	nGhost_ = 1; //Hard coded for now to comply with SORTest.c
	dx_ = xlength/(xSize_);
	dy_ = ylength/(ySize_);
	Array temp_p_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
	Array temp_rhs_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
  Array temp_U_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
  Array temp_V_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
  Array temp_F_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
  Array temp_G_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
  Array temp_user_array_ (xSize_ + 2*nGhost_,ySize_ + 2*nGhost_);
	p_ = temp_p_;	
	rhs_ = temp_rhs_;
	U_ = temp_U_;
	V_ = temp_V_;
	F_ = temp_F_;
	G_ = temp_G_;
	user_array_ = temp_user_array_;
	Array temp_XGridfH_ (xSize_ + 2*nGhost_, ySize_ + 2*nGhost_); //TODO! Need to correct these values at a later date
	Array temp_YGridfH_ (xSize_ + 2*nGhost_, ySize_ + 2*nGhost_);
	Array temp_XGridfV_ (xSize_ + 2*nGhost_, ySize_ + 2*nGhost_);
	Array temp_YGridfV_ (xSize_ + 2*nGhost_, ySize_ + 2*nGhost_);
	Array temp_xGridc_ (xSize_ + 2*nGhost_);
	Array temp_yGridc_ (ySize_ + 2*nGhost_);
	setCoord(temp_xGridc_,dx_,nGhost_,false);
	setCoord(temp_yGridc_,dy_,nGhost_,false);
	Array temp_XGridc_ (xSize_ + nGhost_,ySize_ + nGhost_);
	Array temp_YGridc_ (xSize_ + nGhost_,ySize_ + nGhost_);
	meshgrid(temp_xGridc_,temp_yGridc_,temp_XGridc_,temp_YGridc_);
	XGridc_ = temp_XGridc_;
	YGridc_ = temp_YGridc_;
	setCoord(temp_xGridc_,dx_,nGhost_,true);
	setCoord(temp_yGridc_,dy_,nGhost_,false);
	meshgrid(temp_xGridc_,temp_yGridc_,temp_XGridfH_,temp_YGridfH_);
	XGridfU_ = temp_XGridfH_;
	YGridfU_ = temp_YGridfH_;
	setCoord(temp_xGridc_,dx_,nGhost_,false);
	setCoord(temp_yGridc_,dy_,nGhost_,true);
	meshgrid(temp_xGridc_,temp_yGridc_,temp_XGridfV_,temp_YGridfV_);
	XGridfV_ = temp_XGridfV_;
	YGridfV_ = temp_YGridfV_;			
	ASSERT(XGridc_.getSize(0) == xSize_ + 2*nGhost_);	
	ASSERT(XGridc_.getSize(1) == ySize_ + 2*nGhost_);
	ASSERT(p_.getSize(0) == xSize_ + 2*nGhost_);
	ASSERT(p_.getSize(1) == ySize_ + 2*nGhost_);	
	ASSERT(rhs_.getSize(0) == xSize_ + 2*nGhost_);
	ASSERT(rhs_.getSize(1) == ySize_ + 2*nGhost_);		
}

// Getters / Setters for member variables
Array & StaggeredGrid::p()    { return p_;    }
Array & StaggeredGrid::rhs()  { return rhs_;  }
Array & StaggeredGrid::U()    { return U_;    }
Array & StaggeredGrid::V()  { return V_;  }
Array & StaggeredGrid::F()    { return F_;    }
Array & StaggeredGrid::G()  { return G_;  }
Array & StaggeredGrid::UDA() {return user_array_;}

const Array & StaggeredGrid::p()   const { return p_;   }
const Array & StaggeredGrid::rhs() const { return rhs_; }
const Array & StaggeredGrid::U()   const { return U_;   }
const Array & StaggeredGrid::V() const { return V_; }
const Array & StaggeredGrid::F()   const { return F_;   }
const Array & StaggeredGrid::G() const { return G_; }
const Array & StaggeredGrid::UDA() const {return user_array_;}

StaggeredGrid& StaggeredGrid::operator = (StaggeredGrid &rihs)
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

void StaggeredGrid::meshgrid(Array x, Array y, Array &X, Array &Y) // Similar to matlab's implementation but X(i,j),Y(i,j) gives the coordinates of point i,j
{
	ASSERT((x.getDimen() == 1)&&(y.getDimen() == 1));
  int nx = x.getSize();
  int ny = y.getSize();
  Array X_temp(nx,ny);
  Array Y_temp(nx,ny);
  for (int ii = 0; ii < ny; ++ii)
  		for (int jj = 0; jj < nx; ++jj)
	    X_temp(jj,ii) = x(jj);
  for (int ii = 0; ii < nx; ++ii)
  		for (int jj = 0; jj < ny; ++jj)
	    Y_temp(ii,jj) = y(jj);
	X = X_temp;
	Y = Y_temp;	    
}

void StaggeredGrid::setCoord(Array &x, real dx1, int nGhost,bool stagger) //TODO = Generalized to include U,V staggering
{
	ASSERT(x.getDimen()==1);
	if (!stagger){
		for (int i = 0; i < x.getSize(); ++i)
            x(i) = dx1*(i - nGhost) + dx1/2;
	}
	if (stagger){
		for (int i = 0; i < x.getSize(); ++i)
            x(i) = dx1*i;
	}
}

const Array & StaggeredGrid::XC() const { return XGridc_; }
const Array & StaggeredGrid::YC() const { return YGridc_; }
Array & StaggeredGrid::XFU() { return XGridfU_; }
Array & StaggeredGrid::YFU() { return YGridfU_; }
Array & StaggeredGrid::XFV() { return XGridfV_; }
Array & StaggeredGrid::YFV() { return YGridfV_; }
real StaggeredGrid::dx() const { return dx_; }
real StaggeredGrid::dy() const { return dy_; }
int StaggeredGrid::Nx() const { return xSize_; }
int StaggeredGrid::Ny() const { return ySize_; }
int StaggeredGrid::NGhost() const { return nGhost_; }
