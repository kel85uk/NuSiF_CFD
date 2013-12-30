#include "Geometry2D.hh"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

Geometry2D::Geometry2D(){
	geometry_mesh_.resize(1,0);
}


Geometry2D::Geometry2D( int xSize, int ySize )
{
	length_x = xSize;
	length_y = ySize;
	length_t = length_x*length_y;
	geometry_mesh_.resize(length_t);
}

Geometry2D::Geometry2D(const StaggeredGrid& domain){
	length_x = domain.Nx() + 2*domain.NGhost();
	length_y = domain.Ny() + 2*domain.NGhost();
	length_t = length_x*length_y;
	geometry_mesh_.resize(length_t);
}

void Geometry2D::geom_init(const FileReader& conf,StaggeredGrid& domain){
	int i,j;
	int imax = conf.getIntParameter("imax");
	int jmax = conf.getIntParameter("jmax");

	for (i = 0; i<= imax + 1; ++i){
		geom_setCellToObstacle(i,0);
		geom_setCellToObstacle(i,jmax+1);
	}
	for (j = 1; j<= jmax; ++j){
		geom_setCellToObstacle(0,j);
		geom_setCellToObstacle(imax+1,j);
	}

	for(i=1;i<=imax;++i)
		for(j=1;j<=jmax;++j)
			geom_setCellToFluid(i,j);
			
	// Special problems
  std::string specprob1 = "backwardstep";
  std::string specprob2 = "cylinder2D";                   
  if(conf.getStringParameter("SpecProb") == specprob1 ){
		geom_createRectangle(conf,domain);
	}
  else if(conf.getStringParameter("SpecProb") == specprob2 ){
		geom_createCircle(conf,domain);
	}	
}

void Geometry2D::geom_createRectangle(const FileReader& conf,StaggeredGrid& domain){
	int i,j;
	int imax = conf.getIntParameter("imax");
	int jmax = conf.getIntParameter("jmax");
	real Xmin, Ymin, Xmax, Ymax;
	Xmin = conf.getRealParameter("RectangleX1");
	Ymin = conf.getRealParameter("RectangleY1");
	Xmax = conf.getRealParameter("RectangleX2");
	Ymax = conf.getRealParameter("RectangleY2");
	for (i=1;i<=imax;++i)
		for (j=1;j<=jmax;++j)
			if ((domain.XC()(i,j) <= Xmax)&&(domain.XC()(i,j) >= Xmin)&&(domain.YC()(i,j) <= Ymax)&&(domain.YC()(i,j) >= Ymin))
				geom_setCellToObstacle(i,j);
}

void Geometry2D::geom_createCircle(const FileReader& conf,StaggeredGrid& domain){
	int i,j;
	int imax = conf.getIntParameter("imax");
	int jmax = conf.getIntParameter("jmax");
	real circcentX, circcentY, radius;
  circcentX = conf.getRealParameter("CircleX");
  circcentY = conf.getRealParameter("CircleY");
  radius = conf.getRealParameter("CircleR");
	for (i=1;i<=imax;++i)
		for (j=1;j<=jmax;++j)
			if ((domain.XC()(i,j) - circcentX)*(domain.XC()(i,j) - circcentX) + (domain.YC()(i,j) - circcentY)*(domain.YC()(i,j) - circcentY) <= radius*radius)
				geom_setCellToObstacle(i,j);
}

void Geometry2D::geom_init(){
	int i,j;
	int imax = length_x - 2;
	int jmax = length_y - 2;
	// Initialize outer boundaries to C_B
	for (i = 0; i<= imax + 1; ++i){
		geometry_mesh_[i] = C_B;
		geometry_mesh_[i + (jmax + 1)*length_x] = C_B;
	}
	for (j = 1; j<= jmax; ++j){
		geometry_mesh_[length_x*j] = C_B;
		geometry_mesh_[imax+1 + j*length_x] = C_B;
	}
                   /* all inner cells fluid cells */
                   /*-----------------------------*/
	for(i=1;i<=imax;++i)
		for(j=1;j<=jmax;++j)
			geometry_mesh_[i + length_x*j] = C_F;
}

void Geometry2D::geom_finalize(){
                    /* FLAGs for boundary cells (MUST BE CALLED after INIT) */
                    /*--------------------------*/
	int i,j;
	int imax = length_x - 2;
	int jmax = length_y - 2;                    
  nBlocks = 0;
  for(i=1;i<=imax;i++)
     for(j=1;j<=jmax;j++){
        if (!(geometry_mesh_[i + length_x*j] & C_F))
					nBlocks++;
					geometry_mesh_[i + length_x*j] += ((geometry_mesh_[i-1 + length_x*j] & C_F)*B_W + (geometry_mesh_[i+1 + length_x*j] & C_F)*B_O + (geometry_mesh_[i + length_x*(j-1)] & C_F)*B_S + (geometry_mesh_[i + length_x*(j+1)] & C_F)*B_N)/C_F;
        switch (geometry_mesh_[i + length_x*j]){
           case 0x0003:
           case 0x0007:
           case 0x000b:
           case 0x000c:
           case 0x000d:
           case 0x000e:
           case 0x000f:{           
                     printf("Illegal obstacle cell [%d][%d]\n",i,j);
                     exit(0);
                    		 }  
	 			}
    }
}

bool Geometry2D::geom_isFluid(int i, int j){
	return geometry_mesh_[i + length_x*j] & C_F;
}

bool Geometry2D::geom_isBoundary(int i, int j){
	return !(geometry_mesh_[i + length_x*j] & C_F);
}

int Geometry2D::geom_nFluids(){
	return (length_x-2)*(length_y-2) - nBlocks;
}

void Geometry2D::geom_setCellToObstacle(int i, int j){
	geometry_mesh_[i + length_x*j] = C_B;
}

void Geometry2D::geom_setCellToFluid(int i, int j){
	geometry_mesh_[i + length_x*j] = C_F;
}

void Geometry2D::geom_print(){
	int i,j;
	int imax = length_x - 2;
	int jmax = length_y - 2;
		                   /* Printing the geometry of the fluid domain */
		                   /*-------------------------------------------*/
	std::cout << "\nGeometry of the fluid domain:\n\n";
	for(j = jmax+1; j >= 0; --j)
		{
		 for(i=0;i<=imax+1;++i)
		    if (!(geometry_mesh_[i + length_x*j] & C_F))
		       printf("0");
		    else      
		       printf("*");                                    
		 printf ("\n");
		}
	printf ("\n");
	printf ("\n");	
}
