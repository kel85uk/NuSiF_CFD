#ifndef Geometry_HH
#define Geometry_HH

#include <string>
#include <Types.hh>
#include <iostream>
#include <vector>
#include "Array.hh"
#include "StaggeredGrid.hh"
#include "FileReader.hh"


class Geometry2D{


protected:
	
	std::vector<unsigned int> geometry_mesh_;
	int length_x;
	int length_y;
	int length_t;
	int nBlocks;
	void geom_createRectangle(const FileReader& conf,StaggeredGrid& domain);
	void geom_createCircle(const FileReader& conf,StaggeredGrid& domain);
	
	
public:
	
	Geometry2D();
	Geometry2D( int xSize, int ySize );	
	Geometry2D(const StaggeredGrid& domain);


	void geom_init(const FileReader& conf,StaggeredGrid& domain);
	void geom_init();
	void geom_finalize();
	bool geom_isFluid(int i, int j);
	bool geom_isBoundary(int i, int j);
	void geom_setCellToObstacle(int i, int j);
	void geom_setCellToFluid(int i, int j);
	int geom_nFluids();
	void geom_print();
  inline unsigned int& operator () ( int i ,int j );
   
};   
// Operator() 2D
inline unsigned int& Geometry2D::operator ()(int i,int j)
{
	int index = i + length_x*j;
	return this->geometry_mesh_[index];
}


#endif /* Geometry_HH */
