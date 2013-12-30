/*    Copyright (C) 2013  kklloh

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

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
