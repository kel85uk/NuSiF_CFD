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

#include "Array.hh"
#include "FileReader.hh"

#include <iostream>


int main( int argc, char** argv )
{
	FileReader readParam;
	readParam.registerIntParameter("width");
	readParam.registerIntParameter("height");
	readParam.registerIntParameter("x");
	readParam.registerIntParameter("y");
    readParam.registerStringParameter("name");
	readParam.registerRealParameter("initial");
	
	bool res = readParam.readFile ("ArrayTestInput.txt");
	CHECK(res);
	readParam.printParameters();
	int xlength, ylength, xind, yind;
	xlength = readParam.getIntParameter("width");
	ylength = readParam.getIntParameter("height");
	xind = readParam.getIntParameter("x");
	yind = readParam.getIntParameter("y");
	real initial = readParam.getRealParameter("initial");
	
	Array arr(xlength,ylength);
	arr.fill(initial);
	arr(xind,yind) = 2*arr(xind,yind);
	
	arr.print();
    PROGRESS("Creating discrete Laplacian operator using kron product");
	Array Dxx(3,3);
	Array Ixx(3,3);
	Array Dyy(4,4);
	Array Iyy(4,4);
	Array Dxpm1(2);
	Array Dx0(3);
	Array Dy0(4);
	Array Dypm1(3);
	Dxpm1.fill(-1);
	Dypm1.fill(-1);
	Dx0.fill(2);
	Dy0.fill(2);
	Ixx.eye();
	Iyy.eye();
	Array D1(12,12), D2(12,12);
	Dxx.fill(0);
    PROGRESS("Filling in diagonals of the x, and y-sweep matrices");
	Dxx.diag(-1,Dxpm1,1);
	Dxx.diag(1,Dxpm1,1);
	Dxx.diag(0,Dx0,1);
	Dyy.fill(0);
	Dyy.diag(-1,Dypm1,1);
	Dyy.diag(1,Dypm1,1);
	Dyy.diag(0,Dy0,1);
    PROGRESS("Making the kronecker sum for the final discrete operator");
	D1.kron(Dxx,Iyy);
	D2.kron(Ixx,Dyy);
	D1 += D2;
	D1.print2();
	return 0;
}
