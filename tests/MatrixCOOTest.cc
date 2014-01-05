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
#include "MatrixCOO.hh"


int main( )
{
	MatrixCOO sparseTest, speyeTest;
	Array	testmat(3,3), testvec(3), result(3);
	testmat(1,2) = -4.0;
	testmat(2,1) = 4.0;
	for (int i = 0; i < testvec.getSize(); ++i)
		testvec(i) = i+1;
		
	MatrixCOO sparsearrTest(testmat);
	sparseTest.mat_set(1,1,3);
	sparseTest.mat_set(0,1,2);
	sparseTest.mat_set(1,2,3);
	sparseTest.mat_set(0,0,6);
	sparseTest.mat_set(1,0,2);
	sparseTest.print_sp();
	speyeTest.speye(3);
	speyeTest.print_sp();
	sparsearrTest.print_sp();
	testmat.print2();
	result = speyeTest.mvmult(testvec);
	result.print();
	return 0;
}
