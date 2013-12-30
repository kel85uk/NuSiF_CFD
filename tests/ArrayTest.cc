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
#include "Debug.hh"

#include <iostream>

void copyTest()
{
   const int size = 10;

   // Fill single array...
   Array arr (size);
   for ( int i = 0; i<size; ++i )
      arr(i) = i;
   // check if values where set correctly
   for ( int i = 0; i<size; ++i )
      CHECK( arr(i) == i );

   // create a copy of the array, and check if values are the same
   Array arrCopy (arr);
   CHECK(arrCopy == arr);
   for ( int i = 0; i<size; ++i )
      CHECK( arrCopy(i) == arr(i) );


   // write new values to copy
   for ( int i = 0; i<size; ++i )
      arrCopy(i) = 2*i;

   // and check both arrays again
   for ( int i = 0; i<size; ++i )
      CHECK( arr(i) == i );
   for ( int i = 0; i<size; ++i )
      CHECK( arrCopy(i) == 2*i );

}

void  contiguousMemoryTest()
{
   Array arr ( 10, 5 );

   size_t memDifference =  & ( arr (9,4 ) ) - & ( arr(0,0) );
   CHECK( memDifference == 10*5 - 1 );
}

void dimension_fill_Test()
{
	 int size = 50;
   Array arr4(size);
   arr4.fill(4.5);
   for ( int i = 0; i<size; ++i )
      CHECK( arr4(i) == 4.5 );
      
   Array arr5(50,10,20);
   arr5.fill(1.0);
   CHECK(arr5.getSize(0) == 50);
   CHECK(arr5.getSize(1) == 10);
   CHECK(arr5.getSize(2) == 20);
   CHECK(arr5.getSize() == 50*10*20);
	 for ( int i = 0; i<arr5.getSize(); ++i )
    		CHECK( arr5(i) == 1.0 );
   Array arrvec(2);
   arrvec.fill(-1);
   Array arr7(3,3);
   Array arr8(3,2);
	 for (int i = 0; i < arr7.getSize(); ++i)
	 	arr7(i) = i+1;
	 for (int i = 0; i < arr8.getSize(); ++i)
	 	arr8(i) = 2*(i+1);
	 arr7.col(1,arrvec);
	 arr7.print2();
	 arr8.print2();
   Array arr9(5,5); //Just need to set an array
   arr9.fill(4); 
   arr9 = arr7 * arr8;
   Array arrI(3,2);
   arrI.eye();
   arr8.transpose().print2();  
   arr9 *= arrI.transpose();
   arr9.print2();
   Array arr7cpy(arr7);
   arr7cpy.diag(1,arrvec,1);
   arr7cpy.diag(-1,arrvec,1);
   CHECK(arr7cpy(1,2) == -1);  
   arr7cpy.print2();
   Array arrn1(3);
   arrn1.fill(0);
   Array arrn2(arrn1);
   arrn1(1) = 1;
   arrn2(0) = 1;
//   std::cout << "Pass\n";
   real real9 = arrn1.dot(arrn2);
	 CHECK(real9 == 0);
//	 std::cout << "Pass\n";
	 Array A(2,3), B(3,2);
	 Array A1(2,2), B1(2,2);
	 A(0,0) = 2; A(0,1) = -4; A(0,2) = -3;
	 A(1,0) = 4; A(1,1) = -1; A(1,2) = -2;
	 B(0,0) = 2; B(0,1) = -4;
	 B(1,0) = 2; B(1,1) = -3;
	 B(2,0) = 3; B(2,1) = -1;
	 A1(0,0) = 1; A1(0,1) = -2;
	 A1(1,0) = -1; A1(1,1) = 0;
	 B1(0,0) = 4; B1(0,1) = -3;
	 B1(1,0) = 2; B1(1,1) = 3;
	 Array D2D(4,4), D2D1(4,4), Kres(4,4);
	 D2D.kron(A1,B1);
//	 D2D1.kron(B1,A1);
	 Kres = D2D; // + D2D1;
   Kres.print2();
   Array arr11(3);
   arr11(0) = 1; arr11(1) = 2; arr11(2) = 3;
   std::string norm_val = "2";
	 CHECK(arr11.norm(norm_val) == sqrt(14));
	 std::cout << arr9.getSize(1) << "\t" << arr11.getSize(0) << "\n";
	 arr7 = arr9*arr11.normalize();
	 arr11.normalize().print2();
	 CHECK(arr11.normalize().norm(norm_val) == 1.);
}


int main( )
{
   std::cout << "Copy Test: ";
   copyTest();
   std::cout << "OK" << std::endl;

   std::cout << "Contiguous Memory Test: ";
   contiguousMemoryTest();
   std::cout << "OK" << std::endl;
   
   std::cout << "Dimensions, Fill, and Matrix tools Test: \n";
   dimension_fill_Test();
   std::cout << "OK" << std::endl;

   return 0;
}
