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
#ifndef MATRIXCOO_HH
#define MATRIXCOO_HH

#include <Array.hh>
#include <utility>
#include <vector>
#include <tuple>
#include <algorithm>

using std::tuple;
using std::vector;
using std::get;
using std::sort;
using std::make_tuple;
using std::cout;
using std::endl;

typedef tuple<int,int,double> element;
//*******************************************************************************************************************
/*!  Matrix Class for COO format
*
*    - all elements should be stored in a contiguous chunk of memory ( no vector<vector> ! )
*/
//*******************************************************************************************************************
class MatrixCOO
{
public:
   // Constructors for 1D,2D and 3D
   MatrixCOO();
   MatrixCOO( Array &dense2D ); // Create sparse matrix from dense 2D Array
	 
	 void mat_set(int i, int j, double value, vector<element>& data) // Manually set an element inside the matrix

   // Depending on your implementation you might need the following:
   ~MatrixCOO();
	 MatrixCOO(const MatrixCOO& s);


   // Access Operators for 1D, 2D and 3D
   inline real & operator () ( int i ,int j );
   
   // for const MatrixCOOs the following access operators are required
   inline const real & operator () ( int i ,int j ) const;
  
   //Arithmetic Operators for MatrixCOOs
   MatrixCOO& operator= (const MatrixCOO &rhs);
   MatrixCOO& operator+= (const MatrixCOO &rhs);
   MatrixCOO& operator-= (const MatrixCOO &rhs);   
   MatrixCOO& operator*= (const MatrixCOO &rhs);
   const MatrixCOO operator+ (const MatrixCOO &other) const;
   const MatrixCOO operator- (const MatrixCOO &other) const;
   const MatrixCOO operator* (const MatrixCOO &other) const;      


	 //Convenient functions
	 void speye(); //Create sparse Identity matrix
   
	 void spdiags(int dia,MatrixCOO& Mat);
	 void spdiags(Array& arr,int dia,MatrixCOO& Mat);
	 
	 Array matvec(Array& Y); //Matvec operation

	 MatrixCOO col(int Ncol);
	 void col(int Ncol, Array &vec);
	 MatrixCOO row(int Nrow);
	 void row(int Nrow, Array &vec);

	 void kron(MatrixCOO& X, MatrixCOO& Y); //Kronecker product

 	 MatrixCOO transpose(); //Transpose matrix

   // return total size of the MatrixCOO
   int getSize() const;

   // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
   // other dimension values are not allowed
   int getSize(int dimension ) const;


   // Print the whole MatrixCOO ( for debugging purposes )
   void print();

private:
	vector<element> data;
	bool mycompare (const mytuple &lhs, const mytuple &rhs);
};


//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================


// Operator() 2D
inline real& MatrixCOO::operator ()(int i,int j)
{
	int index = i + length_x*j;
	return *(arr + index);
}

inline const real & MatrixCOO::operator () ( int i ,int j ) const
{
	int index = i + length_x*j;
	return *(arr+index);
}


#endif //MatrixCOO_HH

