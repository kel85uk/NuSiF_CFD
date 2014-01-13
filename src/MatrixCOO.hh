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
#include <double_compares.h>
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
using double_compares::d_equal_str_gen;

typedef tuple<int,int,real> element;

//*******************************************************************************************************************
/*!  Matrix Class for COO format
*/
//*******************************************************************************************************************
class MatrixCOO
{	 
	 struct myclass {
		bool operator() (const element &lhs, const element &rhs){
			if(get<0>(lhs) != get<0>(rhs)) return get<0>(lhs) < get<0>(rhs);
			else return get<1>(lhs) < get<1>(rhs);
  		}
	 } myobject;

public:
   MatrixCOO();
   MatrixCOO( Array &dense2D ); // Create sparse matrix from dense 2D Array
	 
	 void mat_set(int i, int j, double value); // Manually set an element inside the matrix

	 MatrixCOO(const MatrixCOO& s); //Not implemented yet

   inline real& operator () ( int i ,int j );
   
   inline const real & operator () ( int i ,int j ) const;
  
   //Arithmetic Operators for MatrixCOOs
   MatrixCOO& operator= (const MatrixCOO &rhs);

	 //Convenient functions
	 void speye(int n); //Create sparse Identity matrix
	 
	 Array mvmult(Array& Y); //Matvec operation
	 MatrixCOO diags_sp();
	 MatrixCOO tril();
	 MatrixCOO triu();
 	 MatrixCOO transpose(); //Transpose matrix (Not implemented yet)

   // return total size of the MatrixCOO
   int getSize() const;

   // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
   // other dimension values are not allowed
   int getSize(int dimension ) const;
   
   // returns rows, cols, and values of the sparse matrix (Not implemented yet)
   vector<int> rows() const;
   vector<int> cols() const;
   vector<real> values() const;

   // Print the whole MatrixCOO ( for debugging purposes )
   void print_sp();
   void printd(); //Not implemented yet

private:
	vector<element> data;
	int Nnz_;
};


//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================


// Operator() 2D

inline real& MatrixCOO::operator ()(int i,int j)
{
	static real result = 0;
	for(auto iter = this->data.begin(); iter != this->data.end(); ++iter){
		if((get<0>(*iter) == i) && (get<1>(*iter) == j)){
			return get<2>(*iter);
		}
  }
  return result;
}


inline const real & MatrixCOO::operator () ( int i ,int j ) const
{
	static real result = 0;
	for(auto iter = this->data.begin(); iter != this->data.end(); ++iter){
		if((get<0>(*iter) == i) && (get<1>(*iter) == j)){
			return get<2>(*iter);
		}
  }
  return result;
}


#endif //MatrixCOO_HH

