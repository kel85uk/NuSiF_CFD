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

//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================
MatrixCOO::MatrixCOO(){
}

MatrixCOO::MatrixCOO(Array& dense2D){
	length_x = dense2D.getSize(0);
	length_y = dense2D.getSize(1);
	Nnz_ = 0;
	for (int i = 0; i < length_x; ++i)
		for (int j = 0; j < length_y; ++j)
			if (!d_equal_str_gen(dense2D(i,j),0.)){
				MatrixCOO::mat_set(i,j,dense2D(i,j));
			}
  for(auto iter = this->data.begin(); iter != this->data.end(); ++iter){
  		rows_.push_back(get<0>(*iter));
  		cols_.push_back(get<1>(*iter));
  		values_.push_back(get<2>(*iter));
  }
}

//===================================================================================================================
//
//  Essential Functions to make life easier
//
//===================================================================================================================
void MatrixCOO::mat_set(int i, int j, double value){
	this->data.push_back(make_tuple(i,j,value));
	sort(this->data.begin(),this->data.end(),myobject);
	++this->Nnz_;
}

void MatrixCOO::speye(int n){
	for (int i = 0; i < n; ++i)
		this->data.push_back(make_tuple(i,i,1.0));	
}

void MatrixCOO::print_sp(){
  for(auto iter = this->data.begin(); iter != this->data.end(); ++iter){
  		cout << "[" << get<0>(*iter) << "," << get<1>(*iter) << "]\t" << get<2>(*iter) << endl;
  }
}
