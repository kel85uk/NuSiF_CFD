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


void mat_set(int i, int j, double value, vector<element>& data){
	data.push_back(make_tuple(i,j,value));
	sort(data.begin(),data.end(),mycompare);
}

bool mycompare (const element &lhs, const element &rhs){
	if(get<0>(lhs) != get<0>(rhs)) return get<0>(lhs) < get<0>(rhs);
  else return get<1>(lhs) < get<1>(rhs);
}
