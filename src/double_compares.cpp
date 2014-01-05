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
#include <double_compares.h>

namespace double_compares
{
	bool d_equal_str(real i, real j, real tol){
		return (std::abs((i-j)/i) <= tol) && (std::abs((i-j)/j) <= tol);
	}
	
	bool d_equal_abs_str(real i, real j, real tol){
		return (std::abs((i-j)) <= tol) && (std::abs((i-j)) <= tol);
	}	
	
	bool d_equal_wk(real i, real j, real tol){
		return (std::abs((i-j)/i) <= tol) || (std::abs((i-j)/j) <= tol);
	}	
	
	bool d_equal_str(real i, real j){
		return (std::abs((i-j)/i) <= 100*mach_err) && (std::abs((i-j)/j) <= 100*mach_err);
	}
	
	bool d_equal_abs_str(real i, real j){
		return (std::abs((i-j)) <= 10*mach_err) && (std::abs((i-j)) <= 10*mach_err);
	}	
	
	bool d_equal_wk(real i, real j){
		return (std::abs((i-j)/i) <= 100*mach_err) || (std::abs((i-j)/j) <= 100*mach_err);
	}
	
	bool d_equal_str_gen(real i, real j){
		if(i<1.0 || j<1.0)
			return d_equal_abs_str(i,j);
		else
			return d_equal_str(i,j);
	}
}
