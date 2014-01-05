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
/** Functions for comparing two doubles (strong and weak form): Default tolerance = 100*(Machine double precision) */

#ifndef DOUBLE_COMPARES
#define DOUBLE_COMPARES

#define mach_err (std::numeric_limits<real>::epsilon())

#include <limits>
#include <cmath>
#include <Types.hh>

namespace double_compares
{
	/** Function that compares for strong equality between two doubles with default tolerance of 100*(Machine precision) */
	bool d_equal_str(real i, real j);
	/** Function that compares for weak equality between two doubles with default tolerance of 100*(Machine precision) */	
	bool d_equal_wk(real i, real j);
	/** Function that compares for strong equality between two doubles with custom tolerance, tol */	
	bool d_equal_str(real i, real j, real tol);
	/** Function that compares for weak equality between two doubles with custom tolerance, tol */	
	bool d_equal_wk(real i, real j, real tol);
	bool d_equal_abs_str(real i, real j, real tol);
	bool d_equal_abs_str(real i, real j);
	bool d_equal_str_gen(real i, real j);
}

#endif
   
