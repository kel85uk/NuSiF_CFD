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

/** User defined function for initial condition modifications */

switch(specialProblems){
	case 1: // Just modify the cells that were initialized by default
		for (int j = 1; j <= grid_.Ny(); ++j)
			for (int i = 1; i <= grid_.Nx(); ++i){
				if(mesh_.geom_isFluid(i,j) && grid_.YC()(i,j) <= conf_.getRealParameter("RectangleY2")){
					grid_.U()(i,j) = 0.;
					grid_.V()(i,j) = 0.;
					grid_.p()(i,j) = 1.;
				}
			}
		break;
	default: break;
}
