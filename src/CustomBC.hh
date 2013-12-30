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

/** User defined function for boundary condition modifications */
/** Example parabolic profile for Poiseuille flow
for (int j = 1; j <= ny; ++j)
	U_WBCv(j) = 2*U_WBC*(1-(grid_.YC()(1,j)-0.5)*(grid_.YC()(1,j)-0.5)/(0.5*0.5)); */

switch(specialProblems){
	case 1:
		for (int j = 1; j <= ny; ++j)
			if(grid_.YC()(1,j) >= conf_.getRealParameter("RectangleY2"))
				U_WBCv(j) = 1;
			else
				U_WBCv(j) = 0.;
		break;
	default: break;
}
