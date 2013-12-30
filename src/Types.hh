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

#ifndef TYPES_HH
#define TYPES_HH


// This typedef makes it possible to switch between float and double accuracy
// please do not use "float" or "double" directly in your code, but use real instead
typedef double real;



// Enumeration of boundary conditions
typedef enum { NOSLIP, SLIP, OUTFLOW, PERIODIC } BCTYPE;


                          /* Macros for the integer array FLAG      */
#define C_B      0x0000   /* interior obstacle cells                */
#define B_N      0x0001   /* obstacle cells adjacent to fluid cells */
#define B_S      0x0002   /* in the respective direction            */
#define B_W      0x0004   
#define B_O      0x0008 
#define B_NW     0x0005    
#define B_SW     0x0006      
#define B_NO     0x0009    
#define B_SO     0x000a

#define C_F      0x0010   /* fluid cell */

/* Macros for POISSON, denoting whether there is an obstacle cell */
/* adjacent to some direction                                     */
/*#define eps_E (mesh(i+1,j) & C_F)
#define eps_W (mesh(i-1,j) & C_F)
#define eps_N (mesh(i,j+1) & C_F)
#define eps_S (mesh(i,j-1) & C_F) */

#endif //TYPES_HH
