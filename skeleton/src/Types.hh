#ifndef TYPES_HH
#define TYPES_HH


// This typedef makes it possible to switch between float and double accuracy
// please do not use "float" or "double" directly in your code, but use real instead
typedef double real;



// Enumeration of boundary conditions
typedef enum { NOSLIP, SLIP, OUTFLOW, PERIODIC } BCTYPE;

typedef enum {NORTH, EAST, SOUTH, WEST, CENTER} Direction;


                          /* Macros for the integer array FLAG      */
#define C_B      0x0000   /* interior obstacle cells                */
#define B_N      0x0001   /* obstacle cells adjacent to fluid cells */
#define B_S      0x0002   /* in the respective direction            */
#define B_W      0x0004   
#define B_E      0x0008 
#define B_NW     0x0005    
#define B_SW     0x0006      
#define B_NE     0x0009    
#define B_SE     0x000a

#define C_F      0x0010   /* fluid cell */

/* Macros for POISSON, denoting whether there is an obstacle cell */
/* adjacent to some direction                                     */
#define eps_E !(mesh(i+1,j) < C_F)
#define eps_W !(mesh(i-1,j) < C_F)
#define eps_N !(mesh(i,j+1) < C_F)
#define eps_S !(mesh(i,j-1) < C_F)

#endif //TYPES_HH
