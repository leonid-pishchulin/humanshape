//
//	o_functions.h	-	Function-Prototypes
//
//	Oliver Grau,23.2.1993
//


class	o_World;
class	read_file;

/*!
  \relates o_World
  
  reads a STP-file and creates the defined structures

  o_World &w is the world object in which the STP-objects are created.
  STP_Readin creates an object if no object with the name
  specified in the stp-file exists. read_file & rf is an open read file
  assigned to the stp-file.

  Returnvalue: 	1 on success or 0 on error
  */

int     STP_Readin( o_World &w, read_file & rf );

/*!
  \relates o_Camera

  solves the z-equation
*/

void	zequation( double *x, double *y, double *z,
	double &a, double &b, double &c);

