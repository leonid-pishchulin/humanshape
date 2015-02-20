
/* SCCS-ID: $Id: oGeM.h 431 2007-09-26 00:26:07Z rhys $ (c) Oliver Grau */

/** \mainpage notitle
 *  \anchor maindoxygenmain
 * 
 * This is the code reference manual for the <a href="./index.html"><b>Main Geometry Modelling Library</b></a> classes of the
 * <a href="../../../overview.html"><b>Digilab Library documentation</b></a>
 */

#ifndef	_o_incl
#  define	_o_incl
//
//	oGeM.h	- Header for oGeM
//
//	Oliver Grau, Jan.1993
//


//extern	"C"	void	exit(int);

#include <fstream>
#include "o_Vector.h"
#include "o_Basis.h"
#include "o_functions.h"

#include "fileio.h"

#ifndef	MAXOBJECT
#  define	MAXOBJECT	4096
// should be a "high" number
#endif

typedef	void	*o_key;

#include "o_MetaClass.h"

#endif
