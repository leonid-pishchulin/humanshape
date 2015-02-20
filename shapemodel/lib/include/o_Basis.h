#ifndef	_o_basic_incl
#  define	_o_basic_incl
//
//	o_Basis.h	- Basic objects/definitions
//
//	Oliver Grau, Feb. 1993
//

/*!
  \relates o_World

  Initialize oGeM-Library

  Initializes the oGeM-Library. Builds the property-lookup-table with the
  standard properties of liboGeM. This function MUST be called once at
  the beginning of any oGeM-program.
  */

void	init_oGeM();

/*! 
  \relates o_SurfelCreator

  Initialize oGeM-Measure-Library

  Initializes the oGeM-Measure-Library. Adds entries to the
  property-lookup-table. This function MUST be called once at the
  beginning of any oGeM-program using the measure-library.
  */

void    init_oGeM_measure();

/*!
  \relates o_Texturing

  Initialize oGeM-Texturing-Library

  Initializes the oGeM-Texturing-Library. Adds entries to the
  property-lookup-table. This function MUST be called once at the
  beginning of any oGeM-program using the texturing-library.
  */

void	init_oGeM_texturing();

typedef	unsigned char	UByte_;
typedef	short	Short_;
typedef	long	Long_;
typedef	float	Float_;
typedef	double	Double_;

struct	Pixel_value {
	UByte_	lum;
	UByte_	r,g,b;
};

#ifndef WIN32
#ifdef  CENTERLINE_CLPP
// bool.h is not supported in Centerline 2.1
#  define TRUE 1
#  define FALSE 0
#  define bool	char
#else
#  ifndef MacOS
#    include	<bool.h>
#  endif
#endif
#else
#ifndef TRUE
#  define TRUE 1
#endif

#ifndef FALSE
#  define FALSE 0
#endif
#endif

#endif
