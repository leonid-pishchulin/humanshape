#ifndef	_incl_texturebinding
#define	_incl_texturebinding
//
//      o_TextureBinding.h  - Header for oGeM texturebinding
//
//      Oliver Grau, 15-JUL-1993
//

#include	"o_Property.h"
#include <o_Vector.h>
#include <vector>
using std::vector;


/*! 
  \class o_TextureBinding o_TextureBinding.h
  \brief class for texturebinding
  */

class	o_TextureBinding : public o_Property {
	public:
		o_TextureBinding() ;
		~o_TextureBinding();
                virtual const  char  *GetClassName() const;

		//! adds a coordinate at the end of the coordinate list
		void	AddCoordinate( o_Vector2D &coord );

		//! sets a coordinate of entry no
		void	SetCoordinate( int no, o_Vector2D &coord );

		/*! returns Coordinate number no, or NULL if no is not a valid
		index range */
		o_Vector2D	*GetCoordinate( int no );

		//! copy coordinates of texture binding
		void Copy(o_TextureBinding &bind);

		void    Copy(o_Property &org);
		virtual o_Property *MakeInstance();

		//! returns the number of coordinates
		int NoOfCoord();


	private:
    std::vector<o_Vector2D*>	coordlist;
};

#endif
