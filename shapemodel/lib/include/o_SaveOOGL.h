#ifndef	_o_saveoogl_incl_
#  define	_o_saveoogl_incl_
//
//      o_SaveOOGL.h  - Header for o_SaveOOGL class
//
//      Oliver Grau, 18-JAN-1995
//

#include	"fileio.h"
#include	"o_World.h"
#include	"o_SaveBody.h"


/*!
  \class o_SaveOOGL o_SaveOOGL.h
  \brief save Body as OOGL (Geomview) file
  */

class	o_SaveOOGL : public o_SaveBody {
	public:
                /*! creates OOGL-File, \sa <A HREF="http://www.tnt.uni-hannover.de/soft/compgraph/fileformats/">fileformats</A> */
		o_SaveOOGL( );

		// old constructor
		o_SaveOOGL( o_Body &body, write_file &of );
		virtual	void	Save_All();
		virtual	void	Save_Pointlist();
		virtual	void	Save_Polygons();
	private:
		void    dopolygon( write_file *of, o_Surface *s, int cnt );
};

#endif

