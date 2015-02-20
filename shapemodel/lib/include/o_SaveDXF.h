#ifndef	_o_saveDXF_incl_
#  define	_o_saveDXF_incl_
//
//      o_SaveDXF.h  - Header for saveDXF class
//
//      Oliver Grau, 06-ARP-1994
//

#include	"fileio.h"
#include	"o_World.h"
#include	"o_SaveBody.h"


/*!
  \class o_SaveDXF o_SaveDXF.h
  \brief save Body to DXF-file (Autocad)
  */

class	o_SaveDXF : public o_SaveBody {
	public:
                //! creates DXF-file object
		o_SaveDXF( );

		int Save( o_Body &body );

		// old constructor
		o_SaveDXF( o_Body &body, write_file &of );

		virtual	void	Save_All();
		virtual	void	Save_Pointlist();
		virtual	void	Save_Polygons();
	protected:
};

#endif

