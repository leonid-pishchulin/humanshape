#ifndef	_o_savestp_incl_
#  define	_o_savestp_incl_
//
//      o_SaveSTP.h  - Header for savestp class
//
//      Oliver Grau, 16-JUN-1993
//

#include	"fileio.h"
#include	"o_World.h"
#include	"o_SaveBody.h"


/*!
  \class o_SaveSTP o_SaveSTP.h
  \brief save Body to stp-file (TRI)
  */

class	o_SaveSTP : public o_SaveBody {
	public:
                /*! creates object for saving body object to stp file.
		(\sa <A HREF="http://www.tnt.uni-hannover.de/project/eu/panorama/fileformats/stp.html">STP-Fileformat</A> ) */
		o_SaveSTP();
		int Save(o_Body &body);

		// old constructor call
		o_SaveSTP( o_Body &body, write_file &of, int extendedformat=1 );

		virtual	void	Save_All();
		virtual	void	Save_Pointlist();
		virtual	void	Save_Polygons();
	protected:
		int	ext;	// flag for the extended stp format
};

#endif	/* _o_savestp_incl_ */

