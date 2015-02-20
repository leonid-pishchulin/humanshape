#ifndef	_o_SaveTCL_incl_
#  define	_o_SaveTCL_incl_
//
//      o_SaveTCL.h  - Header for SaveTCL class
//
//      Oliver Grau, 16-JUN-1993
//      Gerwin Matysiak, 23-Jan-1995
//
#include	"fileio.h"
#include	"o_World.h"
#include	"o_SaveBody.h"


/*!
  \class o_SaveTCL o_SaveTCL.h
  \brief save Body with all surfaces to TCL-file
  */

class	o_SaveTCL : public o_SaveBody {
        private:
                int aux_swfno; //for numbering the surfaces
	public:
		//! creates TCL-file object.
		o_SaveTCL( );

		o_SaveTCL( o_Body &body, write_file &of );
		virtual	void	Save_All(); //for all surfaces
		virtual	void	Save_Pointlist(o_Body *b);
		virtual	void	Save_Body(o_Body *b);
   
//              virtual void    Save_Surface();

//              virtual void    Save_Surface(o_Body *b);
                virtual void    Save_Surface(o_Surface *s);

	protected:
		int	ext;	// flag for the tcl format
};

#endif	/* _o_SaveTCL_incl_ */

