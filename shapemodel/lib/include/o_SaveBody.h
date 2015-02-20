#ifndef	_o_saveworld_incl_
#  define	_o_saveworld_incl_
//
//      o_SaveBody.h  - Header for saveworld class
//
//      Oliver Grau, 16-JUN-1993
//

#include	"fileio.h"
#include	"o_World.h"


/*!
  \class o_SaveBody o_SaveBody.h
  \brief save o_Body to file 
  */

class	o_SaveBody {
	protected:
		o_Body *body;
		int	externfp;
		write_file *of;
		int	flg_pseudo;
		char	*comment;	// optional comment for file headers
	public:
		//! creates file object.
		o_SaveBody( ) ;

		virtual ~o_SaveBody( ) ;

		//! starts writing body and all sub-surfaces to file
		virtual int	Save(o_Body &body);

		//! open file fn for saving. Returns 0 on error
		virtual int	Open(const char *fn);

		//! close assigned file
		virtual void	Close();

		//! set optional comment for file header
		virtual void	SetComment(char *s);

		//! returns value != 0 is pseudo coloring is enabled
		virtual int     GetPseudoColor() { return flg_pseudo;};

		//! set coloring mode to pseudo
		virtual void     SetPseudoColor(int pf) { flg_pseudo = pf;};

		// old constructor
		// old interface
		o_SaveBody( o_Body &body, write_file &of );
		virtual void Set_Body( o_Body &body );
		virtual void Set_File( write_file &of );

		// auxilliary methods
		virtual	void	Save_All();
		virtual	void	Save_Pointlist();
		virtual	void	Save_Edges();
		virtual	void	Save_Polygons();
};

#endif	/* _o_saveworld_incl_ */

