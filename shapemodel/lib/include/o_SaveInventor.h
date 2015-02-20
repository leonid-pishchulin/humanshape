#ifndef	_o_saveinventor_incl_
#  define	_o_saveinventor_incl_
//
//      o_SaveInventor.h  - Header for o_SaveInventor class
//
//      Oliver Grau, 18-JAN-1995
//

#include	"o_SaveInventorBase.h"


/*!
  \class o_SaveInventor o_SaveInventor.h
  \brief save Body as Inventor file
  */

class	o_SaveInventor : public o_SaveInventorBase {
	public:
                //! creates inventor-file object.
		o_SaveInventor( );

		//!  Saves the world w and cameralist cl to inventor-file
    int SaveWorld(std::vector<o_MetaCamera*> &cl, o_World &w);

		// old constructor
		o_SaveInventor( o_Body &body, write_file &of );
	private:
		void SaveCamera( o_MetaCamera *b );
		virtual void    Save_Body( o_Body *b );
		virtual void    Save_Surface( o_Surface *s );
};

#endif

