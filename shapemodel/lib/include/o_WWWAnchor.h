#ifndef	_incl_o_WWWAnchor
#define	_incl_o_WWWAnchor
//
//      o_WWWAnchor.h  - Header for oGeM WWWAnchor
//
//      Oliver Grau, 27-OCT-1993
//

#include	"o_Property.h"
#include	"o_Color.h"


/*!
  \class o_WWWAnchor o_WWWAnchor.h
  \brief class for WWWAnchor
  */

class	o_WWWAnchor : public o_Property     {
	public:
		o_WWWAnchor() ;
		~o_WWWAnchor();
		virtual const   char    *GetClassName() const;
		o_Property *MakeInstance();
		void Copy(o_Property &org);

		//! set value (name) of the anchor
		void	SetAnchor(const char *an) ;

		//! returns value (name) of the anchor
		const char	*GetAnchor() const;

		//! copy contents of mat
		o_WWWAnchor	&Copy( const o_WWWAnchor &o ) ;
	private:
		char	*anchor;
		int	map;
};



#endif
