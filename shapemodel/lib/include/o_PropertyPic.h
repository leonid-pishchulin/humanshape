#ifndef	_o_PropertyPic_incl
#  define	_o_PropertyPic_incl

#include	"o_Property.h"
#include	"o_Picture.h"

/*!
  \class o_PropertyPic_byte o_PropertyPic.h
  \brief Picture Property class
  */

class	o_PropertyPic_byte : public o_Property {
	public:
                //! creates o_PropertyPic_byte object
		o_PropertyPic_byte() {};

		~o_PropertyPic_byte() {};
		o_Picture<unsigned char>	pic;
		virtual	const char	*GetClassName() const;
		virtual o_Property *MakeInstance();
		virtual void    Copy(o_Property &org);
};


#endif
