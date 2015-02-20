#ifndef	_incl_o_material
#define	_incl_o_material
//
//      o_Material.h  - Header for oGeM material
//
//      Oliver Grau, 27-OCT-1993
//

#include	"o_Property.h"
#include	"o_Color.h"

/*! 
  \class o_Material o_Material.h
  \brief meta class for material

  For a description look into \link mat_para.html Material\endlink-intro.
  */

class	o_Material : public o_Property     {
	public:
		o_Material() ;
		~o_Material();
		virtual const   char    *GetClassName() const;

		/*! returns an infostring describing the object-values.
		  Possible modes are 0: long infostring, 1: short infostring */
		const char *GetInfostring(int mode);

		o_Color_f	ambient;
		o_Color_f	diffuse;
		o_Color_f	specular;
		o_Color_f	emissive;
		float	shininess;
		float	transparency;
		o_Material	&Copy( const o_Material &o ) {
			ambient = o.ambient;
			diffuse = o.diffuse;
			specular = o.specular;
			emissive = o.emissive;
			shininess = shininess;
			transparency = transparency;
			return *this;
		};
		o_Property *MakeInstance();

		//! copy contents of mat
		void    Copy(o_Property &org );

};


extern	o_Material	material_default;

#endif
