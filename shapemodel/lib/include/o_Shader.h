#ifndef	_incl_o_shader
#define	_incl_o_shader
//
//      o_Shader.h  - Header for oGeM shader
//
//      Oliver Grau, 02-NOV-1993
//

#include	"o_Property.h"
#include	"o_Color.h"

class	o_LightSource;
class	o_TriangleClip3D;

/*!
  \class o_Shader o_Shader.h
  \brief class for shader
  */

class	o_Shader : public o_Property     {
	public:
		o_Shader() ;
		~o_Shader();
		virtual const   char    *GetClassName() const;
		void	Init();
		void       GetShade( o_Color_f &in, int x, int y );

		o_Material *mat;
		o_TriangleClip3D *pp;
    const std::list<o_LightSource*> *lights;
		o_Vector3D	normal;
		o_Property * MakeInstance();
		void Copy(class o_Property &);

	private:
		o_Color_f	diffref;
};



#endif
