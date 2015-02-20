#ifndef	_o_lightsource_incl
#  define	_o_lightsource_incl
//
//	o_LightSource.h	- Header for oGeM lightsource class
//
//	Oliver Grau, Jan.1993
//

#include	"oGeM.h"
#include	"o_Color.h"


/*!
  \class o_LightSource o_LightSource.h
  \brief directional light source class

  abstract class do not instantiate
  */

class	o_LightSource : public o_MetaClass 	{
	public:
		virtual const   char    *GetClassName() const;

		//! set light source color (default is rgb(1,1,1))
		void SetColor(const o_Color_f c);
		//! return light source color
		o_Color_f GetColor();

		//! get light source intensity 
		double GetIntensity();
		//! set light source intensity (default: 1)
		void SetIntensity(double intens);

		//! return 0 if lightsource is only directional
		virtual int	HasPosition()=0;	
		//! return 0 if lightsource has no direction (point source)
		virtual int	HasDirection()=0;	


		//! set position of light source. Only applicable is light source has a position (check with HasPosition())
		virtual void SetPosition(const o_Vector3D pos);
		//! get position of light source. Only applicable is light source has a position (check with HasPosition())
		virtual o_Vector3D GetPosition();

		//! set orientation of light source. Only applicable is light source has a orientation (check with HasDirection())
		virtual void SetDirection(const o_Vector3D dir);
		//! get orientation of light source. Only applicable is light source has an orientation (check with HasDirection())
		virtual o_Vector3D GetDirection();

	protected:
		o_Color_f	color;
		double		intensity;
		friend class o_World;

		o_Vector3D position;
		o_Vector3D direction;


		//! create light source with name nam
		o_LightSource( );

		virtual ~o_LightSource();
};



/*!
  \class o_DirectionalLight o_LightSource.h
  \brief directional light source class
  */

class	o_DirectionalLight : public o_LightSource 	{
	public:
                //! create directional light source
		o_DirectionalLight();

		virtual ~o_DirectionalLight();

		virtual const   char    *GetClassName() const;
		virtual int	HasPosition() {return 0;}
		virtual int	HasDirection(){return 1;}
};


/*!
  \class o_PointLight o_LightSource.h
  \brief point light source class
  */

class	o_PointLight : public o_LightSource 	{
	public:
                //! create directional light source
		o_PointLight();

		virtual ~o_PointLight();

		virtual const   char    *GetClassName() const;
		virtual int	HasPosition() {return 1;}
		virtual int	HasDirection(){return 1;}

};


#endif
