#ifndef	_o_camera
# define	_o_camera

//
//	o_VirtualCamera.h	- Camera class definitions
//
//	Oliver Grau, Jan./Feb.1993
//	Hellward broszio Feb. 1994
//		-added VisiblePixel() in class o_LabelCamera 


#include "o_World.h"
#include "o_Picture.h"
#include "o_MetaCamera.h"

class	o_Surface;

/*!
  \class o_VirtualCamera o_VirtualCamera.h
  \brief class for virtual cahv camera - meta class is not allowed to be instantiated
  */

//
//	Meta class do not instantiate
//

class	o_VirtualCamera : public o_MetaCamera {
	public:
                //! Clear Render-Status (image)
		virtual	void	Clear( );

		/*! Render a world-hierachie with all bodies/surfaces
		performs a Clear() before starting the render process */
		virtual	void	Render( o_World & );
		
		//! Render a body-hierachie with all sub-/surfaces
		virtual	void	Render_B( o_World&, o_Body &);

		//! Render a surface with all sub-/surfaces
		virtual	void	Render_S( o_World&, o_Surface & );

		//! Render a polygon
		virtual	void	P_render( const o_World& wp, o_Triangle &pp );

		//void	Describe( write_file &wfp, int bin_flag=0 );
		//virtual int	Readparameter( read_file &rfp );
		virtual const   char    *GetClassName() const;
	protected:
		o_VirtualCamera( );
		~o_VirtualCamera( );
};

/*!
  \class o_LabelCamera o_VirtualCamera.h
  \brief generates label image of scene-objects
  */

class	o_LabelCamera : public o_VirtualCamera {
	public:
                //! creates camera
		o_LabelCamera() ;

		~o_LabelCamera();
		void	P_render( const o_World&, o_Triangle & );

		//! returns label picture
		const o_Picture<long>	&GetPicture() const;

		//! returns z-buffer image
		const o_Picture <float>	&GetZBuffer()  const;

		virtual const   char    *GetClassName() const;
		virtual	void	Clear( );
		
		/*! returns the number of visible pixel from the camera
		  viewpoint. The return value is normalized to maximal
		  potential visible pixel in the range [0.0 1.0]. */
		double VisiblePixel(o_Triangle &polygon);

		//! returns the number of visible pixel from the camera viewpoint.
		int NoOfVisiblePixel(o_Triangle &polygon);

		//! deletes the z-buffer, creates no new one
                void DeleteZBuffer();

		//! deletes the label picture, creates no new one
                void DeleteLabelPic();
	protected:
		o_Picture<long>	*labpic;
		o_Picture <float>	*zbuffer;
};

#ifdef	notdef
/*!
  \class o_GreyscaleCamera o_VirtualCamera.h
  \brief generates greyscale image of scene
  */

class	o_GreyscaleCamera : public o_VirtualCamera {
	public:
                //! creates camera
		o_GreyscaleCamera();

		~o_GreyscaleCamera();

		void	P_render( const o_World&, o_Triangle & );
		virtual const   char    *GetClassName() const;

		//! returns render result
		const UByte_image	&GetPicture() const;

		//! returns z-buffer image
		const Picture <float>	&GetZBuffer()  const;

		virtual	void	Clear( );
	protected:
		UByte_image	*pic;
		Picture <float>	*zbuffer;
};
#endif


#endif
