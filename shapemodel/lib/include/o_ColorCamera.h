#ifndef	_o_colorcamera
# define	_o_colorcamera

//
//	o_ColorCamera.h	- Camera class definitions
//
//	Oliver Grau, Jan./Feb.1993
//      Ralf Toenjes 11.10.94 Background, Min-/Magfilter


#include "o_VirtualCamera.h"

#define	ZBufferType	double

//
// Minification and magnification Filter (for colorcamera and texturemap)
//
#define NO_SPECIFICATION -1
#define MIN_POINT         0    
#define MAG_POINT         1
#define MIN_BILINEAR      2
#define MAG_BILINEAR      3
#define MIN_AVERAGE       4

/*!
  \class o_ColorCamera o_ColorCamera.h
  \brief generates color image of scene
  */


class	o_ColorCamera : public o_VirtualCamera {
	public:
#define DEFAULT_BACKGROUND           50
#define DEFAULT_MIN_FILTER           MIN_POINT
#define DEFAULT_MAG_FILTER           MAG_POINT
#define DEFAULT_3D_PROJECTION_DEPTH  -1

                //! creates color camera
		o_ColorCamera();

		~o_ColorCamera();
		/*! Render a world-hierachie with all bodies/surfaces
		    performs a Clear() before starting the render process */
		virtual void    Render( o_World & );
						 

		void	P_render( const o_World&, o_Triangle & );
		virtual const   char    *GetClassName() const;

		/*! returns render result. The Renderer supports only full
		  oblique and transparent pixels. All other alpha values
		  are thresholded with 0.5 to switch between oblique and
		  transparent. */
		const o_CPicture<unsigned char>	&GetPicture() const;

		//! returns z-buffer image (ZBufferType == double)
		const o_Picture <ZBufferType>	&GetZBuffer()  const;

		virtual	void	Clear( );

		/*! copy picture into camera and render on top of it. 
		    Must be set before the Render() call */
		void SetBackgroundImage(o_CPicture<unsigned char> &pic);

		/*! sets background value for band band_no of colorimage.
                The setting is activated by calling o_ColorCamera:Clear() ,
		e.g. when reading the camera file.
		The DEFAULT_BACKGROUND is 50. */
		void SetBackground(UByte_ value, int band_no=0);

		//! returns background value
		UByte_ GetBackground(int band_no=0);

		/*! Sets switch-depth. Texture of polygons closer than 
		  switch-depth is projected using 3D transforms instead of
		  2D affine transform. 3Ddepth = -1 forces 2D affine 
		  transformations. Attention: The texture may still appear
		  distorted, if the texture mapping did not use
		  3D projections for texturing. */
                void Set3DProjectionDepth( double depth_3D );

		/*! returns switch-depth. Texture of polygons closer than
		  switch-depth is projected using 3D transforms instead
		  of 2D affine transform. */
                double Get3DProjectionDepth();

		/*! Sets minification filter method. \sa \link texture.html
		  texture \endlink. Select one of NO_SPECIFICATION,
		  MIN_POINT, MIN_BILINEAR, MIN_AVERAGE. */
                void SetMinFilter( int filter_method );

		/*! returns minification filter method \sa \link 
		  texture.html texture \endlink */
                int  GetMinFilter(); 

		/*! Sets magnification filter method \sa \link texture.html
		  texture \endlink . Select one of NO_SPECIFICATION,
		  MAG_POINT, MAG_BILINEAR. */

                void SetMagFilter( int filter_method );

		/*! returns magnification filter method \sa \link
		  texture.html texture \endlink  */
                int  GetMagFilter(); 
	protected:
		o_CPicture<unsigned char>	*pic;
		o_Picture <ZBufferType>	*zbuffer;
        private:
		UByte_  background[3];
                int minfilter;
                int magfilter;
                double projection_depth_3D;
		int	bg_img_flag;	// set if background image was set

};


#endif
