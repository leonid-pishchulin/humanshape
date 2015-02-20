#ifndef	_o_MeasureCamera_incl
#  define	_o_MeasureCamera_incl
//
//	o_MeasureCamera.h	- Header for the Measurement camera class
//
//	Georg Brink, Mar.1994
//

#include        "smatrix.h"
#include	"o_Camera.h"
#include        "o_HMatrix.h"
#include        "stddef.h"
#include        "o_BilinInterpolation.h"


/*!
  \class o_MeasureCamera o_MeasureCamera.h
  \brief Measurement Camera class
  */

class	o_Vector3D;
class   leda::matrix;
 
class	o_MeasureCamera : public o_Camera {
	friend class o_StereoSensor;
        public:
                //! creates MeasureCamera 
  	       	o_MeasureCamera() ;

	        ~o_MeasureCamera();
		
		//! set parameter like in cam
                virtual void Set_Camera(o_MetaCamera &cam);

		//! set measure camera parameter
                virtual void Set_Camera(o_Vector3D rvec, o_Vector3D tvec, leda::matrix cmat, double f=1,
                                               double xh=0, double yh=0, int wid=CAM_DEFAULT_WID,
                                               int hei=CAM_DEFAULT_HEI, double rd=0, double sx =1,
                                               double sy=1);

		//! set camera covariance matrix (6*6)
		void SetCovariance( leda::matrix &cmat);

		//! get camera covariance matrix (6*6)
		leda::matrix &GetCovariance();

		//! returns depthmap
                o_Picture<float> &GetDepthmap();

		//! set depthmap
                void         SetDepthmap(o_Picture<float> &dm);

		//! returns segmentationmap
                o_Picture<short>  &GetSegmentmap();

		//! set segmentationmap
                void         SetSegmentmap(o_Picture<short> &sm);

		//! returns certaintymap (certainty of depth)
                o_Picture<float> &GetCertaintymap();

		//! set certaintymap (certainty of depth)
                void SetCertaintymap(o_Picture<float> &cm);

		/*! returns 3D-position (<A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#wcs">WCS</A> ) of an image point x/y ( <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> ). The depth of the image point is taken from the depthmap. */
                virtual const o_Vector3D &Position(double bx, double by);

		/*! returns 3D-position (<A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#wcs">WCS</A> ) of an image point x/y ( <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> ). The depth of the image point is given by tau. */
                virtual const o_Vector3D &Position(double bx, double by, double tau);

		/*! returns certainty of 3D-position of image vector icev
		  ( <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> ). Depth is given by tau, depth variance is given by tau_var. If only ivec is given, tau and tau_var are taken from depth and certainty maps of MeasureCamera. */
                const leda::matrix &PositionCertainty(o_Vector2D ivec, 
						double tau = 0.0,
						double tau_cam_var = 0.0);

		/*! returns the projection of the covariance matrix cmat
		  in the image plane. The space point belonging to cmat is given by svec (<A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#wcs">WCS</A> ). */
                const leda::matrix &ImagePointCertainty(o_Vector3D &svec, leda::matrix &cmat);
		//! clear depthmap, certaintymap and segmentationmap
                virtual void Clear();              

               	virtual const   char    *GetClassName() const;
	
         protected:
	        leda::matrix       covmat;
                o_Picture<float> *depthmap;
                o_Picture<float> *certaintymap;
                o_Picture<short>  *segmentmap;

 };


#endif






