#ifndef	_o_MeasureCameraParallel_incl
#  define	_o_MeasureCameraParallel_incl
//
//	o_MeasureCameraParallel.h	- Header for the Measurement cameraparallel class
//
//	Georg Brink, Mar.1994
//

#include        "smatrix.h"
#include	"o_MeasureCamera.h"
#include        "o_HMatrix.h"
#include        "stddef.h"
#include        "o_BilinInterpolation.h"


/*!
  \class o_MeasureCameraParallel o_MeasureCameraParallel.h
  \brief Measurement Camera class
  */

class	o_Vector3D;
class   matrix;
 
class	o_MeasureCameraParallel : public o_MeasureCamera {
	friend class o_StereoSensor;
        public:
                //! creates MeasureCameraParallel
  	       	o_MeasureCameraParallel() ;

	        ~o_MeasureCameraParallel();

		/*! calculates the position of the projection of point 
		in picture coordinates <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs>PCS</A> of the camera . The 
		results are stored to x and y. */
	       	virtual void Project_PCS( const o_Vector3D &, double &, 
					  double &, double &) const ;

		/*! computes additionally the distance of point to the 
		camera center in direction of optical axis, stored in z. */
	       	virtual void Project_PCS( const o_Vector3D &, double &, 
					  double & ) const ;

		/*! returns the line of sight as a unit vector. x and y
                are given in <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> . */
		const   o_Vector3D &LineOfSight(const double x, const double y);

		/*! returns 3D-position ( <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#wcs">WCS</A> ) of an image point x/y (<A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> ). The depth of the image point is taken from the depthmap. */
		virtual const o_Vector3D &Position(double bx, double by);

		/*! returns 3D-position ( <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#wcs">WCS</A> ) of an image point x/y (<A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> ). The depth of the image point is given by tau. */
                virtual const o_Vector3D &Position(double bx, double by, double tau);

		/*! returns Z-Value for X and Y in <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#wcs">WCS</A> */
                virtual double GetZ(double X, double Y);

		/*! Sets depthmap to  height of camera - heightmap (DEM). 
		  The File may be of other size than specified by the camera.
		  Scaling is assured by Position(x,y). */
		virtual void SetHeightMap(o_Picture<float> &dm);

               	virtual const   char    *GetClassName() const;

         protected:
 };


#endif






