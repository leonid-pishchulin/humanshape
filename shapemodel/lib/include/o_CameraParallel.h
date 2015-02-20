
#ifndef	_o_CameraParallel
# define	_o_CameraParallel

//
//	o_CameraParallel.h	- Camera class definitions
//
//	Ralf Toenjes, Jun.1995
//	


#include "o_Picture.h"
#include "o_World.h"
#include "o_Camera.h"



//
const	double	f_infinite=(+1e38);

/*!
  \class o_CameraParallel o_CameraParallel.h
  \brief class for real world camera with parallel projection
  */

class	o_CameraParallel : public o_Camera {
 public:
  //!  creates camera
  o_CameraParallel( );

  ~o_CameraParallel();
  virtual const   char    *GetClassName() const;

  /*! calculates the position of the projection of point in picture coordinates <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> of the camera . The  results are stored to x and y. */
  virtual void	Project_PCS( const o_Vector3D &, double &, double & ) const ;

  /*! computes additionally the distance of point to the 
    camera center in direction of optical axis, stored in z. */
  virtual void	Project_PCS( const o_Vector3D &, double &, double &, double & ) const ;

  /*! returns the line of sight as a unit vector. x and y are given in <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> */
  const   o_Vector3D &LineOfSight(const double x, const double y);
 protected:

 private:

};

#endif


