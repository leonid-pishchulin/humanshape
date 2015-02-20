#ifndef _Clip3d_
#define _Clip3d_
#define VISIBLE_DEPTH 1e-30
#include "o_TriangleClip2D.h"
#include "o_Triangle.h"

//
//	o_TriangleClip3D.h	- class for walking over a polygon in image plane
//
//	
//

/*!
  \class o_TriangleClip3D o_TriangleClip3D.h
  \brief Performs clipping of 
  */

class o_TriangleClip3D : public o_TriangleClip2D
{
 public:
           ~o_TriangleClip3D();

	   //! creates triangle from given information
           o_TriangleClip3D( o_MetaCamera &,  o_Triangle &);

	   //! returns value unequal 0 if polygon is in front of the camera
           int In_front_of_camera();

	   //! returns value unequal 0 if polygon is definitively not visible
           int Is_not_visible();

	   /*! returns Z position of vertex no ( no < range ) in image
	     plane (camera projection) */
           double VertexZ(int no, int range=3);

	   /*! returns value unequal 0 if polygon normal point not in
	     the camera direction */
           int Surface_to_camera(){ return surf_to_cam; }

	   /*! projection of triangle vertices with given camera,
                including clipping of triangles. (Triangles lieing 
		partiately behind the camera center are cut yielding 
		new triangles or quadrangles lieing in front of the camera.) */
           void projection3d(o_MetaCamera &cam, o_Triangle &pp);

	   /*! Clipping of triangles may result in quadrangles. This function
		returns the location of vertex 4.*/
           const int GetVertex4Location() const { return Vertex4Location; }
 private:
           double z[3];
           int Vertex4Location; 
           int surf_to_cam;
           
           void clipping3d(o_MetaCamera &cam, double *pcH, double *pcV, int VisibleVertices); 
           void correction3dOf2InvisibleVertices(o_MetaCamera &cam, double *pcH, double *pcV, double *x_tmp, double *y_tmp, int *LocationCode);
           void correction3dOf1InvisibleVertex(o_MetaCamera &cam, double *pcH, double *pcV, double *x_tmp, double *y_tmp, int *LocationCode);
           double ViewPlaneOffset_x(int width);
           double ViewPlaneOffset_y(int height);
         
};
          
#endif
