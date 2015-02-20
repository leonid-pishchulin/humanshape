#ifndef _Clip2d_
#define _Clip2d_

//
//	o_TriangleClip2D.h	- class for walking over a triangle in image plane
//
//	
//

#include        "o_VirtualCamera.h"

/*!
  \class o_TriangleClip2D o_TriangleClip2D.h
  \brief Performs stepping over all points in image plane belonging to a triangle
  */

class o_TriangleClip2D 
{
 public:
           o_TriangleClip2D();
           ~o_TriangleClip2D();

	   //! calculation of shape equation parameters
           void clipping2d( int dx, int dy, double * x, double * y);

	   //! init iterator
           void Init();

	   /*! returns value unequal 0 if there is an other row in the image
		in which the triangel lies. In this case y contains the image
		row, col the start column and endcol the endcolumn plus one. */
           int Nextline( int& y, int& col, int& endcol);

	   /*! returns value 0 if all 3 vertices of polygon are not
	     in the image plane */
           int In_camera_frame();

	   //! returns value unequal 0 if polygon is definitively not visible
           int Is_not_visible();

	   /*! returns position of vertex no ( no < range ) in image
	     plane (camera projection) */
           double VertexX( int no, int range=3 );
           double VertexY( int no, int range=3 );

	   //! returns number of polygon vertices
           int Vertices() { return 3; } 
 
 protected:
           double x[4],y[4];
 
 private:
           double as[3], bs[3], ds[3];
           double ae[3], be[3], de[3];
  	   int image_dx, image_dy;
           double startline,endline;
           int startedges,endedges;
           int scanline, scanend;
           double startcol, endcol;
            
};

#endif
