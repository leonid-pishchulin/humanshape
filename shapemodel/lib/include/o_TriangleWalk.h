#include	"o_Triangle.h"
#include	"o_TriangleClip2D.h"
#include        <math.h>
#include        "o_World.h"
#include        "o_Surface.h"


#ifndef o_TriangleWalk_h
#define o_TriangleWalk_h

/*!
  \class o_TriangleWalk o_TriangleWalk.h
  \brief Steps over all points of a triangle in 3D space.

  Steps over all points of a triangle in 3D space. The resolution can be set arbitrarily. The triangle is stepped in lines parallel to the base line (0/1).
  */

class o_TriangleWalk : public o_TriangleClip2D {
  o_Triangle * tr;
  double l1;
  int res_01;
  int res_h;
  double w;
  double h;
  int col, endcol, scanline;

public:
  /*! The Object is created and initialized. t is the triangle to walk on,
    res10 the number of steps in base direction of the triangle, resh the
    number of steps in the direction of the height of the triangle. */
  o_TriangleWalk(o_Triangle & t, int res10, int resh);

  ~o_TriangleWalk ();
//  int Init();

  /*! Tests if there any points left on the triangle surface to step to. 
    If there arent any left returns 0. */
  int NotFinished();

  /*! Increments internal column counter and tests whether to increment row
    counter. Afterwards current spatial coordinates on triangle surface are
    returned to point3d, while triangle related 2D coordinates are returned
    to point2d. If all points on the triangle surface have been visited,
    the method returns 0. */
  int NextPoint(o_Vector3D & point3d, o_Vector2D & point2d);


};





#endif

