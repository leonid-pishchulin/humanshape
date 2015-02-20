//
//	Oliver Grau, May.1994
//	Oliver Grau, ARP-1996
//

/* $Id: o_Triangle.h 433 2007-09-27 06:10:41Z thormae $ */


#ifndef	_o_Triangle
#  define	_o_Triangle


#include	"o_PlanarFace.h"
//#include	"o_Picture.h"


//forward declaration
class o_Camera;

/*! 
  \class o_Triangle o_Triangle.h
  \brief 3-D triangle class

  triangles can build up using three points as constructor arguments or
  by adding them with AddPoint(). Note that a triangle must be unique.
  That means a triangle can only inserted one time in a body. The
  insertion of a triangle into the parent surface is done during the
  instantiation. On the deletion the triangle is removed from the
  parent surface. So do not remove it explicitly! Further on the 
  points included are stripped from the body point list and deleted
  if they are not referenced by other triangles (reference count is 
  equal 2)
 */

class	o_Triangle : public o_PlanarFace, public o_MetaSurface {
  public:
    const   char    *GetClassName() const;

		/*! appends a point to the triangle.
    The order of the points must be counterclockwise when
    locking on the surface.
    Takes care that the point will be
    inserted in the body point list if necassary */
    void	AddPoint( o_Point &p );

		/*! copies triangle contents. The method tries to use the
    o_Point object with the same index in the associated body
    as in the source triangle. If the points with the same indices
    have diffent coordinates a copy of the point object is
    made !!! */
    void	Copy(o_Triangle &org);

		//! create triangle as sub-surface of surf
    o_Triangle( o_Surface & );

		/*! create Triangle. The order of the points must be
    counterclockwise when locking on the surface */
    o_Triangle( o_Surface &s,o_Point &p1, o_Point &p2, o_Point &p3);
	o_Triangle( o_Surface &s,o_Point &p1, o_Point &p2, o_Point &p3, const std::string &l);

		/*! Deletes triangle. DON'T USE THIS DIRECTLY !! Instead use
    Method "SubTriangle(o_Triangle &)" of o_Surface, which handles
    proper registration of the triangle within the surface 
    structure */
    ~o_Triangle( );
		
		/*! gets texture for triangle, the texture map is a scrap of
    the image of the camera o_Camera . If a map and binding
    exitist before these are stripped and deleted if the reference
    count is zero. */
    void  Get_Texture(o_Camera &cam);

    //
		//      virtual methods for parent-class o_LocalCoordinateSystem
    //
    virtual const o_Vector3D  &Origin();
    virtual const o_Vector3D  &Axis(int n);

		//! returns area of the defined surface
    double  GetArea();

		/*! add binding of projected vertices in image plane
    of the camera */
    void SetBinding(o_Camera &cam);

		//! moves a triangle to the surface 
    void Move( o_Surface &surf );

		//! returns the number of same vertices.
    int NumberOfSameVertices(o_Triangle &poly );

		/*! writes triangles into the list of neighbors, which have at 
    least 'num' common vertices */
    int Neighbors(std::vector<o_Triangle*> &neighborlist, int number);

		/*! appends points to the list of points, which belongs to both
    triangles */
    void SetSamePoints(std::vector<o_Point*> &point_list, o_Triangle &tri);
		//! returns the index of point
    int PointIndex(o_Point &point);

    /*! returns a reference to the triangle which has points p1 and p2 in common with this one, null-pointer if no neighbour was found */
    o_Triangle *GetNeighbour(o_Point &p1, o_Point &p2);

		/*! returns the area of projected triangle in image plane of 
    the camera  o_Camera (mm x mm) */
    double AreaInCameraPlane(o_Camera &camera);  

    double Resolution(o_Camera &camera);

		//! connect the object to a new surface or body
    void    ConnectRootNode(o_Surface &s);

		/*! Returns the underlying plane or NULL, if plane could not be
    computed. */
    const o_Plane *GetPlane();

		// -- old calls only inserted for compatiblity
		// do not use the next call, because it is obsolete
    void	AppendPoint( o_Point &p );
	void setLabel(const std::string &l) {label = l;}
	const char *getLabel(){ return label.c_str(); }

  protected:
	  std::string label;

};


#endif
