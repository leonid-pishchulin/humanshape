#ifndef	_o_o_PlanarFace
#  define	_o_o_PlanarFace
//
//	o_o_PlanarFace.h	- Header for oGeM face classes (triangle,..)
//
//	Oliver Grau, MAY-1994
//


#include	"o_MetaSurface.h"
#include	"o_Plane.h"


class	o_Triangle;
class o_TextureMap;
class o_Camera;
class o_Point;
class o_Surface;
class o_Body;
class o_World;

/*!
 \class o_PlanarFace o_PlanarFace.h
 \brief abstract class
 
 abstract class - do not instantiate
 */

class	o_PlanarFace {
  public:
     /*! calls AppendPoint() (must be implemented in the derived calls )
     and takes care that the point will be
    inserted in the body point list if necessary */
    void	AddPoint( o_Point &p, o_Body &body);

		/*! appends a point to the polygon.
    The order of the points must be counterclockwise when
    locking on the surface. The Point p must also be inserted in
    the body point list */
    virtual void	AppendPoint( o_Point &p ) = 0;
		
		//! returns the points number of the polygon
    int	NumberOfPoints() const ;
    
    //! returns a reference to point n of the polygon
    o_Point &Point( int n );
    
		//! returns a pointer to point n of the polygon
    o_Point *GetPoint( int n );
    
    /*! return index [0..NumberOfPoints()] if p is a point
    point or -1 otherwise */
    int PointIndex(  const o_Point &p ) ;
    
    virtual const   char    *GetClassName() const;
    
    virtual const o_Plane *GetPlane() = 0;
  protected:
    o_PlanarFace( );
    virtual ~o_PlanarFace( );

    o_Vector3D        *origin;        // must be up to date
    o_Vector3D        *orient_vec[3]; // must be up to date
    virtual int    UpdateGeometry();
    o_Plane *GetHiddenPlane() const;
    int SetHiddenPlane ();

    //!  inserts point item before the list member succ
    void InsertBoundingPoint ( o_Point & itm, o_Point & succ );
    //! adds point to the end of the point list
    void AddBoundingPoint(o_Point & itm);
    //!  substract point itm from the point list
    bool SubBoundingPoint( o_Point & itm );
    
  private:
    o_Plane *ppl;		//! hidden plane
    std::vector<o_Point*>	bounding_point_list;	//! contains the bounding points
};


#endif
