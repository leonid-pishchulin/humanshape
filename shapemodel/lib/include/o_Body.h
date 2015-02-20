#ifndef	_o_body_incl
#define	_o_body_incl
//
//	o_Body.h	- Header for o_Body class
//
//	Oliver Grau, Nov.1995
//

#include "o_TopologyNode.h"
#include "o_HashTable.h"
#include "o_VertexInfo.h"

#include <list>
using std::list;
#include <vector>
using std::vector;
#include <map>
using std::map;

//
//      Forward declaration of backwardchained classes
//

class   o_Surface;
class   o_NurbsSurface;
class   o_NurbsCurve;
class   o_ObjectTransformation;
class	  o_Polygon;
class	  o_HMatrix;
class   o_Property;
class   o_VirtualCamera;
class   o_Camera;

/*!
  \class o_Body o_Body.h
  \brief Class for body 

  A o_Body object contains a list of (hierachical) o_Surface objects
  and a unique list of o_Point objects. The o_Point objects are used
  as vertices for the surface primitives (e.g. o_Triangle ).
 */

class	o_Body : public o_LocalCoordinateSystem, public o_TopologyNode {
  public:
                //! creates body with name nam 
    o_Body(char *nam = NULL);

		//! creates body with reference to w and name nam
    o_Body(o_World &, char *nam = NULL);

    ~o_Body( );

    virtual	void	Copy( o_Body &org);
    virtual o_Body *MakeInstance();

		/*! Performs a homogeneous transfomation on all points. Calls
    ApplyTransform(1) first. */
    void	ApplyTransform(o_HMatrix &mat);

		/*! Applies the transformations, which have been added to the body
    and its subsurfaces as properties. If remove_prop is not 0, the
    properties will be removed after transformation.
    The transformation properties must be instances of
    o_ObjectTransformation . */
    void	ApplyTransform (int remove_prop = 1);

		/*! 	determines the box surrounding the surface. After
    the call p1 contains the box vertex with the lowest
    x, y and z-coordinate anf p2 the biggest coordinates. */
    void	SurBox( o_Vector3D &p1, o_Vector3D &p2 );

		//! returns the null vector
    const o_Vector3D     &Origin();

		//! returns coordinate axis n of the world coordinates
    const o_Vector3D     &Axis(int n);

		//! returns the mean of the points in point_list
    o_Vector3D Center( );

    /*! fills array containing covariance matrix of body
    array m needs to be allocated before.  */
    void CalcCoVariance(double ** m); 
		 
    //! returns area of the defined surface
    double	GetArea();

		/*! returns reference of the first triangle in body that has
    a reference to 'p' or NULL if point was not found */
    o_Triangle	*Search_Point( o_Point &p );

		/*!
    returns a hashtable pts containing all points that are are 
    common to both surfaces. If recursive==true subsurfaces are 
    also scanned
     */
    void RetrieveCommonPoints(o_HashTableSmpl<o_Point *, 
                              o_Point> & pts, o_Surface & s1, o_Surface & s2, 
                              int recursive); 

		/*!
    Fills hashtable pts with all points that are connected to 
    point pt by a triangle edge. 
     */
    void GetAdjacentPoints(o_HashTableSmpl<o_Point *, o_Point> & pts, 
                           o_Point & pt);
		
		/*!
    Search surface list and recursively all subsurfaces for a surface
    with the given name. If none is found NULL is returned.
     */
    o_Surface * SearchSurfaceByName(const char * s_name);

    virtual const   char    *GetClassName() const;

    //! returns the vector of points in the body
    const std::vector<o_Point*> & GetPointList() const { return point_list; } 
		/*! inserts point and returns reference. If point (coordinates)
    already exists n the point list the point will not inserted
    and the reference of the already inserted point will be returned */
    o_Point *AddPointChecked( o_Point &point, bool &success );
    //! adds a point to the end of the point list
    void AddPoint ( o_Point &itm );
    //! inserts point item before the list member succ
    void InsertPoint ( o_Point &itm , o_Point &succ);
    //!  substract point itm from the point list
    bool SubPoint ( o_Point &itm );
		//! returns true if p is already in the point list
    bool IsInPointList( o_Point &p );
		/*! returns pointer to first point object in the point list of
    the body with distance <= tol (that means equal coordinates in
    practice) to the given vector p. If no point exists a NULL
    pointer is returned. */
    o_Point *CheckCoordinate( o_Vector3D &p, double tol=1e-38 );
    //! returns the index of the point in the point list 
    int GetPointIndex( o_Point &p);
    //! returns the point with index n
    o_Point *GetPoint( int n );
    //! returns the number of points
    int GetNumPoints() const;
    
    /*! returns reference to o_VertexInfo object or Null if p is
    not in the point list. This contains topology information of the triangles */ 
    o_VertexInfo	*GetVertexInfo( o_Point &p) ;
    
    //! returns the surface list of the body
    const std::list<o_Surface*> & GetSurfaceList() const { return surface_list; } 
    //!  inserts surface item before the list member succ
    void InsertSurface ( o_Surface & itm, o_Surface & succ );
    //! adds surface to the end of the surface list
    void AddSurface(o_Surface & itm);
    //!  substract surface itm from the surface list
    bool SubSurface( o_Surface & itm );

		//! Sets the geometry-update-flag for each triangle of the body
    void SetTriangleGeometryUpdateFlag();

    std::list<o_NurbsSurface*> nurbssurface_list; //! the nurbssurface list of the body
    std::list<o_NurbsCurve*> nurbscurve_list; //! the nurbscurve list of the body

  protected:
    //! the surface list of the body
    std::list<o_Surface*>  surface_list; 
    //! the vector of points in the body
    std::vector<o_Point*> point_list;
    //! for quick look up this map keeps the index of the point in the point_list
    std::map<o_Point*,int> point_list_index;
    //! each point in point_list gets an o_VertexInfo
    o_HashTableSmpl<o_VertexKey,o_VertexInfo>	vertex_info;
};


#endif
