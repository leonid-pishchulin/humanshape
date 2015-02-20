
#ifndef	_o_surface_incl
#  define	_o_surface_incl
//
//	o_Surface.h	- Header for oGeM
//
//	Oliver Grau, Jan./Mar.1993
//


#include	"o_MetaSurface.h"
#include	"o_World.h"
#include  "o_HashTable.h" 


//
//	Forward declaration of backwardchained classes
//
class o_Camera;
class	o_Triangle;
class	o_Surface;
class	o_Body;
class	o_World;
class	o_PoseParameterSurface;

/*!
  \class o_Surface o_Surface.h
  \brief Surface which can contain sub-surfaces and triangles
 */

class	o_Surface : public o_MetaSurface	{
  public:
    //! returns the list of sub surfaces in the surface
    const std::list<o_Surface*> & GetSurfaceList() const { return surface_list; } 
    //!  inserts surface item before the list member succ
    void InsertSurface ( o_Surface & itm, o_Surface & succ );
    //! adds surface to the end of the surface list
    void AddSurface(o_Surface & itm);
    //!  substract surface itm from the surface list
    bool SubSurface( o_Surface & itm );
    
    //! returns the list of triangles in the surface
    const std::list<o_Triangle*> & GetTriangleList() const { return triangle_list; } 
    //!  inserts triangle item before the list member succ
    void InsertTriangle ( o_Triangle & itm, o_Triangle & succ );
    //! adds triangle to the end of the triangle list
    void AddTriangle(o_Triangle & itm);
    //!  substract triangle itm from the triangle list
    bool SubTriangle( o_Triangle & itm );

		// creators

    //! create surface without root node spezified
    o_Surface( ) ; 

    //! assign surface with name nam to body b
    o_Surface( o_Body &b, char *nam ) ; 

    //! creates surface as a sub-surface
    o_Surface( o_Surface &b, char *nam ) ;
 
    ~o_Surface( );

		// copy + instancing
    virtual void    Copy( o_Surface &org);
    virtual void    CopyWSS( o_Surface &org );
    virtual o_Surface *MakeInstance( o_Surface &parent);
    virtual o_Surface *MakeInstance( o_Body &parent);

    virtual const   char    *GetClassName() const;

    /*! determines the bounding box of the surface. After
    the call p1 contains the box vertex with the lowest
    x, y and z-coordinate and p2 the biggest coordinates. */
    void	SurBox( o_Vector3D &p1, o_Vector3D &p2 );

		//! returns area of the defined surface
    double	GetArea();

    /*! like Surbox() but does not initialize p1 and p2. Is used
    for recursive determination of bounding box */
    void	getSurBox( o_Vector3D &p1, o_Vector3D &p2 );

    /*! returns reference of the first polygon in the surface or
    subsurfaces that has a reference to 'p' or NULL if point was not found */
    o_Triangle *Search_Point( o_Point &p );

    /*! connect the surface to a body,
    inserts surface item before the list member succ (if given)
    or appends to the surfacelist otherwise */
    void	ConnectRootNode(o_Body &b, o_Surface *succ=NULL);

    /*! connect the sub-surface to a surface,
    inserts surface item before the list member succ (if given)
    or appends to the surfacelist otherwise */
    void	ConnectRootNode(o_Surface &s, o_Surface *succ=NULL);

		/*! Sets the geometry-update-flag for each triangle of the
    surface and its subsurfaces */
    void SetTriangleGeometryUpdateFlag();

    /*! returns 0 for normal o_Surface object and 1 for switch
    objects (e.g. o_LevelOfDetail ) */
    virtual int IsSwitch();
		
		//! set the rigid flag
    void SetRigidConnectedFlag(int val);

    /*! indicates whether the surface (and it's sub surfaces) are
    connected rigidly or not. */
    int GetRigidConnectedFlag();

		/*! returns hashtable  containing all points referenced by the 
    surface mesh. If recursive==true subsurfaces are also scanned */
    void RetrievePoints(o_HashTableSmpl<o_Point *, 
                        o_Point> & pts, int recursive);

		/*! takes all triangles from surface ng_surf that are connected
    to this surface and moves them to surface dest_surf. Does not
    duplicate any points, used for flexible connections between
    adjacent surfaces
     */
    void MakeFlexibleConnection(o_Surface & ng_surf, 
                                o_Surface & dest_surf);

		/*!
    Search subsurface list and recursively all subsurfaces for a surface
    with the given name. If none is found NULL is returned.
     */
    o_Surface * SearchSurfaceByName(const char * s_name);

    //! fill tri_list with the triangles from the surface and all sub-surfaces
    void SampleTriangles(std::list<o_Triangle*> *tri_list);

  protected:
    friend class o_Polygon;
    friend class o_Triangle;
    friend class o_Body;
    friend class o_PoseParameterSurface;
    int	rigid_flg;
    std::list<o_Surface*> surface_list;	//! list of sub-surfaces
    std::list<o_Triangle*> triangle_list;	//! list of triangles
    void ApplyTransformIncr (int remove_prop);
};

/*!
  \relates o_Surface

  apply all (sub-)/surfaces to a function. 

  The function afuncp is applied to all surfaces of bp.
  (int *afuncp(o_Surface *)) : afunc is the pointer of a user
  defined function. The function
  is called for each surface, which is parameter of that
  function. If the function returns a zero ApplyToAllSurfaces()
  steps not further in deeper subsurfaces
 */

void    ApplyToAllSurfaces( o_Body *bp, int (*afuncp)(o_Surface *sp) );

/*!
  \relates o_Surface
  \overload
 */

void    ApplyToAllSurfaces( o_Surface *sp, int (*afuncp)(o_Surface *sp) );


#endif
