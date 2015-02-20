//
//	Oliver Grau, APR-1996
//
//	$Id: o_Polygon.h 469 2007-10-11 02:44:31Z ben $

#ifndef	_o_Polygon_incl
#  define	_o_Polygon_incl

#include	"o_PlanarFace.h"
#include	"o_Surface.h"
#include	<vector>
#include	"o_MetaCamera.h" 

#ifndef PI
#define PI 3.14159265358979323846264338327950288419716939937511
#endif

/*!
  \class o_Polygon o_Polygon.h
  \brief 3-D Polygon class

  polygonal object in 3-D. Is internally handled like a surface
  with triangles.
 */

class	o_DataBaseOut;
class	o_DataBaseIn;

class	o_Polygon : public o_PlanarFace, public o_Surface {
  public:
    virtual	const char	*GetClassName() const;

		//! create a polygon as a sub object of surf
    o_Polygon( o_Surface &s );

    ~o_Polygon();

		/*! calls AppendPoint and takes care that the point will be
    inserted in the body point list if necassary */
    void    AddPoint( o_Point &p );

    void    AppendPoint( o_Point &p );

    /*! adds a list of point to polygon and takes care that
    the point will be inserted in the body point list if
    necessary */
    bool    AddPointList( std::vector<o_Point*> &pl );
    
		// Access functions
    virtual const o_Vector3D  &Origin() const;
    virtual const o_Vector3D  &Axis(int n) const;

		/*! Returns the underlying plane or NULL, if plane
    could not be computed. */
    const o_Plane *GetPlane();

    virtual void    Copy( o_Polygon &org);
    virtual o_Surface *MakeInstance( o_Surface &parent);

    //! constraints the normal for tesselation
    void ConstraintNormal(const o_Vector3D &normal) { m_constraint_normal=normal; m_use_constraint_normal=true; }
    //! project the polygon into the camera image plane for tesselation
    void ConstraintCamera(const o_MetaCamera *camera) { m_constraint_camera=camera; m_use_constraint_camera=true; }
    //! uses the normal of the internal calculated plane for tesselation
    void FreeConstraints() { m_use_constraint_normal=false; m_use_constraint_camera=false; }
    
  private:
    // can not handle a self-intersecting polygon or a polygon with a polygonal hole
    bool TesselatePoly();
		bool DelaunayTriangulation(vector<int> &init_tri_list, vector<int> &constrained_edge, int num_verts, vector<o_Vector2D> &vert_list);
		//helper funtions for TesselatePoly()
    
		bool  IsAnyPointInside(std::vector <o_Vector2D> &points, int i, int j, int k) const;
    bool IsPointInside(std::vector <o_Vector2D> &points, int pt, int i, int j, int k)  const;
		
		bool  IsAnyPointInside(std::vector <o_Vector3D> &points, int i, int j, int k, int nVertex) const;
    bool  IsPointInside(const o_Vector3D &point, const o_Vector3D &q2)  const;
    int TriangleArea(std::vector <o_Vector3D> &points,int i,int j,int k,const o_Vector3D& normal);
    bool LineSegmentsIntersect2D(o_Vector3D &l11, o_Vector3D &l12, o_Vector3D &l21, o_Vector3D &l22);
    //helper variables for TesselatePoly()
    enum  { pos_normal = 1, degenerate = 0, neg_normal = -1 };
    std::vector < int > m_nIndex;
    o_Vector3D m_e0;
    o_Vector3D m_e1;
    o_Vector3D m_N;
    o_Vector3D m_constraint_normal;
    bool m_use_constraint_normal;
    const o_MetaCamera *m_constraint_camera;
    bool m_use_constraint_camera;
    double m_A;
};


#endif
