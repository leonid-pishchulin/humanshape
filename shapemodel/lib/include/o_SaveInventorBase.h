#ifndef	_o_saveivbase_incl_
#  define	_o_saveivbase_incl_
//
//      o_SaveInventorBase.h  - Header for o_SaveInventorBase class
//
//      Oliver Grau, 18-JAN-1995
//      Oliver Grau, 20-FEB-1997
//

#include	"fileio.h"
#include	"o_World.h"
#include	"o_SaveBody.h"


/*!
  \class o_SaveInventorBase o_SaveBody
  \brief save Body. Base class for inventor derivates (Inventor, VRML)

  abstract class do not instantiate
  */

class	o_MetaCamera;
class	o_PointTable;
class	o_TextureMap;

class	o_SaveInventorBase : public o_SaveBody {
	protected:
		o_SaveInventorBase( );
		virtual ~o_SaveInventorBase( );

		// old constructor
		o_SaveInventorBase( o_Body &body, write_file &of );
	public:
		virtual	void	Save_All();

		// parameters for texture mapping

		/*! set a name mapping format (c-print format). The default is:
		"map_%d.%s". Where %d is the current map number starting with 0.
		%s is the file extension. */
		void	SetMapNameFormat( char *path );

		char	*GetMapNameFormat( );

		//! set path for external modules (texture maps)
		void	SetPath( char *path );

		//! returns path for external modules 
		char	*GetPath( );

    //! returns filename for external modules 
    char	*GetFileName( );
    
		enum	TextureMapFormat { NoMap, JPEG, Tiff };

		/*! set texture map file format. Supported formats:
		{ NoMap, JPEG, Tiff }. On NoMap no texture map is
		written out. */
		void    SetTextureFileFormat(TextureMapFormat ff);

		/*! if set true the filter writes one separate point list
		for each (sub-)surface. If false only one common list for
		each body is written. Default: true */
		void    SetSurfaceWisePointListFlag(int m);

		// ShapeHints
		enum VERTEX_ORDERING { UNKNOWN_ORDERING,CLOCKWISE,COUNTERCLOCKWISE };
		enum SHAPE_TYPE { UNKNOWN_SHAPE_TYPE,SOLID };
		enum FACE_TYPE { UNKNOWN_FACE_TYPE, CONVEX };

		/*! set vertex ordering.
		  <P>
		  enum VERTEX_ORDERING { UNKNOWN_ORDERING,CLOCKWISE,
		  COUNTERCLOCKWISE }:</P>
                  <P>
		  <UL>
		  <LI>UNKNOWN_ORDERING    Ordering of vertices is unknown</LI>
		  <LI>CLOCKWISE           Face vertices are ordered clockwise
		  (from the outside)</LI>
		  <LI>COUNTERCLOCKWISE    (default) Face vertices are ordered
		  counterclockwise (from the outside)</LI>
		  </UL>
		  </P> */
		void	SetVertexOrdering (VERTEX_ORDERING vo);

		VERTEX_ORDERING GetVertexOrdering();

		/*! set shape type. enum FACE_TYPE { UNKNOWN_FACE_TYPE, CONVEX }.
		default: CONVEX */
		void	SetShapeType (SHAPE_TYPE vo);

		SHAPE_TYPE GetShapeType();

		/*! set the crease angle.The crease angle is the angle between
		surface normals on adjacent polygons.
		For example, a crease angle of .5 radians (the default value)
		means that an edge between two adjacent
		polygonal faces will be smooth shaded if the normals to the two
		faces form an angle that is less than .5 radians
		(about 30 degrees). Otherwise, it will be faceted */
		void	SetCreaseAngle (double vo);
		
		//! get values of the shape hints (see Update functions)
		double GetCreaseAngle();


		// obsolete
		int	redundant_coord;
	protected:
		virtual void    SaveCamera( o_MetaCamera *c );
		virtual void    Save_Body( o_Body *b );
		virtual void    Save_Surface( o_Surface *s );
		virtual void    Save_Polygon( o_Polygon *s );
		virtual void    SaveShapeHints();
		virtual void    Check_Material( o_PropertyList *s );
		virtual void    Check_ObjectTransformation( o_PropertyList *s );
		virtual void    Check_Texture( o_PropertyList *s );
		// auxliary
		char	*SaveTextureMap(o_TextureMap *map);
    void  Save_IndexedFaceSet(std::vector<o_Triangle*> &tl, o_PointTable &plist);
    void  Save_TextureBindings(std::vector<o_Triangle*> &tl, std::vector<o_Vector2D*> &blist);
		virtual void	Tab();
		virtual void	Info(const char *);

    char	*MapFileName(char *ext, bool only_filename=false); // returned string must be
					// deleted after use !!!!!

		// ShapeHints
		VERTEX_ORDERING	ordering;
		SHAPE_TYPE	shapetype;
		FACE_TYPE	facetype;
		double	creaseangle;

		int	level;	// for tabulators
		int	index;	// for number for texture/point coordinates
		int	texmapno;	// for number for texture maps
		char	*epath;
		char	*mapnamefmt;
		TextureMapFormat     texfileformat;  // 0=jpg, 1=tiff, 2=rgb;
		int     surface_wise_coords;
};

#endif

