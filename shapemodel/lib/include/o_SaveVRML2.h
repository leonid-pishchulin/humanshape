#ifndef	_o_savevrml2_incl_
#  define	_o_savevrml2_incl_
//
//      o_SaveVRML2.h  - Header for o_SaveVRML2 class
//
//      Oliver Grau, 20-FEB-1997
//

#include	"fileio.h"
#include	"stdlib.h"
#include	"o_World.h"
//#include	"o_SaveBody.h"
#include	"o_SaveInventorBase.h"
#include	"o_TransformationList.h"
#include 	"o_MetaCamera.h" 
#include <vector>
using std::vector;
//
//	helper class for animations
//
class   o_AnimationList ;
#include <vector>

/*!
  \class o_SaveVRML2 o_SaveVRML2.h
  \brief save Body as VRML2 (Virtual Reality Modeling Language) file
  */

class	o_MetaCamera;
class	o_Material;
class	o_NurbsCurve;

class	o_SaveVRML2 : public o_SaveInventorBase {
	public:
                /*! creates <A HREF="http://www.tnt.uni-hannover.de/subj/vrml/overview.html">VRML2-File</A> object. \sa <A HREF="http://www.tnt.uni-hannover.de/soft/compgraph/fileformats/overview.html">fileformats</A> */
		o_SaveVRML2( );

		~o_SaveVRML2( );

		/*! save a list of viewpoints generated from the list of cameras
		and a world object to VRML 2 file */

    int Save(std::vector< o_MetaCamera* >  &cl, o_World &w);

		// translate name
		const	char *VRMLName(o_MetaClass &o);
		void	SetRedundantMaterial(int m);

		virtual	void	Save_All();
	protected:
		virtual void    Info(const char *);
		virtual void    Tab();
		void	StartTransform(const char *nam);
		void	EndTransform(const char *nam=0);

		// return 1 if local texture map/ material
		int	UpdateAppearance(o_PropertyList *b);

    void	Save_TextureBindings(std::vector<o_Triangle*> &tlist,std::vector<o_Vector2D*> &blist );

		void    SaveCamera( o_MetaCamera *b );
		void    Save_Body( o_Body *b );
		void    Save_Surface( o_Surface *s );
		void    Save_Mesh( o_Surface *s );
		void    Check_Anchor( o_LocalCoordinateSystem *s );
		void    Check_AnimationList( o_LocalCoordinateSystem *s );
		void    Check_Texture( o_LocalCoordinateSystem *s );
		void    Check_Material( o_Body *s );
		void    Check_ObjectTransformation( o_PropertyList *s );
		void    Check_Material( o_Surface *s );
		void    SaveShapeHints( o_Surface *s ); // in IFS for meshes
		void    SaveMaterial( o_PropertyList *s, o_Material      *mat );
		void    Check_Texture( o_Surface *s );
    void    Save_NurbsCurve( const o_NurbsCurve *n );
		void	StoreAnimationParameter( );
		void	RouteAnimationParameters( );

		int	camidcnt;
		//int	texfileformat;	// 0=jpg, 1=tiff, 2=rgb;
		int	insertdefmaterial;	// if true always include default
						// material

    std::vector<o_AnimationList*>	alist;
};

/*!
  \class o_AnimationList o_SaveVRML2.h
  \brief helper class for animations
  */

class   o_AnimationList {
	public:
		o_TransformationList    *tl;
		o_LocalCoordinateSystem *key;
		char	*idnam;
		o_SaveVRML2	*sp;
		void CreateIdName(o_LocalCoordinateSystem *s ) {
			int l= strlen(sp->VRMLName(*s));
			idnam = new char[l+3];
			strcpy(idnam,"AL");
			strcat(idnam,sp->VRMLName(*s));
		}
		o_AnimationList(o_SaveVRML2     *savp,
				o_LocalCoordinateSystem *s) { 
			sp = savp;
			CreateIdName(s);
			key=s;tl=0;
		}
		~o_AnimationList() { if(idnam) delete [] idnam;}
};

#endif

