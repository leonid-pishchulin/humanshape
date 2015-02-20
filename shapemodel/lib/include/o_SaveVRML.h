#ifndef	_o_savevrml_incl_
#  define	_o_savevrml_incl_
//
//      o_SaveVRML.h  - Header for o_SaveVRML class
//
//      Oliver Grau, 18-JAN-1995
//      Oliver Grau, 20-FEB-1997
//

#include	"fileio.h"
#include	"o_World.h"
//#include	"o_SaveBody.h"
#include	"o_SaveInventorBase.h"
#include 	"o_MetaCamera.h"

/*!
  \class o_SaveVRML o_SaveVRML.h
  \brief save Body as VRML (Virtual Reality Modeling Language) file
  */

class	o_MetaCamera;

class	o_SaveVRML : public o_SaveInventorBase {
	public:
                /*! creates <A HREF="http://www.tnt.uni-hannover.de/subj/vrml/overview.html">VRML-File</A> object. \sa <A HREF="http://www.tnt.uni-hannover.de/soft/compgraph/fileformats/overview.html">fileformats</A> */
		o_SaveVRML( );

    int	Save(std::vector<o_MetaCamera*> &cl, o_World &w);
    int	Save_MultipleVertices(std::vector<o_MetaCamera*> &cl, o_World &w);
	
		// old constructor
		o_SaveVRML( o_Body &body, write_file &of );

		virtual	void	Save_All();
	private:
		void    SaveCamera( o_MetaCamera *b );
		void    Save_Body( o_Body *b );
		void    Save_Surface( o_Surface *s );
		void    Save_Polygon( o_Polygon *s );
		void    Check_Anchor( o_PropertyList *s );
		void    Check_Material( o_PropertyList *s );
		void    Check_Texture( o_PropertyList *s );

		void    Save_Surface_MultipleVertices( o_Surface *s );
		void    Save_Body_MultipleVertices( o_Body *b );
    void    Save_NurbsCurve( const o_NurbsCurve *n );
		int	camidcnt;
};

#endif

