#ifndef	_o_world_incl
#  define	_o_world_incl
//
//	o_World.h	- Header for oGeM
//
//	Oliver Grau, Jan.1993
//

#include	<string>
#include	"oGeM.h"
#include	"o_Point.h"
#include	"o_LocalCoordinateSystem.h"
#include	"o_TopologyNode.h"

#include <list>
using std::list;
//
//      Forward declaration of backwardchained classes
//

class   o_Triangle;
class   o_Surface;
class   o_Body;
class   o_World;
class	  o_LightSource;
class   o_Camera;
class   SpxModel;

#include	"o_Body.h"

class	o_ColorCamera;
class   o_DoubleZBufCamera;

/*! 
  \class o_World o_World.h
  \brief class for world model of a scene
  
  top node of a virtual scene. Keeps a list of o_Body objects,
  o_LightSource objects and o_Camera objects.
*/

class	o_World  : public o_LocalCoordinateSystem, 
		   public o_TopologyNode
{
	public:
    
    //! returns the body list of the world
    const std::list<o_Body*> & GetBodyList() const { return body_list; } 
    //!  inserts body item before the list member succ
    void InsertBody ( o_Body & itm, o_Body & succ );
    //! adds body to the end of the body list
    void AddBody(o_Body & itm);
    //!  substract body itm from the body list
    bool SubBody( o_Body & itm );

    //! returns the lightsource list of the world
    const std::list<o_LightSource*> & GetLightSourceList() const { return lightsource_list; } 
    //!  inserts lightsource item before the list member succ
    void InsertLightSource ( o_LightSource & itm, o_LightSource & succ );
    //! adds lightsource to the end of the lightsource list
    void AddLightSource(o_LightSource & itm);
    //!  substract lightsource itm from the lightsource list
    bool SubLightSource( o_LightSource & itm );

    //! returns the lightsource list of the world
    const std::list<o_Camera*> & GetCameraList() const { return camera_list; }
    //!  inserts camera item before the list member succ
    void InsertCamera ( o_Camera & itm, o_Camera & succ );
    //! adds camera to the end of the camera list
    void AddCamera(o_Camera & itm);
    //!  substract camera itm from the camera list
    bool SubCamera( o_Camera & itm );
  
		void	Copy( o_World &org );

		o_World(char *nam ) ;
		~o_World( );
      
		/*! search body with given name. returns reference on
		  success or NULL otherwise */
		o_Body	* Search_Body ( const char *bn );

		/*! search surface with given name. returns reference on
		  success or NULL otherwise. Address: body.surf[.surfsub] */
		o_Surface	* Search_Surface ( const char *sn );

		char	*BodyNameList();
		virtual const   char    *GetClassName() const;
    
  protected:
    std::list <o_Body*> body_list;
    std::list <o_Camera*> camera_list;
    std::list <o_LightSource*> lightsource_list;
};

/*!
  \relates o_World

  Tries to read in a scene description from stp- or PANORAMA- file.
  The PANORAMA- file can contain a o_Body or a o_World object.
  If a o_Body could be read a dummy world object is created and returned.
  
  Returnvalue: pointer to allocated object on success or 0 else
*/

o_World  *o_ReadWorld( const char *fn ); 
#endif
