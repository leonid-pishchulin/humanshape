
#ifndef	_o_metasurface_incl
#  define	_o_metasurface_incl
//
//	o_MetaSurface.h	- Header for oGeM
//
//	Oliver Grau, Jan./Mar.1993
//


#include	"o_TopologyNode.h"
#include	"o_World.h"

//
//	Forward declaration of backwardchained classes
//
class	o_Triangle;
class	o_Surface;
class	o_Body;
class	o_World;
class	o_Property;
class	o_Line;
/*!
  \class o_MetaSurface o_MetaSurface.h
  \brief Meta class for surface classes (do not instantiate)
*/

class	o_MetaSurface : public o_LocalCoordinateSystem, public o_TopologyNode {
#ifdef WIN32
	protected:
#else
	private:
#endif
		o_Body	*bodyref;
		o_Surface	*surfaceref;
		int geometry_update_flag;
	public:
		virtual ~o_MetaSurface();

		//! copy + instancing
		void    Copy( o_MetaSurface &org);

		/*! returns a value unequal 0 if line is a line segment
		of the surface boundary curve */
		virtual	int	IsBoundary( o_Line &lin );

		/*! get surface property 'propertyname'. If a surface contains not
		the property then the function searches in the parent surfaces.
		If no property with the name is found it returns NULL.
		So there is a inheritance of properties. See also:
		o_PropertyList::GetLocalProperty() */
		o_Property	*GetProperty( char *propertyname );

		//! returns reference to parent surface if any or Null otherwise
		o_Surface *GetSurfaceRef() const;

		//! returns reference to body, object belongs to
		o_Body *GetBodyRef() const;

		//! 	Returns 1 if geometry of the object has changed, otherwise 0.
		int GetGeometryUpdateFlag() const 
		  	{return geometry_update_flag;}

		/*! Sets an internal flag, indicating that the geometry of the
		object has changed. */
		void SetGeometryUpdateFlag() {geometry_update_flag = 1;}

	protected:
		o_MetaSurface( o_Body &b, char *nam ) ; 
		o_MetaSurface( o_Surface &s, char *nam ) ; 
		void	Set_SurfaceRef( o_Surface &b );
		void	Set_BodyRef( o_Surface &b );
		void	Set_BodyRef( o_Body &b );
		void ResetGeometryUpdateFlag() {geometry_update_flag = 0;}
};


#endif
