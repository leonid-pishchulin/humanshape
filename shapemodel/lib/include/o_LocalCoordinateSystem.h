#ifndef	_o_spatial_incl
#  define	_o_spatial_incl
//
//	o_LocalCoordinateSystem.h	- Header for oGeM
//
//	Oliver Grau, May.1993
//

#include	"oGeM.h"
#include	"o_PropertyList.h"

class o_HMatrix;

/*!
  \class o_LocalCoordinateSystem o_LocalCoordinateSystem.h
  \brief meta class for extended (non point) datatypes - contains offset and orientation

*/

class	o_LocalCoordinateSystem : public o_PropertyList {
	public:
                //! returns the offset of the object in the world coordinate system
		virtual	const	o_Vector3D	&Origin() const;

		//! returns coordinate axis n of the object. n must be in {0,1,2}
		virtual	const	o_Vector3D	&Axis(int n) const;

		/*! returns the position of a point given in world coordinates
		in a relative coordinates (given by Origin() an Axis(n)) */
		o_Vector3D	RelPosition( const o_Vector3D &worldcoord) const;

		/*! returns the position of a point in world coordinates
		  relative to the coordinate system given by Origin() an Axis(n) */
		o_Vector3D	WorldPosition( const o_Vector3D &localcoord ) const;

		/*! Get the homogeneous matrix describing the transformation
		from relative to world coordinates. */
		void 		WorldTransformation (o_HMatrix &) const;
		
		/*! Get the homogeneous matrix describing the transformation
		from world to relative coordinates. */
		void		RelTransformation (o_HMatrix &) const;

		void    Copy( o_LocalCoordinateSystem &org );

	protected:
		o_LocalCoordinateSystem();
		//virtual	~o_LocalCoordinateSystem();
};


#endif




