#ifndef	_o_point_incl
#  define	_o_point_incl
//
//	o_Point.h	- Header for the point class
//
//	Oliver Grau, Jan.1993
//

#include	"oGeM.h"
#include	"o_PropertyList.h"

/*!
  \class o_Point o_Point.h
  \brief 3-D Point class

  represents point in 3-D euclidean space for the use as vertices 
  ( o_Triangle , o_Line ). The coordinates are given in [mm].
  Objects of the class o_Point must instantiated dynamicly. The 
  destructor is not public. Use ReqDelete instead of for deleting.
*/

class   o_HMatrix;

class	o_Point : public o_Vector3D, public o_PropertyList {
	public:
                //! creates Point with cartesian coordinates (x,y,z)
		o_Point(double x=0, double y=0, double z=0 ) ;

		//! creates Point from (3 dimensional) o_Vector3D
		o_Point(const o_Vector3D &vec ) ;

		o_Point(o_Point &vec ) ;
		void Copy(o_Point &org);

		/*! sends a request for deleting the point. If there
		are no links to other oGeM-objects the delete
		operator will be called. The delete operator must not
		be called for this class objects directly. The
		method returns 1 if the point was deleted */
		int ReqDelete( o_MetaClass * optr );

		//! transforms a point to new position with homogenious matrix
		virtual void ApplyTransform(o_HMatrix &matrix);

		friend class o_HMatrix;
		~o_Point() ;
  
		private:   // Disabled operator=
		  o_Point &operator=( const o_Point & );
	
};

#endif
