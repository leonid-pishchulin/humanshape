#ifndef	_o_Plane
#  define	_o_Plane
//
//	o_Plane.h	- Header for oGeM Plane class
//
//	Oliver Grau, Jan./Mar.1993
//


#include	"o_StraightLine.h"

/*!
  \class o_Plane o_Plane.h
  \brief Plane class
*/

class	o_Plane  {
	public:
                /*! note that until the Plane parameters are successfully 
		  set the Plane is set to: x*(0,0,1)T = 1 (Hesse form) */
		o_Plane( );

		~o_Plane( );

		/*! returns coefficients of the coordinate representation
		  of the plane: d = a*x + b*y + c*z */
		void	Get_CoordinateCoef( double &a, double &b,
					double &c, double &d ) const;

		//! returns normal vector of the plane
		const o_Vector3D	&Normal() const;
		
		/*! returns distance from coordinate origin to the plane. Note
		that X * Normal() = D() forms the Hesse form of the plane.
		D() must not be positive */
		double	D() const;

		//! returns the distance of a point p to the plane
		double	Distance(const o_Vector3D &p) const;

		/*! returns the point of intersection between a line 
		  and the plane */
		o_Vector3D	Point_of_Intersection( const o_StraightLine &l)
					const;

		//! returns the line of intersection between the plane and p2
		o_StraightLine	Straight_of_Intersection( const o_Plane &p2)
					const;
		// update methods

		//! set Plane parameter equal to these of pl
		o_Plane& operator=(const o_Plane &pl );

		//! set plane in coordinate form ( a*x + b*y + c*z = d )
		int Set_Plane( double a, double b, double c, double d );

		/*! set plane with norm the new normal vector and d the 
		new distance */
		int Set_Plane( const o_Vector3D &norm, double d );

		/*! Set plane parameters to new plane spanned by the
		  tree points */
		int Set_Plane( const o_Vector3D &p1, const o_Vector3D &p2,
			const o_Vector3D &p3 );

	private:
		o_Vector3D	n;	// normal vector of Plane
		double	d;		// distance of the origin

};


#endif
