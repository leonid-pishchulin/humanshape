#ifndef	_o_StraightLine
#  define	_o_StraightLine
//
//	o_StraightLine.h	- Header for oGeM Straight Line class
//
//	Oliver Grau, MAY-1994
//



/*!
  \class  o_StraightLine o_StraightLine.h
  \brief StraightLine class
  */

class	o_Vector3D;

class	o_StraightLine  {
	public:
                /*! create line. Note that until parameters are set otherwise
		  the default after creation is: X = (0,0,0)T + r* (0,0,1)T */
		o_StraightLine( );

		/*! create line and set it to straight line that runs thru the
		both given points */
		o_StraightLine( const o_Vector3D &p1, const o_Vector3D &p2);

		~o_StraightLine( );

		// access methods

		/*! returns vector a and b that give parametric form of the
		  line:	X = A() + r * B() */
		const o_Vector3D  &A() const { return a; }
		const o_Vector3D  &B() const { return b; }

		//! returns the point on the line perpendicular from p
		o_Vector3D FootPointofPerpendicular( const o_Vector3D & p );

		// update methods

		//! set line equal to l
		o_StraightLine& operator=(const o_StraightLine &l );

		//! set line so that it runs thru p1 and p2
		int Set_StraightLine( const o_Vector3D &p1, const o_Vector3D &p2);
	private:
		o_Vector3D	a,b;	// two point on the line

};


#endif
