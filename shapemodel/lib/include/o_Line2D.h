#ifndef _o_Line2D_incl
#define _o_Line2D_incl

#include "o_Vector.h"

/*!
  \class o_Line2D o_Line2D.h
  \brief 2-D line segment
  */

class	o_Line2D {
	public:
		enum CClass {
			Undetermined=0, Silhouette, C0Discontinuous,
			C0Continuous, C1Continuous
		};

		//! create line segment
		o_Line2D( ) {segclass=Undetermined;};

		//! create line segment with values
		o_Line2D( const o_Vector2D &start, const o_Vector2D &end ) {
			s = start; e = end;
		};

		~o_Line2D() {};

		//! return coordinates of start point of the line
		o_Vector2D  GetStartPoint() const { return s; }
		//! return coordinates of end point of the line
		o_Vector2D  GetEndPoint() const { return e; }

		/*! return continuous class. CClass is a enum type with the
		  following meaning:
		  Undetermined, Silhouette, C0Discontinuous,
		  C0Continuous, C1Continuous */
		CClass  GetCClass() const { return segclass; }

		//! set the continuous class of the line segment
		void SetCClass(CClass isegclass) { segclass = isegclass; }

		//! set coordinates of start point of the line
		void SetStartPoint(o_Vector2D &start) { s = start; }
		//! set coordinates of end point of the line
		void SetEndPoint(o_Vector2D &end) { e = end; }
	private:
		CClass	segclass;
		o_Vector2D	s;
		o_Vector2D	e;
};

#endif
