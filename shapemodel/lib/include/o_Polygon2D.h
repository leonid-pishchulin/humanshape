#ifndef	_o_Polygon2D_incl
#  define	_o_Polygon2D_incl

#include	"oGeM.h"
#include <vector>
using std::vector;

/*!
  \class o_Polygon2D o_Polygon2D.h
  \brief 2d Polygon class
  */

class	o_Polygon2D : public o_MetaClass {
	public:
		virtual	const char	*GetClassName() const;

		//! create a 2d polygon
		o_Polygon2D();

		~o_Polygon2D();

    o_Vector2D *GetPoint(int n);
		o_Vector2D &Point(int n) ;
		o_Vector2D &operator [](int n) ;

		//! returns the number of points in the point list
		int NumberOfPoints() const;

		/*! appends a point to the polygon. The point p will be
		  copied by this routine and deleted by the destructor */
		void AppendPoint( const o_Vector2D &p );

		//! deletes point n from the point list
		void DeletePoint( int n );
	private:
    std::vector<o_Vector2D*>	plist;
};


#endif
