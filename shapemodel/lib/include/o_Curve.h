#ifndef	_o_curve_incl_
#  define	_o_curve_incl_
//
//      o_Curve.h  - Header for curve class
//
//      Oliver Grau, 10-JUN-1993
//	Oliver Grau, OCT-1996	list interface for points
//

#include	"oGeM.h"
#include	"o_Line.h"
#include <vector>
using std::vector;

/*!
  \class o_Curve o_Curve.h
  \brief curve in 3D approximated by line segments
  */

class	o_Curve : public o_MetaClass {
	public:

		//! creates curve with optional name
		o_Curve(char *name=NULL);

		const   char    *GetClassName() const;

    //! returns the vector of points in the curve
    const std::vector<o_Point*> & GetPointList() const { return point_list; } 
    void AddPoint ( o_Point &itm );
    //! inserts point item before the list member succ
    void InsertPoint ( o_Point &itm , o_Point &succ);
    //!  substract point itm from the point list
    bool SubPoint ( o_Point &itm );
    //! returns the point with index n
    o_Point * GetPoint(int n);
    
	protected:
    std::vector<o_Point*> point_list;
};

#endif	/* _o_curve_incl_ */

