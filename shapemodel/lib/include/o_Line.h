#ifndef	_o_line_incl_
#  define	_o_line_incl_
//
//      o_line.h  - Header for line class
//
//      Oliver Grau, 16-JUN-1993
//

#include	"oGeM.h"
#include	"o_Point.h"

/*!
  \class o_Line o_Line.h
  \brief line of regions or surfaces
  */

class	o_Curve;

class	o_Line : public o_MetaClass {
	public:
                /*! creates line from start and end point. The point must
		  be instantiated by new. The instantiated line object
		  connects to the point reference */
		o_Line( o_Point	&start, o_Point &end ) :o_MetaClass(NULL){
			start.Connect( this );
			end.Connect( this );
			p1 = &start;
			p2 = &end;
		};

		/*! on deletion of the object the point references will be
		  disconnected. This means the points corresponding to the 
		  references will be deleted also if no other reference is
		  established ! See the o_Point class documentation */
		~o_Line() ;

		//! returns start point of the line
		o_Point	&Start_Point() { return *p1; }
		//! returns start point of the line
		o_Point	&End_Point() { return *p2; }

	protected:
		friend class o_Curve;
		o_Point	*p1;
		o_Point	*p2;

		//! only for internal use
		o_Line( ) :o_MetaClass(NULL){
			p1=p2=NULL;
		};
};



#endif	/* _o_line_incl_ */

