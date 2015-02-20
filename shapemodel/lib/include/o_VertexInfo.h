#ifndef	_o_vertexinfo_incl
#  define	_o_vertexinfo_incl
//
//	o_VertexInfo.h	- Header for o_VertexInfo class
//
//	Oliver Grau, Jan.1996
//

#include <vector>
using std::vector;

/*!
  \class o_VertexInfo o_VertexInfo.h
  \brief Topology information class for a body vertex
  */

class	o_Body;

typedef	o_Point	*	o_VertexKey;

class	o_VertexInfo {
	public:
                //! return number of triangle which have the common vertex
		int TriangleListLen ( ) { return trilist.size(); }

		//! return triangle 
		o_Triangle	*GetTriangle( int n) {
			if(n>=0 && n<TriangleListLen() ) return trilist[n];
			return NULL;
		}

	private:
		friend class o_Body;
		friend class o_Triangle;

    std::vector<o_Triangle*>	trilist;

};

#endif
