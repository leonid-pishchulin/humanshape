
/* SCCS-ID: $Id: svector.h 431 2007-09-26 00:26:07Z rhys $ (c) Oliver Grau */



#ifndef	_SVECTOR_
#  define	_SVECTOR_

using namespace std;
#include <iostream>

#include <stdio.h>
//#include <math.h>

# include	"o_Basis.h"
# include "o_Notify.h"
# define        error_handler(l,m)      { \
			o_Notify(O_NFY_WARN,"%s line:%d :%s",__FILE__,__LINE__ , m); }

#ifndef	nil	
# define	nil	(0)
#endif
#ifndef	true	
# define	true	(1)
#endif
#ifndef	false	
# define	false	(0)
#endif

//forward declarations
class o_Vector3D;
class o_Vector2D;

namespace leda {
  
  class vector
    {
    protected:
      friend class matrix;
      
      double* v;
      int d;
      
      void check_dimensions(const vector&) const;
      
    public:
      
      vector(int=0); 
      vector(double, double);
      vector(double, double, double);
      vector(const vector&);
      vector(const o_Vector3D &);
      vector(const o_Vector2D &);
      ~vector() { if (v) delete[] v; }
      
      
      double length() const;
      
      int    dim()    const { return d; }
      vector norm()   const { return *this/length(); }
      
      double angle(const vector&) const; 
      
      vector& operator=(const vector&);
      
      double& operator[](int);
      
      double  operator[](int) const;
      
      vector  operator+(const vector&) const;
      vector  operator-(const vector&) const;
      vector  operator*(double)        const;
      vector  operator/(double)        const;
      double  operator*(const vector&) const;
      
      int     operator==(const vector&) const;
      int     operator!=(const vector& w)  const { return !(*this == w); }
      
      /*
	friend vector operator*(double f, const vector& v);
	friend vector operator/(const vector& v, double f);
      */
      
      
      friend std::ostream& operator<<(std::ostream& o, const vector& v);
      friend std::istream& operator>>(std::istream& i, vector& v);
      
#ifdef	notdef
      friend void Print(const vector& v, std::ostream& out=cout) { out << v; } 
      friend void Read(vector& v, std::istream& in=cin)  { in >> v; }
      
      friend void   Init(vector& x)          { x = vector(0); }
      friend void   Clear(vector& x)         { delete (vector*)&x; }
      friend int compare(const vector&, const vector&) { return 0; }
#endif
      
    };
  
}   /* eof namespace leda */

#endif



