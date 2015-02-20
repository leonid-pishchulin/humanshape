#ifndef	_oGeMvector_incl
#  define	_oGeMvector_incl
//
//	o_Vector.h	- Header for oGeM -vectors
//
//	Oliver Grau, Mar.1993
//
#ifdef WIN32
#include <stdlib.h>
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#endif
#include <iostream>
#include <fstream>
#include "svector.h"
#include "smatrix.h"
#include "o_Notify.h"
using namespace std;  

class	o_HVector;

/*!
  \class o_Vector3D o_Vector.h
  \brief 3-D Vector class
*/

class	o_Vector3D {
 public:
  //! creates 3d vector with optional content (x,y,z).
  o_Vector3D(double x=0, double y=0, double z=0 );
  
  //! creates 3d vector from (LEDA)-vector.
  o_Vector3D(const leda::vector &vec );
  
  o_Vector3D(const o_HVector &vec );
  
  //! the copy constructor.
  o_Vector3D (const o_Vector3D &);
  
  ~o_Vector3D() ;
  
  //! returns x-component of the vector.
  double	GetX() const { return data[0]; }
  
  //! returns y-component of the vector.
  double	GetY() const { return data[1]; }
  
  //! returns z-component of the vector.
  double	GetZ() const { return data[2]; }
  
  //! assigns x-component of the vector with the value
  void	SetX( const double x) { data[0] = x; }
  
  //! assigns y-component of the vector with the value
  void	SetY( const double y) { data[1] = y; }
  
  //! assigns z-component of the vector with the value
  void	SetZ( const double z) { data[2] = z; }
  
  //! returns the length of the vector.
  double Length() const;
  
  //! returns the normalized vector.
  o_Vector3D Norm() const;
  
  //! returns 3.
  int Dimension() const;
  
  //! returns the angle between *this and vec in RAD.
  double angle (const o_Vector3D &) const;
  
  //! returns the ith element as an LValue.
  double &operator[](int i);
  
  double operator[](int i) const;
  
  o_Vector3D operator+ (const o_Vector3D &) const;
  o_Vector3D operator- (const o_Vector3D &) const;
  o_Vector3D operator* (double) const;
  o_Vector3D operator/ (double) const;
  double operator* (const o_Vector3D &) const;
  int operator== (const o_Vector3D &) const;
  int operator!= (const o_Vector3D &) const;
  
  friend std::ostream& operator<< (std::ostream &, const o_Vector3D &);
  friend std::istream& operator>> (std::istream &, o_Vector3D &);
  
  //! assigns vector the contens of an other vector
  o_Vector3D& operator=(const o_Vector3D&);

  o_Vector3D& operator=(const leda::vector&);   
  
  //! add vector to this
  o_Vector3D& operator +=(const o_Vector3D&);   
  
    //! substract vector to this
  o_Vector3D& operator -=(const o_Vector3D&);
  
  //! returns vector product of *this and vec.
  o_Vector3D Cross(const o_Vector3D&) const;
  
  /*!
  \brief returns perpendicular vector
    
    In three dimensions, there are an infinite number of vectors perpendicular to a given vector.
    This function returns a reasonable choice.
  */
  o_Vector3D Perp() const;
  
  double	data[3];
};

/*!
  \class o_HVector o_Vector.h
  \brief homogenous 3-D vector class
  
  class for homogenous 3-D vectors of the form: (x,y,z,w)T. The
  component w is normalized to 1.
*/

class	o_HVector {
 public:
  //! creates homogenous 3d vector with optional content (x,y,z,w).
  o_HVector(double x=0, double y=0, double z=0, double w=1 ) ;
  
  //! creates homogenous 3d vector from (LEDA)-vector.
  o_HVector(leda::vector &vec, double w=1 ) ;
  
  //! creates homogenous 3d vector from a 3d vector.
  o_HVector(const o_Vector3D &vec, double w=1 ) ;
  
  //! the copy constructor
  o_HVector (const o_HVector &vec);
  
  //! returns x-component of the vector
  double	GetX() const { return data[0]; }
  
  //! returns y-component of the vector
  double	GetY() const { return data[1]; }
  
  //! returns z-component of the vector
  double	GetZ() const { return data[2]; }
  
  //! returns w-component of the vector
  double	GetW() const { return data[3]; }
  
  //! assigns x-component of the vector with value x
  void	SetX( const double x) { data[0] = x; }
  
  //! assigns y-component of the vector with value y
  void	SetY( const double y) { data[1] = y; }
  
  //! assigns z-component of the vector with value z
  void	SetZ( const double z) { data[2] = z; }
  
  //! assigns w-component of the vector with value w
  void	SetW( const double w) { data[3] = w; }
  
  ~o_HVector() ;
  o_HVector& operator=(const o_HVector&);
  o_HVector operator+ (const o_HVector &) const;
  o_HVector operator- (const o_HVector &) const;
  o_HVector operator* (double) const;
  o_HVector operator/ (double) const;
  double operator* (const o_HVector &) const;
  o_HVector& operator +=(const o_HVector&);
  o_HVector& operator -=(const o_HVector&);
  double &operator[](int i);
  double operator[](int i) const;
  
  friend std::ostream& operator<< (std::ostream &, const o_HVector &);
  friend std::istream& operator>> (std::istream &, o_HVector &);
  
 protected:
  double	data[4];
};

/*!
  \class o_Vector2D o_Vector.h
  \brief 2-D Vector class
*/

class	o_Vector2D {
 public:
  //! creates 2d vector with optional content (x,y).
  o_Vector2D(double x=0, double y=0 );
  
  //! creates 2d vector from (LEDA)-vector.
  o_Vector2D(const leda::vector &vec );
  
  //! the copy constructor.
  o_Vector2D (const o_Vector2D &);
  
  ~o_Vector2D() ;
  
  //! returns x-component of the vector
  double	GetX() const { return data[0]; }
  
  //! returns y-component of the vector
  double	GetY() const { return data[1]; }
  
  //! assigns x-component of the vector with the value x
  void	SetX( const double x) { data[0] = x; }
  
  //! assigns y-component of the vector with the value y
  void	SetY( const double y) { data[1] = y; }
  
  //! returns the length of the vector.
  double Length() const;
  
  //! returns the normalized vector.
  o_Vector2D Norm() const;
  
  //! returns 2.
  int Dimension() const;
  
   //! returns perpendicular vector
  o_Vector2D Perp() const;
  
  //! returns the angle between *this and vec in RAD.
  double angle (const o_Vector2D &) const;

  /*!
     \brief returns the angle of vec relative to *this

     acos will usually return a value between 0 and PI (in RAD).
     If we want a + or - value to indicate which vector is ahead,
     then we need to use the atan2 function
  */
  double RelativeAngle (const o_Vector2D &) const;
  
  double &operator[] (int i);

  //! returns the ith element as an LValue.
  double operator[] (int i) const;

  o_Vector2D operator+ (const o_Vector2D &) const;
  o_Vector2D operator- (const o_Vector2D &) const;
  o_Vector2D operator* (double) const;
  o_Vector2D operator/ (double) const;
  double operator* (const o_Vector2D &) const;
  int operator== (const o_Vector2D &) const;
  int operator!= (const o_Vector2D &) const;

  friend std::ostream &operator<< (std::ostream &, const o_Vector2D &);
  friend std::istream &operator>> (std::istream &, o_Vector2D &);

		//! assigns vector the contens of an other vector
  o_Vector2D &operator=(const o_Vector2D&);

  o_Vector2D &operator= (const leda::vector &);
		
  double	data[2];
};


/*!
  \relates vector

  global functions dealing with vectors
*/
void	Normalize( leda::vector &v );

/*!
  \relates o_Vector3D

  global functions dealing with vectors
*/
void    Normalize(o_Vector3D &);

/*!
  \relates o_Vector2D

  global functions dealing with vectors
*/
void    Normalize(o_Vector2D &);





///////////////////////////////////////////////////////////////////////////
//
//Implementation of o_Vector3D
//
inline o_Vector3D::o_Vector3D(double x, double y, double z )
{
  data[0] = x;
  data[1] = y;
  data[2] = z;
}


inline o_Vector3D::o_Vector3D(const leda::vector &vec )
{
  if( vec.dim() != 3 )
    o_Notify(O_NFY_FATAL,
		    "o_Vector3D: Constructing vector has dimensionality unequal 3\n"
		    );
  data[0] = vec[0];
  data[1] = vec[1];
  data[2] = vec[2];
}


inline o_Vector3D::o_Vector3D(const o_HVector &vec )
{
  data[0] = vec[0];
  data[1] = vec[1];
  data[2] = vec[2];
}

inline o_Vector3D::o_Vector3D (const o_Vector3D &vec)
{
  data[0] = vec.data[0];
  data[1] = vec.data[1];
  data[2] = vec.data[2];
}

inline o_Vector3D::~o_Vector3D()
{
  //nothing to do
}

inline double o_Vector3D::operator[] (int i) const
{
  if (i < 0 || i > 2)
    o_Notify(O_NFY_FATAL, "o_Vector3D::operator[]: index out of range\n");
  return data[i];
}

inline double o_Vector3D::operator* (const o_Vector3D &vec) const
{
  double result = 0.0;
  result += data[0] * vec.data[0];
  result += data[1] * vec.data[1];
  result += data[2] * vec.data[2];
  return result;
}

inline double &o_Vector3D::operator[] (int i)
{
  if (i < 0 || i > 2)
    o_Notify(O_NFY_FATAL, "o_Vector3D::operator[]: index out of range\n");
  return data[i];
}

inline double o_Vector3D::Length() const
{
  return sqrt (*this * *this);
}

inline o_Vector3D o_Vector3D::operator/ (double a) const
{
  return o_Vector3D (data[0] / a,
		   data[1] / a,
		   data[2] / a);
}

inline o_Vector3D o_Vector3D::Norm() const
{
  return *this/Length();
}

inline int o_Vector3D::Dimension() const
{
  return 3;
}

inline o_Vector3D o_Vector3D::operator* (double a) const
{
  return o_Vector3D (data[0] * a,
		   data[1] * a,
		   data[2] * a);
}


inline double o_Vector3D::angle (const o_Vector3D &vec) const
{
  double l = Length();
  double yl = vec.Length();
  
  if ( l==0 || yl==0)
    o_Notify(O_NFY_FATAL,"angle: zero argument\n");
  
  // catch acos: DOMAIN error
  double	arg=  ((*this)*vec/(l*yl));  
  int	isone = (arg > 1.0-1e-38 );
  if(isone) return 0;
  int   isminusone = (arg < -1.0+1e-38);
  if (isminusone)
    return acos(-1.0);
  else
    return  acos( arg );
}

inline o_Vector3D o_Vector3D::operator+ (const o_Vector3D &vec) const
{
  return o_Vector3D (data[0] + vec.data[0],
		   data[1] + vec.data[1],
		   data[2] + vec.data[2]);
}

inline o_Vector3D o_Vector3D::operator- (const o_Vector3D &vec) const
{
  return o_Vector3D (data[0] - vec.data[0],
		   data[1] - vec.data[1],
		   data[2] - vec.data[2]);
}
  
inline o_Vector3D o_Vector3D :: Cross(const o_Vector3D &b) const
{
  return o_Vector3D( 
		  data[1] * b.data[2] - data[2] * b.data[1],
		  data[2] * b.data[0] - data[0] * b.data[2],
		  data[0] * b.data[1] - data[1] * b.data[0]
		  );
}



inline int o_Vector3D::operator== (const o_Vector3D &vec) const
{
  int	ret;

  if (data[0] != vec.data[0]) ret = 0;
  else if (data[1] != vec.data[1]) ret = 0;
  else if (data[2] != vec.data[2]) ret = 0;
  else ret = 1;
  return ret;
}

inline int o_Vector3D::operator!= (const o_Vector3D &vec) const
{
  return !(*this == vec);
}

inline o_Vector3D &o_Vector3D::operator= (const leda::vector &vec)
{
  if (vec.dim() == 3)
    {
      data[0] = vec[0];
      data[1] = vec[1];
      data[2] = vec[2];
    }
  else
    o_Notify(O_NFY_WARN, "o_Vector3D::operator=(): dimension mismatch\n");
  return *this;
}

// quick copy O.Grau 5.2.99
inline o_Vector3D &o_Vector3D::operator= (const o_Vector3D &vec)
{
  double	*v1 = &data[0];
  double	const *v2 = &(vec.data[0]);
  	
  *v1++ = *v2++;
  *v1++ = *v2++;
  *v1 = *v2;
  return *this;
}

// quick add O.Grau 5.2.99
inline o_Vector3D &o_Vector3D::operator += (const o_Vector3D &vec)
{
  double	*v1 = &data[0];
  double	const *v2 = &(vec.data[0]);
  	
  *v1++ += *v2++;
  *v1++ += *v2++;
  *v1 += *v2;
  return *this;
}

inline o_Vector3D &o_Vector3D::operator -= (const o_Vector3D &vec)
{
  double	*v1 = &data[0];
  double	const *v2 = &(vec.data[0]);
  
  *v1++ -= *v2++;
  *v1++ -= *v2++;
  *v1   -= *v2;
  return *this;
}


///////////////////////////////////////////////////////////////////////////
//
//Implementation of o_Vector2D
//
inline o_Vector2D::o_Vector2D(double x, double y)
{
  data[0] = x;
  data[1] = y;
}


inline o_Vector2D::o_Vector2D(const leda::vector &vec )
{
  if( vec.dim() != 2 )
    o_Notify(O_NFY_FATAL,
		    "o_Vector3D: Constructing vector has dimensionality unequal 2\n"
		    );
  data[0] = vec[0];
  data[1] = vec[1];
}


inline o_Vector2D::o_Vector2D (const o_Vector2D &vec)
{
  data[0] = vec.data[0];
  data[1] = vec.data[1];
}

inline o_Vector2D::~o_Vector2D()
{
  //nothing to do
}

inline double o_Vector2D::operator[] (int i) const
{
  if (i < 0 || i > 1)
    o_Notify(O_NFY_FATAL, "o_Vector3D::operator[]: index out of range\n");
  return data[i];
}

inline double &o_Vector2D::operator[] (int i)
{
  if (i < 0 || i > 1)
    o_Notify(O_NFY_FATAL, "o_Vector3D::operator[]: index out of range\n");
  return data[i];
}


inline o_Vector2D o_Vector2D::operator/ (double a) const
{
  return o_Vector2D (data[0] / a,
		     data[1] / a);
}

inline double o_Vector2D::operator* (const o_Vector2D &vec) const
{
  double result = 0.0;
  result += data[0] * vec.data[0];
  result += data[1] * vec.data[1];
  return result;
}

inline double o_Vector2D::Length() const
{
  return sqrt (*this * *this);
}

inline o_Vector2D o_Vector2D::Norm() const
{
  return *this/Length();
}

inline int o_Vector2D::Dimension() const
{
  return 2;
}

inline double o_Vector2D::angle (const o_Vector2D &vec) const
{
  double l = Length();
  double yl = vec.Length();
  
  if ( l==0 || yl==0)
    o_Notify(O_NFY_FATAL,"angle: zero argument\n");
  
  // catch acos: DOMAIN error
  double	arg=  ((*this)*vec/(l*yl));  
  int	isone = (arg > 1.0-1e-38 );
  if(isone) return 0;
  int   isminusone = (arg < -1.0+1e-38);
  if (isminusone)
    return acos(-1.0);
  else
    return acos( arg );
}

inline double o_Vector2D::RelativeAngle (const o_Vector2D &vec) const
{
  double l = Length();
  double yl = vec.Length();
  
  if ( l==0 || yl==0)
    o_Notify(O_NFY_FATAL,"angle: zero argument\n");
  
  double tmp = atan2(vec.GetY(),vec.GetX()) - atan2(this->GetY(),this->GetX());
        
  if(tmp >  M_PI) tmp -= 2.0*M_PI;
  if(tmp < -M_PI) tmp += 2.0*M_PI;
  return tmp;
}

inline o_Vector2D o_Vector2D::operator+ (const o_Vector2D &vec) const
{
  return o_Vector2D (data[0] + vec.data[0],
		     data[1] + vec.data[1]);
}

inline o_Vector2D o_Vector2D::operator- (const o_Vector2D &vec) const
{
  return o_Vector2D (data[0] - vec.data[0],
		     data[1] - vec.data[1]);
}

inline o_Vector2D o_Vector2D::operator* (double a) const
{
  return o_Vector2D (data[0] * a,
		     data[1] * a);
}

inline o_Vector2D o_Vector2D::Perp() const
{
  return o_Vector2D (data[1], -data[0]);
}

inline int o_Vector2D::operator== (const o_Vector2D &vec) const
{
  int	ret;

  if (data[0] != vec.data[0]) ret = 0;
  else if (data[1] != vec.data[1]) ret = 0;
  else ret = 1;
  return ret;
}

inline int o_Vector2D::operator!= (const o_Vector2D &vec) const
{
  return !(*this == vec);
}

inline o_Vector2D &o_Vector2D::operator= (const leda::vector &vec)
{
  if (vec.dim() == 2)
    {
      data[0] = vec[0];
      data[1] = vec[1];
    }
  else
    o_Notify(O_NFY_WARN, "o_Vector3D::operator=(): dimension mismatch\n");
  return *this;
}

inline o_Vector2D &o_Vector2D::operator= (const o_Vector2D &vec)
{
  data[0] = vec.data[0];
  data[1] = vec.data[1];
  return *this;
}

///////////////////////////////////////////////////////////////////////////
//
//Implementation of o_HVector
//
inline o_HVector o_HVector::operator+ (const o_HVector &vec) const
{
  return o_HVector (data[0] + vec.data[0],
                   data[1] + vec.data[1],
                   data[2] + vec.data[2],
                   data[3] + vec.data[3]);
}

inline o_HVector o_HVector::operator- (const o_HVector &vec) const
{
  return o_HVector (data[0] - vec.data[0],
                   data[1] - vec.data[1],
                   data[2] - vec.data[2],
                   data[3] - vec.data[3]);
}

inline double o_HVector::operator* (const o_HVector &vec) const
{
  double result = 0.0;
  result += data[0] * vec.data[0];
  result += data[1] * vec.data[1];
  result += data[2] * vec.data[2];
  result += data[3] * vec.data[3];
  return result;
}

inline o_HVector o_HVector::operator* (double a) const
{
  return o_HVector (data[0] * a,
                   data[1] * a,
                   data[2] * a,
                   data[3] * a);
}

inline o_HVector o_HVector::operator/ (double a) const
{
  return o_HVector (data[0] / a,
                   data[1] / a,
                   data[2] / a,
                   data[3] / a
                  );
}

inline o_HVector &o_HVector::operator += (const o_HVector &vec)
{
  double	*v1 = &data[0];
  double	const *v2 = &(vec.data[0]);
  
  *v1++ += *v2++;
  *v1++ += *v2++;
  *v1++ += *v2++;
  *v1   += *v2;
  return *this;
}

inline o_HVector &o_HVector::operator -= (const o_HVector &vec)
{
  double	*v1 = &data[0];
  double	const *v2 = &(vec.data[0]);
  
  *v1++ -= *v2++;
  *v1++ -= *v2++;
  *v1++ -= *v2++;
  *v1   -= *v2;
  return *this;
}

inline double o_HVector::operator[] (int i) const
{
  if (i < 0 || i > 3)
    o_Notify(O_NFY_FATAL, "o_HVector::operator[]: index out of range\n");
  return data[i];
}

inline double &o_HVector::operator[] (int i)
{
  if (i < 0 || i > 3)
    o_Notify(O_NFY_FATAL, "o_Vector3D::operator[]: index out of range\n");
  return data[i];
}

inline o_HVector &o_HVector::operator= (const o_HVector &vec)
{
  double	*v1 = &data[0];
  double	const *v2 = &(vec.data[0]);
  
  *v1++ = *v2++;
  *v1++ = *v2++;
  *v1++ = *v2++;
  *v1   = *v2;
  return *this;
}

#endif
