//
// o_NurbSurface.cc
//
// Copyright (C) 2006 University of Adelaide
//
// Permission is granted to any individual or institution to use, copy, modify,
// and distribute this software, provided that this complete copyright and
// permission notice is maintained, intact, in all copies and supporting
// documentation.
//
// University of Adelaide provides this software "as is" without
// express or implied warranty.
//
//
// History:
//
// Created: T. Thormaehlen  September 2006      Initial design and implementation.
// Updated:
//



#ifndef	_o_NurbsSurface
#define	_o_NurbsSurface

#include "o_PlanarFace.h"
#include <iostream> 
#include <vector>
#include "o_NurbsCurve.h"

/*! 
  \class o_NurbsSurface o_NurbsSurface.h
  \brief A class to represent a NURBS surface

  Most algorithms are taken from 'The NURBS Book' by Piegl and Tiller or from the
  NURBS++ library http://libnurbs.sourceforge.net/
*/

class o_NurbsSurface : public o_MetaSurface {

public:
  const char *GetClassName() const;
  //! create a NURBS surface without root node spezified
  o_NurbsSurface(unsigned knot_count_U, unsigned knot_count_V, unsigned degreeU, unsigned degreeV); 
  
  //! assign NURBS surface with name nam to body b
  o_NurbsSurface(unsigned knot_count_U, unsigned knot_count_V, unsigned degreeU, unsigned degreeV, o_Body &b, char *nam ) ; 
    
  //! copy constructor
  o_NurbsSurface(const o_NurbsSurface& in);
  
  //! = operator
  o_NurbsSurface& operator=(const o_NurbsSurface& in);
  
  ~o_NurbsSurface( );
  
  std::vector< o_NurbsCurve* >& GetTrimmingCurves() {return trimming_curves;}
  
  //! get size of U knot vector
  unsigned GetUknotCount() const {return m_U.size();} 
    //! get size of V knot vector
  unsigned GetVknotCount() const {return m_V.size();} 
  //! get degree of the surface in U
  unsigned GetUDegree() const {return m_degU;}
  //! get degree of the surface in V
  unsigned GetVDegree() const {return m_degV;}
  //! get the order of the curve in U, which is degree+1
  unsigned GetUOrder() const {return m_degU+1;}
  //! get the order of the curve in V, which is degree+1
  unsigned GetVOrder() const {return m_degV+1;}
  //! get strinde in U
  unsigned GetUStride() const {return m_u_stride;}
  //! get strinde in V
  unsigned GetVStride() const {return m_v_stride;}

   //! set a control point
  void SetCP(unsigned u, unsigned v, unsigned z, float value) { *(m_P+u*m_u_stride+v*m_v_stride+z)=value; }
   //! set a control point
  void SetCP(unsigned u, unsigned v, const o_Vector3D& vec) const { for(unsigned i=0; i<3; i++) *(m_P+u*m_u_stride+v*m_v_stride+i)=float(vec[i]); *(m_P+u*m_u_stride+v*m_v_stride+3) = 1.0; }
  //! set a control point
  void SetCPH(unsigned u, unsigned v, const o_HVector& vec) const { for(unsigned i=0; i<4; i++) *(m_P+u*m_u_stride+v*m_v_stride+i)=float(vec[i]); }
  //! get a control point
  float GetCP(unsigned u, unsigned v, unsigned z) const { return *(m_P+u*m_u_stride+v*m_v_stride+z); }
  //! get a control point vector
  float* GetCPv(unsigned u, unsigned v) const { return m_P+u*m_u_stride+v*m_v_stride; }
  //! get a control point o_Vector3D
  const o_Vector3D GetCP(unsigned u, unsigned v) const { float* off=m_P+u*m_u_stride+v*m_v_stride; return o_Vector3D(*(off), *(off+1), *(off+2)); } 
  //! get a control point o_HVector
  const o_HVector GetCPH(unsigned u, unsigned v) const { float* off=m_P+u*m_u_stride+v*m_v_stride; return o_HVector(*(off), *(off+1), *(off+2),*(off+3)); }
  
  
  //! get number of control points in U
  unsigned GetUNumCP(){ return m_P_u_count; } 
  //! get number of control points
  unsigned GetVNumCP(){ return m_P_v_count; } 
  
   //! get control points data pointer for opengl rendering
  float * GetOpenGL_CP(){return m_P;}
  //! get U knot vector data pointer for opengl rendering
  float * GetOpenGL_Uknot(){return &m_U[0];}
  //! get V knot vector data pointer for opengl rendering
  float * GetOpenGL_Vknot(){return &m_V[0];}
  
  //! get U knot vector 
  const std::vector< float >& GetUKnot() const {return m_U;}
  //! get V knot vector 
  const std::vector< float >& GetVKnot() const {return m_V;}
  //! set U knot vector 
  bool SetUKnot(const std::vector< float >& u);
  //! set V knot vector 
  bool SetVKnot(const std::vector< float >& v);
  
  //! connect the object to a body.
  void ConnectRootNode(o_Body &b);
  
  /*! 
    \brief Generates a surface using global interpolation

    \param Q  a matrix of 3D points
    \param pU  the degree of interpolation in the U direction
    \param pV  the degree of interpolation in the V direction
  */
  void globalInterp(const std::vector < std::vector < o_Vector3D > >& Q, int pU, int pV);
  
  /*! 
  \brief Computes the parameters for global surface interpolation

  Computes the parameters for global surface interpolation. 
  For more information, see A9.3 on p377 on the NURBS book.

  \param Q   the matrix of 3D points
  \param uk  the knot coefficients in the U direction
  \param vl  the knot coefficients in the V direction            

  \return 0 if an error occurs, 1 otherwise
   */
  int surfMeshParams(const std::vector < std::vector < o_Vector3D > >& Q, std::vector <float >& uk, std::vector <float >& vl);
  
  /*! 
  \brief generates a surface using global least squares approximation.

  \param   Q  a matrix of 3D points
  \param  pU  the degree of interpolation in the U direction
  \param  pV  the degree of interpolation in the V direction
  \param  nU  the number of points in the U direction
  \param  nV  the number of poitns in the V direction
  */
  void leastSquares(const std::vector < std::vector < o_Vector3D > >& Q,  int pU, int pV, int nU, int nV);
  
  /*! 
  \brief Generates a surface using global approximation

  \param  Q  a matrix of 3D points
  \param pU  the degree of interpolation in the U direction
  \param pV  the degree of interpolation in the V direction
  \param  S  the interpolated surface
   */
  void globalApproxErrBnd(const std::vector < std::vector < o_Vector3D > >& Q, int pU, int pV, double error);
  
  /*! 
  \brief Returns the point on the surface at \a u,v

  Returns the point on the surface at \a u,v

  \param u  the u parametric value
  \param v  the v parametric value

  \return The homogenous point at \a u,v
   */
  o_HVector hpointAt(float u, float v) const;
  
    /*! 
  \brief Returns the point on the surface at \a u,v

  Returns the point on the surface at \a u,v

  \param u  the u parametric value
  \param v  the v parametric value

  \return The 3D point at \a u,v
     */
  o_Vector3D pointAt(float u, float v) const;
  
  
  /*! 
  \brief finds the span in the U and V direction

  Finds the span in the U and V direction. The spanU is the index
  \a k for which the parameter \a u is valid in the \a [u_k,u_{k+1}]
  range. The spanV is the index \a k for which the parameter \a v is 
  valid in the \a [v_k,v_{k+1}] range.

  \param u  find the U span for this parametric value 
  \param v  find the V span for this parametric value
  \param spanU  the U span
  \param spanV  the V span
   */
  void findSpan(float u, float v, int& spanU, int& spanV) const;
  
  /*! 
  \brief finds the span in the U direction
  
  Finds the span in the U direction. The span is the index
  \a k for which the parameter \a u is valid in the \a [u_k,u_{k+1}]
  range.

  \param u --> find the span for this parametric value

  \return the span for \a u
   */
  int findSpanU(float u) const;
  
  /*! 
  \brief finds the span in the V direction

  Finds the span in the V direction. The span is the index
  \a k for which the parameter \a v is valid in the \a [v_k,v_{k+1}]
  range. 

  \param v  find the span for this parametric value    

  \return the span for \a v
   */
  int findSpanV(float v) const;
  
  /*! 
  \brief Find the non-zero basis functions in the U and V direction

  \param   u  the u parametric value
  \param   v  the v parametric value 
  \param spanU  the span of u
  \param spanV  the span of v
  \param    Nu  the vector containing the U non-zero basis functions
  \param    Nv  the vector containing the V non-zero basis functions
   */
  void basisFuns(float u, float v, int spanU, int spanV, std::vector < float >& Nu, std::vector < float >&Nv) const;
  
  /*! 
  \brief Finds the non-zero basis function in the U direction
        
  \param   u  the u parametric value
  \param span  the span of u
  \param    N  the vector containing the basis functions
   */
  void basisFunsU(float u, int span, std::vector < float >& N) const;
  
  /*! 
  \brief Finds the non-zero basis function in the V direction

  \param    v  the v parametric value
  \param span  the span of v 
  \param    N  the vector containing the basis functions
   */
  void basisFunsV(float v, int span, std::vector < float >& N) const;
  
  /*! 
  \brief Computes the point and the derivatives of degree 
  \a d and below at \a (u,v)

  Computes the matrix of derivatives at \a u,v . 
  The value of skl(k,l) represents the 
  derivative of the surface \a S(u,v) with respect to 
  \a u, \a k times and to \a v, \a l times.

  \param   u  the u parametric value
  \param   v  the v parametric value
  \param   d  the derivative is computed up to and including to  this value
  \param skl  the matrix containing the derivatives
   */
  void deriveAt(float u, float v, int d, std::vector < std::vector < o_Vector3D > > &skl) const;
  
  /*! 
  \brief computes the point and the derivatives of degree 
  \a d and below at \a (u,v)

  Computes the matrix of derivatives at \a u,v . 
  The value of skl(k,l) represents the 
  derivative of the surface \a S(u,v) with respect to 
  \a u \a k times and to \a v \a l times.

  \param  u  the u parametric value
  \param  v  the v parametric value
  \param  d  the derivative is computed up to and including to this value
  \param skl  the matrix containing the derivatives
   */
  void deriveAtH(float u, float v, int d, std::vector < std::vector < o_HVector > > &skl) const;
  
  /*!
  \brief Compute the derivatives functions 

  For information on the algorithm, see A2.3 on p72 of 
  the NURBS book. The result is stored in the ders matrix, where 
  ders is of size (n+1,deg+1) and the derivative 
  C'(u) = ders(1,span-deg+j) where j=0...deg+1.

  \param n  the degree of the derivation
  \param u  the parametric value
  \param  span  the span for the basis functions
  \param deg  the degree of the curve
  \param U  the knot vector on which the Basis functions must be computed.
  \param ders  A matrix containing the derivatives of the curve.
   */
  void dersBasisFuns(int n, float u, int span,  int deg, const std::vector < float > &U, leda::matrix& ders) const; 
  
    /*!
  \brief project a point onto the surface
    
  It finds the closest point in the surface to a point $p$.
  For more information, see chapter 6.1 'Point Inversion and Projection for Curves and Surfaces' of the 
  NURBS book
  
  \param p  the point \a p being projected
  \param guess_u  initial guess for u
  \param guess_v  initial guess for v
  \param u  the parametric value for which \a S(u,v) is closest to \a p.
  \param V  the parametric value for which \a S(u,v) is closest to \a p.
  \param proj  the point on \a S(u,v) )
  \param e1  the maximal error threshold \a e1
  \param e2  the maximal error threshold \a e2
  \param maxTry  the maximum number of times to try before returning from the function
  
     */
  bool projectTo(const o_Vector3D& P, float guess_u,  float guess_v, float& u, float& v, o_Vector3D& proj, float e1=0.001, float e2=0.001, int maxTry=100) const;
  
  
  /*!
  \brief finds the closest point on the surface sampling only at the knots positions
    
  It finds the closest point in the surface to a point $p$ by sampling the knots positions.
  The result can be used to initialize the projectTo() function
  
  \param p  the point \a p being projected
  \param u  the parametric value for which \a S(u,v) is closest to \a p.
  \param V  the parametric value for which \a S(u,v) is closest to \a p.
  \param proj the point on \a S(u,v) )
  \param oversample oversample between the knots positions
  */
  void sampleKnotsForClosestPoint(const o_Vector3D& P, o_Vector3D& proj, float& u, float& v, unsigned oversample=0) const;
  
  /*! 
  \brief Computes the normal of the surface at \a (u,v)

  \param  u  the u parametric value
  \param  v  the v parametric value

  \return the normal at \a (u,v) .
   */
  o_Vector3D normal(float u, float v) const;
  
  /*! 
  \brief Computes the intersection of a ray with the surface (Ray Tracing)

  \param o origin of ray
  \param d unit direction 
  \param guess_u  initial guess for u
  \param guess_v  initial guess for v
  \param u  the parametric value for the intersection with \a S(u,v)
  \param V  the parametric value for the intersection with \a S(u,v) 
  \param proj  the point on \a S(u,v) )
  \param e1  the maximal error threshold
  \param e2  threshold to determine singular Jacobian matrices
  \param maxTry  the maximum number of times to try before returning from the function
  */
  bool traceRay(const o_Vector3D& o, const o_Vector3D& d, float guess_u,  float guess_v, float& u, float& v, o_Vector3D& proj, float e1=0.00001, float e2=1e-11, int maxTry=100) const;
  
  /*!
  \brief finds an initial point for ray tracing by sampling at the knots positions
    
  It finds an initial point for ray tracing by sampling the knots positions.
  The result can be used to initialize the traceRay() function
  
  \param o origin of ray
  \param d unit direction 
  \param u  the parametric value for the initial point
  \param V  the parametric value forthe initial point
  \param proj  the initial point on \a S(u,v) )
   */
  void sampleKnotsForRayTracing(const o_Vector3D& o, const o_Vector3D& d, o_Vector3D& proj, float& u, float& v) const;
  
  /*!
  \brief Generates a surface of revolution

  Generates a surface of revolution of a profile around an
  arbitrary axis (specified by a starting point S and a 
  tangent T). 
  
  \param profile  the curves to rotate around the axis
  \param S  a point on the axis
  \param T  the tangent vector of the rotation axis

   */
  void makeFromRevolution(const o_NurbsCurve& profile, const o_Vector3D& S, const o_Vector3D& Tvec);
  
  /*! 
  \brief Generates a NURBS surface from skinning

  The NURBS surface is generated from skinning. The skinning is performed 
  in the V direction.

  \param   ca  an array of NURBS curves
  \param  degV  the degree to skin in the V direction

  \return 0 if an error occurs, 1 otherwise
  */
  int  skinV(o_NurbsCurveArray& ca, int dV);
  
  /*! 
  \brief Generates a NURBS surface from skinning

  The NURBS surface is generates from skinning. The skinning is performed 
  in the U direction.

  \param    ca  an array of NURBS curves
  \param   degU  the degree to skin in the Udirection
  
  \return 0 if an error occurs, 1 otherwise
  */
  int skinU(o_NurbsCurveArray& ca, int degU);
  
  //! transform the curve 
  void ApplyTransform(o_HMatrix &hmat);
  
  friend std::ostream& operator<< (std::ostream &, const o_NurbsSurface &);
  
private:
    // initalize members with correct size
  void InitMem(unsigned knot_count_U, unsigned knot_count_V, unsigned degreeU, unsigned degreeV);
  // resize members with correct size
  void ResizeMem(unsigned knot_count_U, unsigned knot_count_V, unsigned degreeU, unsigned degreeV);
  // resize members with correct size
  void ResizeMem();
  // copy data from one surface to the other
  void CopyMem(const o_NurbsSurface& in);
  // setup a matrix containing binomial coefficients
  void binomialCoef(leda::matrix& Bin) const;
  
protected:
  std::vector< float > m_U ; //!< the U knot vector
  std::vector< float > m_V ; //!< the V knot vector
  float* m_P ; //!< The matrix of rational control points (4 floats per vertex)
  unsigned m_P_u_count;
  unsigned m_P_v_count;
  unsigned m_degU ; //!< the degree of the surface in U
  unsigned m_degV ; //!< the degree of the surface in V
  unsigned m_u_stride;
  unsigned m_v_stride;
  std::vector< o_NurbsCurve* > trimming_curves; //!< a list of trimming curves
};


#endif
