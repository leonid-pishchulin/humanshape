//
// o_NurbCurve.cc
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



#ifndef	_o_NurbsCurve
#define	_o_NurbsCurve

#include "o_PlanarFace.h"
#include <iostream> 
#include <vector>
#include "o_HMatrix.h"

/*! 
  \class o_NurbsCurve o_NurbsCurve.h
  \brief A class to represent a NURBS curve

  This class is used to represent and manipulate NURBS curves, which 
  can have any degree and have any number of control points.
  Most algorithms are taken from 'The NURBS Book' by Piegl and Tiller or from the
  NURBS++ library http://libnurbs.sourceforge.net/
  
*/

class o_NurbsCurve:public o_MetaClass {

public:
  const char *GetClassName() const;
  
  //!default constructor
  o_NurbsCurve();
  
  //! create a NURBS curve without root node spezified
  o_NurbsCurve(unsigned knot_count_U, unsigned degreeU); 
  
  //! assign NURBS curve with name nam to body b
  o_NurbsCurve(unsigned knot_count_U, unsigned degreeU, o_Body &b, char *nam ) ; 
  
  //! copy constructor
  o_NurbsCurve(const o_NurbsCurve& in);
  
  //! = operator
  o_NurbsCurve& operator=(const o_NurbsCurve& in);
  
  ~o_NurbsCurve( );

  //! get size of U knot vector
  unsigned GetUknotCount() const {return m_U.size();}
  //! get degree of the curve in U
  unsigned GetUDegree() const {return m_degU;}
  //! get the order of the curve which is degree+1
  unsigned GetUOrder() const {return m_degU+1;}
  //! get u stride
  unsigned GetUStride() const {return m_u_stride;}

  //! set a control point
  void SetCP(unsigned u, unsigned z, float value) const { *(m_P+u*m_u_stride+z)=value; }
   //! set a control point
  void SetCP(unsigned u, const o_Vector3D& vec) const { for(unsigned i=0; i<3; i++) *(m_P+u*m_u_stride+i)=float(vec[i]); *(m_P+u*m_u_stride+3)=1.0;}
  //! set a control point
  void SetCPH(unsigned u, const o_HVector& vec) const { for(unsigned i=0; i<4; i++) *(m_P+u*m_u_stride+i)=float(vec[i]); }
  //! get a control point
  float GetCP(unsigned u, unsigned z) const { return *(m_P+u*m_u_stride+z); }
  //! get a control point vector
  float* GetCPv(unsigned u) const { return m_P+u*m_u_stride; }
  //! get a control point o_Vector3D
  const o_Vector3D GetCP(unsigned u) const { return o_Vector3D(*(m_P+u*m_u_stride), *(m_P+u*m_u_stride+1), *(m_P+u*m_u_stride+2)); } 
  //! get a control point o_HVector
  const o_HVector GetCPH(unsigned u) const { return o_HVector(*(m_P+u*m_u_stride), *(m_P+u*m_u_stride+1), *(m_P+u*m_u_stride+2),*(m_P+u*m_u_stride+3)); }
  
  //! get number of control points
  unsigned GetNumCP() const { return m_P_count; } 
  
   //! get control points data pointer for opengl rendering
  float * GetOpenGL_CP(){return m_P;}
  //! get U knot vector data pointer for opengl rendering
  float * GetOpenGL_Uknot(){return &m_U[0];}

    //! get U knot vector 
  const std::vector< float >& GetUKnot() const {return m_U;} 
  
  //! connect the object to a body.
  void ConnectRootNode(o_Body &b);
  
   /*! \brief global curve interpolation with points in 3D
   
     \param Q the points in 3D, 
     \param deg the degree of the interpolation 
   */
  void globalInterp(const std::vector< o_Vector3D >& Q, int deg);
  
  /*! \brief global curve interpolation with points in 3D
  
    \param Q the points in 3D, 
    \param ub parametric values
    \param deg the degree of the interpolation 
   */
  void globalInterp(const std::vector< o_Vector3D >& Q, const  std::vector < float >& ub, int deg);
  
  /*!
  \brief global curve interpolation with 4D points, 
  a knot vector defined and the parametric value
  vector defined.

  Global curve interpolation with 4D points, a knot vector
  defined and the parametric value vector defined.

  \param Q  the 3D points to interpolate
  \param ub  the parametric values vector
  \param Uc  the knot vector to set the curve to
  \param d  the degree of the interpolation

   */
  void globalInterpH(const std::vector< o_HVector >& Q, const std::vector < float >& ub, const std::vector < float >& Uc, int d);
  
  
  //! perform chord length parameterization from a vector of points Q
  float chordLengthParam(const std::vector< o_Vector3D >& Q, std::vector < float >& ub);
  
  /*! \brief generates a knot vector using the averaging technique
  
      \param uk are the knot coefficients
      \param deg is the degree of the curve associated with the knot vector
      \param U is the average knot vector
   */
  void knotAveraging(const std::vector < float >& uk, int deg, std::vector < float >& U);
  
  /*!
  \brief determines the knot span index
    
   Determines the knot span for which their exists non-zero basis 
   functions. The span is the index for which the parameter 
  \a u is valid in the [u_k,u_{k+1}] range.
  
  \param   u  the parametric value
  \return the span index at \a u.
   */
  int findSpan(float u) const;
    
  
  /*!
  \brief Computes the basis function of the curve
    
  Computes the \a i basis function of degree \a p  of the curve at 
  parameter \a u. 
  \param u  the parametric variable
  \param i  specifies which basis function to compute
  \param p  the degree to which the basis function is computed
  \return the value of  \a N_{ip}(u)
 */
  float basisFun(float u, int i, int p=-1) const;
  
  /*!
  \brief computes the non-zero basis functions of the curve
    
   Computes the non-zero basis functions and puts the result 
   into \a N. \a N has a size of deg+1. To relate \a N to the basis 
   functions, Basis[span -deg +i] = N[i] for i=0...deg.
   see chapter 2 'B-Spline Basis Functions' of 'The NURBS Book' by Piegl and Tiller for details
   \param u  the parametric value
   \param i  the non-zero span of the basis functions
   \param N  the non-zero basis functions
  */
  void basisFuns(float u, int i, std::vector <float>& N) const;
  
  /*!
  \brief A least squares curve approximation
    
  This routines generates a curve that approrimates the points in the 
  Least square sense, you can find more details in the LaTeX version.
  
  For more details, see section 9.4.1 on page 491 of the NURBS 
  book.
  
  \param Q  the vector of 3D points
  \param degC  the degree of the curve
  \param n  the number of control points in the new curve.
  
  \return 1 if succesfull, 0 if the number of points to approximate 
          the curve with is too big compared to the number of points.
  
  */
  int leastSquares(const std::vector< o_Vector3D >& Q, int degC, int n);
  
  /*!
  \brief  A least squares curve approximation
    
  This routines generates a curve that approrimates the points in the 
  least square sense, you can find more details in the LaTeX version.
  
  For more details, see section 9.4.1 on page 491 of the NURBS 
  book.
  
  \param Q  the vector of 3D points
  \param degC  the degree of the curve
  \param n  the number of control points in the new curve
  \param ub  the knot coefficients
  
  \return 1 if succesfull, 0 it the number of points to approximate 
          the curve with is too big compared to the number of points.
  */
  int leastSquares(const std::vector< o_Vector3D >& Q, int degC, int n, std::vector < float >& ub);
  
  /*!
  \brief  A least squares curve approximation
    
  This routines generates a curve that approrimates the points in the 
  least square sense, you can find more details in the LaTeX version.
  
  For more details, see section 9.4.1 on page 491 of the NURBS 
  book.
  
  \param Q  the vector of 3D points
  \param degC  the degree of the curve
  \param n  the number of control points in the new curve
  \param ub  the knot coefficients
  \param knot  the knot vector to use for the curve
  
  \return 1 if succesfull, 0 it the number of points to approximate 
          the curve with is too big compared to the number of points.
  */
  int leastSquares(const std::vector< o_Vector3D >& Q, int degC, int n, const std::vector < float >& ub, const std::vector < float >& knot);
  
  /*!
  \brief Approximation of a curve bounded to a certain error
    
  It is a type II approximation: it starts with a lot of control
  points then tries to eliminate as much as it can as long as
  the curve stays within a certain error bound.
  
  The method uses least squares fitting along with knot
  removal techniques. It is the algorithm A9.10 on p 431 of 
  the NURBS book.
  
  \param  Q  the points to approximate
  \param degree  the degree of the approximation curve
  \param E  the maximum error allowed
  */
  void globalApproxErrBnd(std::vector< o_Vector3D >& Q, int degC, float E);
  
  /*!
  \brief Approximation of a curve bounded to a certain error
    
  It is a type II approximation: it starts with a lot of control
  points then tries to eliminate as much as it can as long as
  the curve stays within a certain error bound.
  
  The method uses least squares fitting along with knot
  removal techniques. It is the algorithm A9.10 on p 431 of 
  the NURBS book.
  
  \param  Q  the points to approximate
  \param ub  the vector of parameters where the points are located
  \param degree  the degree of the approximation curve
  \param E  the maximum error allowed
  */
  void globalApproxErrBnd(std::vector< o_Vector3D >& Q, std::vector < float >& ub, int degC, float E);
  
  /*!
  \brief Remove knots from a curve without exceeding an error bound

  For more information about the algorithm, see A9.9 on p429 of the NURB book.

  \param ub  the knot coefficients
  \param ek  the error after removing
  \param E  the maximum error allowed
   */
  void removeKnotsBound(const std::vector < float >& ub, std::vector < float >& ek, float E);
  
  /*!
  \brief Get the knot removal error bound for an internal knot 
    
  Get the knot removal error bound for an internal knot r 
  (non-rational). For more information on the algorithm, see 
  A9.8 from the Nurbs book on page 428.
  
  
  \param curve  a NURBS curve
  \param r  the index of the internal knot to check
  \param s  the multiplicity of that knot
  
  \return The maximum distance between the new curve and the old one
  */
  float getRemovalBnd(int r, int s ) const;
  
  /*!
  \brief Removes an internal knot from a curve.
  This is A5.8 on p185 from the NURB book modified to not check for 
  tolerance before removing the knot.
  
  \param r  the knot to remove
  \param s  the multiplicity of the knot
  \param num  the number of times to try to remove the knot
  */
  void removeKnot(int  r, int s, int num);
  
  /*!
  \brief degree elevate a curve a number of times
    
  For more information, see A5.9 on p 206 of the NURBS book
  
  \param t  the number of times to increase the degree of the curve
  */
  void degreeElevate(int t);
  
  /*!
  \brief project a point onto the curve
    
  It finds the closest point in the curve to a point $p$.
  For more information, see chapter 6.1 'Point Inversion and Projection for Curves and Surfaces' of the 
  NURBS book
  
  \param p  the point \a p being projected
  \param guess  initial guess
  \param u  the parametric value for which \a C(u) is closest to \a p.
    \param r  the point at \a C(u)
    \param e1  the minimal error
  \param e2  the maximal error
  \param maxTry  the maximum number of times to try before returning from the function
  
  */
  bool projectTo(const o_Vector3D& p, float guess, float& u, o_Vector3D& r, float e1=0.001, float e2=0.001, int maxTry=100) const;
  
  /*!
  \brief Evaluates the curve in 3D at parameter \a u

  For more details on the algorithm, see A4.1 on page 124 of 
  the Nurbs Book.
  
  \param u the parametric value at which the curve is evaluated
  
  \return the 3D point at \a C(u)
  */
  o_Vector3D pointAt(float u) const;
  
  /*!
  \brief Computes the derivative of degree \a d of the 
         curve at parameter \a u in the homonegeous domain
  
  For more information on the algorithm used, see A3.2 on p 93 
  of the NurbsBook.
  
  \param u  the parametric value to evaluate at
  \param d  the degree of the derivative
  \param ders a vector containing the derivatives of the curve at \a u.
    
  */
  void deriveAtH(float u,int d, std::vector < o_HVector >& ders) const;
  
  /*!
  \brief Compute the derivatives functions at \a u of the NURBS curve
    
  For information on the algorithm, see A2.3 on p72 of the NURBS 
  book.
  
  The result is stored in the ders matrix, where ders is of 
  size \a (n+1,deg+1) and the derivative 
  N'_i(u) = ders(1,i=span-deg+j) where j=0...deg+1.
  
  \param  n   the degree of the derivation
  \param  u   the parametric value
  \param span  the span for the basis functions
  \param ders  A matrix containing the derivatives of the curve.
  
  */
  void dersBasisFuns(int n, float u, int span, leda::matrix& ders) const;
  
  /*!
  \brief Computes the derivative at the parameter \a u
    
    \param u  the parameter at which the derivative is computed
  \param d  the degree of derivation
  \param ders  the vector containing the derivatives of the point at \a u.
    
  */
  void deriveAt(float u, int d, std::vector< o_Vector3D >& ders) const;
    
  /*!
  \brief Merges the knot vector of a curve with another knot vector
  
  Will merge the Knot vector U with the one from the curve  
  and it will refine the curve appropriately.
  
  \param Um  the knot vector to merge with
   */
  void mergeKnotVector(const  std::vector <float> &Um);
  
  /*!
  \brief Refine the curve knot vector

  For more information, see A5.4 on page 164 of the NURBS book

  \param X the new knots to insert in the knot vector
   */
  void refineKnotVector(const std::vector <float>& X);
  
  //! transform the curve 
  void ApplyTransform(o_HMatrix &hmat);
  
  /*!
  \brief Reads a NurbsCurve from a file.
  \param filename the filename to read the curve from 
  \return 0 if an error occurs, 1 otherwise

   */
  int read(const char* filename);
  
  /*!
  \brief reads a NurbsCurve from a file
  \param fin  an input file stream
  \return 0 if an error occurs, 1 otherwise
   */
  int read(ifstream &fin);
  
  /*!
  \brief Writes a NurbsCurve<T,N> to a file.
  
  \param filename  the filename to write to.
  \return 0 if an error occurs, 1 otherwise
   */
  int write(const char* filename) const;
  
  /*!
  \brief Writes a NurbsCurve to an output stream.

  \param fout  the output stream
  \return 0 if an error occurs, 1 otherwise

   */
  int write(ofstream &fout) const;
  
  
  friend std::ostream& operator<< (std::ostream &, const o_NurbsCurve &);
  
private:
  // initalize members with correct size
  void InitMem(unsigned knot_count_U, unsigned degreeU);
  // resize members with correct size
  void ResizeMem(unsigned knot_count_U, unsigned degreeU);
  // resize members with correct size
  void ResizeMem();
  // copy data from one curve to the other
  void CopyMem(const o_NurbsCurve& in);
  // setup a matrix containing binomial coefficients
  void binomialCoef(leda::matrix& Bin) const;
  
protected:
  std::vector< float > m_U ; //!< the U knot vector
  float* m_P ; //!< The matrix of rational control points (4 floats per vertex)
  unsigned m_P_count;
  unsigned m_degU ; //!< the degree of the curve in U
  unsigned m_u_stride;
};

/*!
    \brief an array of NurbsCurve
    
    This class represents an array of NurbsCurve.
 */
class o_NurbsCurveArray {
  public:
    //!returns the size
    unsigned size() const { return sze ; }
    //! Contructor
    o_NurbsCurveArray() { C = 0 ; sze = 0 ; rsize = 0 ;}
    //!Constructor from a pointer to an array of curves
    o_NurbsCurveArray(o_NurbsCurve* Ca, int size) ;
    //! Destructor
    virtual ~o_NurbsCurveArray(){ if(C){ for(int i=0;i<rsize;i++) delete C[i];  delete []C ; }}
      
    virtual o_NurbsCurve& operator[](int i) { return *(C[i]) ; }
    virtual o_NurbsCurve operator[](int i) const { return *(C[i]) ; }
      
    virtual void resize(int s) ;
    void init(o_NurbsCurve* Ca, int size) ;

    /*!
    \brief  Generate compatible curves from an array of curves
    
    This routine will force the same degree to all the curves in 
    the array and it will ensure that they have the same knot 
    vector.
     */
    void generateCompatibleCurves();
    
  protected:
    o_NurbsCurve& curve(int i) { return *(C[i]) ; }
    o_NurbsCurve curve(int i) const { return *(C[i]) ; }
      
    unsigned sze ; // the number of NURBS curves in the array
    int rsize ; // the number of space allocated for the array
    o_NurbsCurve** C ; // An array of pointers to NURBS curves
  private:
    void knotUnion(const std::vector < float >& Ua, const std::vector < float >& Ub, std::vector < float >& U);

};


#endif
