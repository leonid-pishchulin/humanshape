#ifndef	_o_HMatrix_incl
#define	_o_HMatrix_incl

//
//
//       o_HMatrix.h    -Header for oGeM-matrix
//
//       Neubearbeitung Juni 1997, Oliver Stahlhut
//

#include <iostream>
#include <fstream>
#include "o_World.h"
#include "smatrix.h"
#include "math.h"

using namespace std;

/*!
  \class o_HMatrix o_HMatrix.h
  \brief Class for homogeneous matrix
*/

class  o_HMatrix {
  
 public:
  //!  creates homogeneous (4x4)-matrix
  o_HMatrix();
  ~o_HMatrix();

  //! clear homogeneous matrix, sets all diagonal elements 1
  void Clear();

  //! inserts matrix m into homogenous matrix at position i,j
  void Insert(int i, int j, const leda::matrix &m);

  //! extracts a (smaller) LEDA-matrix from homogenous matrix
  leda::matrix Sub(int i, int j, int m, int n) const;

  //! returns a 3x3-LEDA-matrix with the rotational component of hom. matrix
  leda::matrix GetRotComponent() const { return Sub(0, 0, 3, 3); }
 
  //! returns value of matrix-element i, j 
  double & operator ()(int i, int j);

  //! returns value of matrix-element i, j 
  double operator ()(int i, int j) const;

  //! sets homogenous matrix to hom. matrix m
  o_HMatrix operator =(const o_HMatrix m);

  //! sets homogenous matrix to LEDA-matrix m, dimensions of m must be 4x4
  o_HMatrix & operator =(const leda::matrix &);

  //! checks for identity
  int operator ==(const o_HMatrix m) const;

  o_HMatrix operator + (const o_HMatrix m) const;
  o_HMatrix operator - (const o_HMatrix m) const;
  void operator += (const o_HMatrix m);
  void operator -= (const o_HMatrix m);

  //! multiplies homogenous matrix with vector vec
  o_Vector3D operator *(const o_Vector3D &vec) const;

  //! multiplies homogenous matrix with point p
  o_Point operator *(const o_Point &point) const;

  //! multiplies all elements of matrix with val
  o_HMatrix Mul(const double &val);

  //! applies homogenous matrix to vector
  void ApplyTransform(o_Vector3D &vec) const;
 
  //! applies homogenous matrix to point
  void ApplyTransform(o_HVector &hvec) const;
  
  //! applies homogenous matrix to point
  void ApplyTransform(o_Point &point) const;

  //! multiplies homogenous matrix with hom. matrix
  o_HMatrix operator *(const o_HMatrix &hmat) const;

  //! scales hom. matrix with vector s_vec
  void Scale(const o_Vector3D &s_vec);

  //! adds translation t_vec to matrix
  void AddTranslation(const o_Vector3D &t_vec);

  //! returns the translation vector of the Matrix. 
  void GetTranslation(o_Vector3D &vec) const;

  //! sets translation component of hom. matrix to t_vec
  void SetTranslation(const o_Vector3D &vec);

  o_HMatrix GaussJordan(double *vfield, int m) const;

  //! returns inverted hom. matrix using Gauss-Jordan elimination
  o_HMatrix Invert(void) const;

  //! method for manually setting the homogenous matrix
  void SetValue (double m00, double m01, double m02, double m03,
		 double m10, double m11, double m12, double m13,
		 double m20, double m21, double m22, double m23,
		 double m30, double m31, double m32, double m33);

  /*! sets the rotation-component of the hom. matrix,
    &vec1 and &vec2 constitute the rotational-axis, angle [radiant]*/
  void AddRotation(const o_Vector3D &vec1, const o_Vector3D &vec2, double angle);

  /*! Sets this o_HMatrix to describe the rotation around 
    an axis with the given direction in a fixed coordinate system.*/
  void SetRotation(const o_Vector3D &axis, double angle=0);

  /*! returns the angle and the axis-direction of the
    rotation described by this o_HMatrix.
    In Order to arrange a unique solution the angle 
    is forced to be within the intervall 0..PI. */
  double GetRotation (o_Vector3D &axis) const;

  /*! Sets this o_HMatrix to describe a rotation around the
    axes of a fixed coordinate system. x: 1st, y: 2nd, z: 3rd */
  void SetRotation (double xangle, double yangle, double zangle);

  /*!  calculates the angles for three subsequent rotations around the axes
    of a fixed coordinate aystem from homogenous matrix. The six angle
    values of the two possible solutions are returned via the references */
  void GetRotation (double &x_1, double &y_1, double &z_1,
		    double &x_2, double &y_2, double &z_2) const;

  /*! uses the routine above. It's possible to set boundaries for the three
    rotational angles (default are max. angles -PI to PI). Only the set with
    all three angles within the boundary-ranges is returned via references */
  void GetLimitedAngles (double &rx, double &ry, double &rz,
			 const double rx_min = -M_PI, const double rx_max = M_PI,
			 const double ry_min = -M_PI, const double ry_max = M_PI,
			 const double rz_min = -M_PI, const double rz_max = M_PI);

  friend std::ostream& operator<<(std::ostream &, const o_HMatrix &);
  
  //! sets the precision and width if matrix is output via cout
  void SetPrecision(short pre, short wid);	
				       
  //! returns the matrix in OpenGL matrix array
  void GetOpenGLMatrix(float _mat[16]) const;
  //! set H-matrix with  OpenGL matrix array
  void SetOpenGLMatrix(const float _mat[16]);

  //! reads a H-Matrix from a file
  int Read(const char* filename);
  //! reads a H-Matrix from a file
  int Read(ifstream &fin);
  //! writes a H-Matrix to m a file
  int Write(const char* filename) const;
  //! writes a H-Matrix to m a file
  int Write(ofstream &fout) const;
  
  //! old routine-call! (for backward-compatibility only ...)
  void rotation(o_Vector3D &vec1, o_Vector3D &vec2, double angle); 

  /*! 
    old routine-call! (for backward-compatibility only ...). 
    sets hom. matrix to rotate points around the origin (world-coordinates);
    x-, y-, z-angle [radiant], mode = 1 (default) assumes a fixed coordinate
    system (rotational order is x/y/z), mode = 0 sets a moving coordinate
    system  (rotational order z/y/x : this means for example that the
    rotation of the z-axis depends on the x- and y-angle), mode = 2 uses
    the coordinate definition of BUNOR (Photogrammetry package, rotational
    order y/x/z)*/
  void rotation(double xangle, double yangle, double zangle, int mode=1);   

  /*!
    old routine-call! (for backward-compatibility only ...).
    sets hom. matrix to rotate points around the centrepoint of a
    body. centrepoint should be a vector pointing to the centrepoint
    of the body (use method Center of o_Body to calculate such a
    vector) */
  void rotation(o_Vector3D &centrepoint, double X_angle=1, double Y_angle=0, double Z_angle=0, int mode=1);

  /*! 
    old routine-call! (for backward-compatibility only ...)
    sets hom. matrix to rotate around a relative coordinate system
    given by rel */
  void rotation(o_LocalCoordinateSystem &Cam, double A_angle=1, double H_angle=0, double V_angle=0, int mode=1);

  /*!
    old routine-call! (for backward-compatibility only ...)
    returns the components of the rotational vector that described
    the matrix using mode as defined in SetRotation() . This is an
    inverse function to SetRotation() . Mode 0: moving cs (xyz),
    Mode 1: static cs (xyz), Mode 2: moving (yxz - BUNOR) */
  void GetAngles(double &rx, double &ry, double &rz, int mode=1);
  
  //! old routine-call! (for backward-compatibility only ...) 
  void clear() { this->Clear(); }
  //! old routine-call! (for backward-compatibility only ...) 
  void translation(const o_Vector3D &t_vec) { AddTranslation(t_vec); }
  //! old routine-call! (for backward-compatibility only ...) 
  void scale(const o_Vector3D &s_vec) { Scale(s_vec); }
  //! old routine-call! (for backward-compatibility only ...) 
  o_HMatrix inv(void) const { return Invert(); }
  
  
  

 protected:

 private:
  double hmat[4][4];
  short precision; //used for output of matrix
  short width; //used for output of matrix

};               

#endif












