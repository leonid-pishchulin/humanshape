
/* SCCS-ID: $Id: smatrix.h 431 2007-09-26 00:26:07Z rhys $ (c) Oliver Grau */


/*******************************************************************************
+
+  LEDA  3.0
+
+
+  matrix.h
+
+
+  Copyright (c) 1992  by  Max-Planck-Institut fuer Informatik
+  Im Stadtwald, 6600 Saarbruecken, FRG     
+  All rights reserved.
+ 
*******************************************************************************/


#ifndef LEDA_MATRIX_H
#define LEDA_MATRIX_H

static const double O_NULL_DET_EPSILON=1e-12;

//------------------------------------------------------------------------------
//  matrices
//------------------------------------------------------------------------------

#ifdef WIN32
#include <stdlib.h>
#endif

#include <iostream>
using namespace std;  

#include "sbasic.h"
#include "svector.h"


namespace leda { 

  class matrix
    {
      vector** v;
      int  d1;
      int  d2;
      
      void     flip_rows(int,int);
      void     check_dimensions(const matrix&) const; 
      double&  elem(int i, int j) const { return v[i]->v[j]; }
      double** triang(const matrix&, int&, bool allow_null_det=false) const;
      
    public:
      
      matrix(int=0, int=0);
      matrix(const matrix&);
      matrix(const vector&);
      matrix(int,int,double**);
      
      matrix& operator=(const matrix&);
      
      ~matrix();
      
      //LEDA_MEMORY(matrix)
      
      //! returns the number of rows
      int     dim1()  const  {  return d1; }
      //! returns the number of cols
      int     dim2()  const  {  return d2; }
      
      //! returns the number of rows
      int     rows()  const  {  return d1; }
      //! returns the number of cols
      int     cols()  const  {  return d2; }
      
      vector& row(int) const;
      vector  col(int i) const;
      matrix  trans() const;
      
      //! allow null determinante, if true: don't exit with o_error_handler
      matrix  inv(bool allow_null_det=false)   const;
      //! allow null determinante, if true: don't exit with o_error_handler
      double  det(bool allow_null_det=false)   const;
      
      matrix solve(const matrix&, bool allow_null_det=false) const;
      vector solve(const vector& b, bool allow_null_det=false) const { 
	return vector(solve(matrix(b), allow_null_det)); 
      }
      
      matrix Sub(int row, int col, int dim1, int dim2);
      void   Insert(int row, int col, matrix& M);
      
      /*! Cholesky LR decomposition  A= LL^T. The result ist the Matrix L, 
	returns true on succeed otherwise false and leave the matrix as it is) */
      bool LRDecomposeCholesky();

      operator vector() const; 
      
      int     operator==(const matrix&)    const;
      int     operator!=(const matrix& x)  const { return !(*this == x); }
      
      vector& operator[](int i)    const { return row(i); }
      double& operator()(int, int);
      double operator() (int, int) const;
      
      //double  operator()(int i, int j) const { return operator()(i,j); };
      
      matrix operator+(const matrix&);
      matrix operator-(const matrix&);
      
      matrix operator*(double);
      matrix operator*(const matrix&);
      vector operator*(const vector& v) { return vector(*this * matrix(v)); }
      
      friend std::ostream& operator<<(std::ostream&, const matrix&);
      friend std::istream& operator>>(std::istream&, matrix&);
      
    };
  
  inline void Print(const matrix& m, std::ostream& out=cout) { out << m; }
  inline void Read(matrix& m, std::istream& in=cin)          { in >> m;  }
  
#ifdef	notdef
  inline int compare(const matrix&, const matrix&) 
    { 
      error_handler(1,"compare not defined for type `matrix`"); 
      return 0;
    }
  
  LEDA_TYPE_PARAMETER(matrix)
#endif
    
    }  /* eof namespace leda */

#endif
