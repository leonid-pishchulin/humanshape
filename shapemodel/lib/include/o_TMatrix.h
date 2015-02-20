#ifndef	_o_TMatrix_incl
#define	_o_TMatrix_incl

#ifdef WIN32
#include <stdlib.h>
#endif
#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <o_TVector.h>

/*!
  \class o_TMatrix o_TMatrix.h
  \brief Template Matrix class
*/

template <class T>
class o_TMatrix {
private:
  int nn;
  int mm;
  T **v;
public:
  o_TMatrix();
  o_TMatrix(int n, int m);			// Zero-based array
  o_TMatrix(int n, int m, const T &a);	//Initialize to constant
  o_TMatrix(int n, int m, const T *a);	// Initialize to array
  o_TMatrix(const o_TMatrix &rhs);		// Copy constructor
  o_TMatrix & operator=(const o_TMatrix &rhs);	//assignment
  typedef T value_type; // make T available externally
  inline T* operator[](const int i);	//subscripting: pointer to row i
  inline const T* operator[](const int i) const;
  inline int nrows() const;
  inline int ncols() const;
  inline unsigned size1() const;
  inline unsigned size2() const;
  void resize(int newn, int newm); // resize (contents not preserved)
  void assign(int newn, int newm, const T &a); // resize and assign a constant value
  ~o_TMatrix();
};

template <class T>
o_TMatrix<T>::o_TMatrix() : nn(0), mm(0), v(NULL) {}

template <class T>
o_TMatrix<T>::o_TMatrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
  int i, nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1;i<n;i++) v[i] = v[i-1] + m;
  // int j;
  // for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = T(0);
}

template <class T>
o_TMatrix<T>::o_TMatrix(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
  int i,j,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< n; i++) v[i] = v[i-1] + m;
  for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

template <class T>
o_TMatrix<T>::o_TMatrix(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
  int i,j,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< n; i++) v[i] = v[i-1] + m;
  for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

template <class T>
o_TMatrix<T>::o_TMatrix(const o_TMatrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL)
{
  int i,j,nel=mm*nn;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
o_TMatrix<T> & o_TMatrix<T>::operator=(const o_TMatrix<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
  if (this != &rhs) {
    int i,j,nel;
    if (nn != rhs.nn || mm != rhs.mm) {
      if (v != NULL) {
        delete[] (v[0]);
        delete[] (v);
      }
      nn=rhs.nn;
      mm=rhs.mm;
      v = nn>0 ? new T*[nn] : NULL;
      nel = mm*nn;
      if (v) v[0] = nel>0 ? new T[nel] : NULL;
      for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
    }
    for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
  }
  return *this;
}

template <class T>
inline T* o_TMatrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
  if (i<0 || i>=nn) {
    throw("o_TMatrix subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
inline const T* o_TMatrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
  if (i<0 || i>=nn) {
    throw("o_TMatrix subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
inline int o_TMatrix<T>::nrows() const
{
  return nn;
}

template <class T>
inline int o_TMatrix<T>::ncols() const
{
  return mm;
}

template <class T>
inline unsigned o_TMatrix<T>::size1() const
{
  return nn;
}

template <class T>
inline unsigned o_TMatrix<T>::size2() const
{
  return mm;
}

template <class T>
void o_TMatrix<T>::resize(int newn, int newm)
{
  int i,nel;
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[] (v[0]);
      delete[] (v);
    }
    nn = newn;
    mm = newm;
    v = nn>0 ? new T*[nn] : NULL;
    nel = mm*nn;
    if (v) v[0] = nel>0 ? new T[nel] : NULL;
    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  }
}

template <class T>
void o_TMatrix<T>::assign(int newn, int newm, const T& a)
{
  int i,j,nel;
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[] (v[0]);
      delete[] (v);
    }
    nn = newn;
    mm = newm;
    v = nn>0 ? new T*[nn] : NULL;
    nel = mm*nn;
    if (v) v[0] = nel>0 ? new T[nel] : NULL;
    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  }
  for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
}

template <class T>
o_TMatrix<T>::~o_TMatrix()
{
  if (v != NULL) {
    delete[] (v[0]);
    delete[] (v);
  }
}


#endif
