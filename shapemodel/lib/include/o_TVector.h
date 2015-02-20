#ifndef	_o_TVector_incl
#define	_o_TVector_incl

#ifdef WIN32
#include <stdlib.h>
#endif
#include <iostream>
#include <fstream>


/*!
  \class o_TMatrix o_TMatrix.h
  \brief Template Matrix class
*/

template <class T>
class o_TVector {
private:
  int nn;	// size of array. upper index is nn-1
  T *v;
public:
  o_TVector();
  explicit o_TVector(int n);	// Zero-based array
  o_TVector(int n, const T &a);	//initialize to constant value
  o_TVector(int n, const T *a);	// Initialize to array
  o_TVector(const o_TVector &rhs);	// Copy constructor
  o_TVector & operator=(const o_TVector &rhs);	//assignment
  typedef T value_type; // make T available externally
  inline T & operator[](const int i);	//i'th element
  inline const T & operator[](const int i) const;
  inline int size() const;
  void resize(int newn); // resize (contents not preserved)
  void assign(int newn, const T &a); // resize and assign a constant value
  ~o_TVector();
};

// o_TVector definitions

template <class T>
o_TVector<T>::o_TVector() : nn(0), v(NULL) {}

template <class T>
o_TVector<T>::o_TVector(int n) : nn(n), v(n>0 ? new T[n] : NULL) {
  //for(int i=0; i<n; i++) v[i] = T(0);
}

template <class T>
o_TVector<T>::o_TVector(int n, const T& a) : nn(n), v(n>0 ? new T[n] : NULL)
{
  for(int i=0; i<n; i++) v[i] = a;
}

template <class T>
o_TVector<T>::o_TVector(int n, const T *a) : nn(n), v(n>0 ? new T[n] : NULL)
{
  for(int i=0; i<n; i++) v[i] = *a++;
}

template <class T>
o_TVector<T>::o_TVector(const o_TVector<T> &rhs) : nn(rhs.nn), v(nn>0 ? new T[nn] : NULL)
{
  for(int i=0; i<nn; i++) v[i] = rhs[i];
}

template <class T>
o_TVector<T> & o_TVector<T>::operator=(const o_TVector<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
  if (this != &rhs)
  {
    if (nn != rhs.nn) {
      if (v != NULL) delete [] (v);
      nn=rhs.nn;
      v= nn>0 ? new T[nn] : NULL;
    }
    for (int i=0; i<nn; i++)
      v[i]=rhs[i];
  }
  return *this;
}

template <class T>
inline T & o_TVector<T>::operator[](const int i)	//subscripting
{
#ifdef _CHECKBOUNDS_
  if (i<0 || i>=nn) {
    throw("o_TVector subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
inline const T & o_TVector<T>::operator[](const int i) const	//subscripting
{
#ifdef _CHECKBOUNDS_
  if (i<0 || i>=nn) {
    throw("o_TVector subscript out of bounds");
  }
#endif
  return v[i];
}

template <class T>
inline int o_TVector<T>::size() const
{
  return nn;
}

template <class T>
void o_TVector<T>::resize(int newn)
{
  if (newn != nn) {
    if (v != NULL) delete[] (v);
    nn = newn;
    v = nn > 0 ? new T[nn] : NULL;
  }
}

template <class T>
void o_TVector<T>::assign(int newn, const T& a)
{
  if (newn != nn) {
    if (v != NULL) delete[] (v);
    nn = newn;
    v = nn > 0 ? new T[nn] : NULL;
  }
  for (int i=0;i<nn;i++) v[i] = a;
}

template <class T>
o_TVector<T>::~o_TVector()
{
  if (v != NULL) delete[] (v);
}

#endif
