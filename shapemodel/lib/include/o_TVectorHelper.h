#ifndef	_o_TVectorHelper_incl
#define	_o_TVectorHelper_incl

#include <vector>
#include <o_TVector.h>

//
// o_TVector helper functions
// 
// operator+       , vector vector addition
// operator-       , vector vector substraction
// operator+=      , vector vector addition
// operator-=      , vector vector substraction
// operator*       , vector vector multiply
// operator*       , vector scalar multiply
// operator*=      , vector scalar multiply
// compare         , comparison of vectors
// homogenize      , homogenize vector
// norm            , vector norm
// norm2           , vector norm^2
// cross product   , cross product
// operator<<      , vector pretty output
// write_ascii     , vector ascii output
// read_ascii      , vector ascii input
// cast2std        , cast to std:::vector


///////////////////////////////////////////////////////////////////////////////////////////////////
// vector helper functions
///////////////////////////////////////////////////////////////////////////////////////////////////

/**
  vector vector addition
*/
template <class T>
o_TVector<T> operator+(const o_TVector<T> &A, const o_TVector<T> &B)
{
  int n = A.size();
  
  if (B.size() != n ){
    throw("operator+: vector dimensions do not match");
    return o_TVector<T>();
  } else {
    o_TVector<T> C(n);
    
    for (int i=0; i<n; i++)
    {
      C[i] = A[i] + B[i];
    }
    return C;
  }
}

/**
  vector vector substraction
*/
template <class T>
o_TVector<T> operator-(const o_TVector<T> &A, const o_TVector<T> &B)
{
  int n = A.size();
  
  if (B.size() != n ){
    throw("operator-: vector dimensions do not match");
    return o_TVector<T>();
  } else {
    o_TVector<T> C(n);
    
    for (int i=0; i<n; i++)
    {
      C[i] = A[i] - B[i];
    }
    return C;
  }
}

/**
  vector vector addition
*/
template <class T>
o_TVector<T>&  operator+=(o_TVector<T> &A, const o_TVector<T> &B)
{
  int n = A.size();
  
  if (B.size() == n)
  {
    for (int i=0; i<n; i++)
    {
      A[i] += B[i];
    }
  } else {
    throw("operator+=: vector dimensions do not match");
  }
  return A;
}

/**
  vector vector substraction
*/
template <class T>
o_TVector<T>&  operator-=(o_TVector<T> &A, const o_TVector<T> &B)
{
  int n = A.size();
  
  if (B.size() == n)
  {
    for (int i=0; i<n; i++)
    {
      A[i] -= B[i];
    }
  } else {
    throw("operator-=: vector dimensions do not match");
  }
  return A;
}

/**
  vector vector multiply, A is assumed transposed, therefore the result is a scalar
*/
template <class T>
T operator*(const o_TVector<T> &A, const o_TVector<T> &B)
{
  if (A.size() != B.size()) {
    throw("operator*: vector dimensions do not match");
    return T(0);
  }
  T C = T(0);
  
  for (int i=0; i<A.size(); i++)
  {
    C += A[i] * B[i];
  }
  return C;
}

/**
   vector scalar multiply
*/
template <class T>
o_TVector<T> operator*(const o_TVector<T> &A, const T &B)
{
  int n = A.size();
  
  o_TVector<T> C(n);
  
  for (int i=0; i<n; i++)
  {
    C[i] = A[i] * B;
  }
  return C;
}

/**
  vector scalar multiply
*/
template <class T>
o_TVector<T>&  operator*=(o_TVector<T> &A, const T &B)
{
  int n = A.size();
  
  for (int i=0; i<n; i++)
  {
    A[i] *= B;
  }
  return A;
}

/**
  comparison of vectors
*/
template <class T>
bool compare(const o_TVector<T> &A, const o_TVector<T> &B,  T error_threshhold=T(1e-8))
{
  if ( A.size() != B.size() )
    return false;
  
  int N=A.size();
  for (int i=0; i<N; i++) {
    if ( A[i] > B[i] && A[i] - B[i] > error_threshhold ) return false;
    if ( B[i] > A[i] && B[i] - A[i] > error_threshhold ) return false;
  }
  return true;
}

/**
  homogenize vector 
*/
template <class T>
o_TVector<T>&  homogenize(o_TVector<T> &A)
{
  int n = A.size();
  
  for (int i=0; i<n-1; i++)
  {
    A[i] /= A[n-1];
  }
  A[n-1] = T(1);
  return A;
}

/**
  vector norm
*/
template <class T>
T norm(const o_TVector<T> &A)
{
  return sqrt(norm2(A));
}

/**
  vector norm^2
*/
template <class T>
T norm2(const o_TVector<T> &A)
{
  T out = T(0);
  for (int i=0; i<A.size(); i++)
  {
    out += A[i] * A[i];
  }
  return out;
}

/**
  cross product
*/
template< class T >
o_TVector< T > crossproduct(const o_TVector<T> &left, const o_TVector<T> &right)
{
  o_TVector<T> result(3);
  if (left.size() != 3 || right.size() != 3)
  {
    throw ("crossproduct: dimension of vectors must be 3");
    return result;
  }
  result [0] = left[1]*right[2] - left[2]*right[1];
  result [1] = left[2]*right[0] - left[0]*right[2];
  result [2] = left[0]*right[1] - left[1]*right[0];
  return result;
}


/**
  vector pretty output
*/
template <class T>
std::ostream& operator<<(std::ostream &s, const o_TVector<T> &A)
{
  s.setf(std::ios::fixed);
  s << "[ ";
  for (int j=0; j<A.size(); j++)
  {
    s << A[j] << " ";
  }
  s << "]" << endl;
  s.unsetf(std::ios::fixed);
  return s;
}

/**
  vector output
*/
template <class T>
std::ostream& write_ascii(std::ostream &s, const o_TVector<T> &A)
{
  int N=A.size();
  streamsize old_prec = s.precision(16);
  s << N << "\n";
  for (int j=0; j<N; j++)
  {
    s << A[j] << "\n";
  }
  s << "\n";
  s.precision(old_prec);
  return s;
}

/**
  vector input
*/
template <class T>
std::istream& read_ascii(std::istream &s, o_TVector<T> &A)
{
  int N;
  s >> N;
  
  o_TVector<T> B(N);
  for (int i=0; i<N; i++)
    s >> B[i];
  A = B;
  return s;
}

/**
  vector binary output
*/
template <class T> bool write_binary(std::ostream &s, const o_TVector< T > &m) {
  unsigned col;
  col = m.size();
  if(!s.write((char*)&"nr_matrix",sizeof(char)*9)) return 0 ;
  if(!s.write((char*)&col,sizeof(unsigned))) return 0 ;
  for(unsigned c = 0; c < col; c++){
    T tmp = m[c];
    if(!s.write((char*)&tmp,sizeof(T))) return 0 ;
  }
  return 1;
}

/**
  matrix binary input
*/
template <class T> 
bool read_binary(std::istream &s, o_TVector< T > &m) {
  
  unsigned col;
  char *type ;
  type = new char[9] ;
  if(!s.read(type,sizeof(char)*9)) return false;
  int suc = strncmp(type,"nr_matrix",9) ;
  if(suc) return 0 ;
  if(!s.read((char*)&col,sizeof(unsigned))) return false;
  m.resize(col) ;
  for(unsigned c = 0; c < col; c++){
    if(!s.read((char*)&m[c],sizeof(T))) return false;
  }
  delete []type ;
  return true;
}
/**
  cast to std::vector
*/
template <class T>
std::vector<T> cast2std(const o_TVector<T> &B)
{
  std::vector<T> A(B.size());
  
  for (int i=0; i<B.size(); i++)
  {
      A[i] = B[i];
  }
  return A;
}



#endif
