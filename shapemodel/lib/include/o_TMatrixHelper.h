#ifndef	_o_TMatrixHelper_incl
#define	_o_TMatrixHelper_incl

#include <o_TMatrix.h>
#include <iostream>
#include <fstream>

///////////////////////////////////////////////////////////////////////////////////////////////////
// matrix helper functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// MRmatrix helper functions 
//
// operator+       , matrix matrix addition
// operator-       , matrix matrix substraction
// operator+=      , matrix matrix addition
// operator-=      , matrix matrix substraction
// operator*       , matrix matrix multiply
// operator*       , matrix vector multiply
// operator*       , matrix scalar multiply
// operator*=      , matrix scalar multiply
// compare         , comparison of matrices
// transpose       , transpose matrix
// transpose       , transpose vector
// crosspro_mat    , crossproduct as matrix
// submatrix       , get sub matrix
// getcolumn       , get column
// setcolumn       , set column
// getrow          , get row
// setrow          , set row
// frobenius_norm  , frobenius_norm
// operator<<      , matrix pretty output
// write_ascii     , matrix ascii output
// read_ascii      , matrix ascii input
// identity        , generate an identity matrix
// det             , determinate of matrix
// copycolumn      , copy a colum into another column of matrix
/**
  matrix matrix addition
*/
template <class T>
o_TMatrix<T> operator+(const o_TMatrix<T> &A, const o_TMatrix<T> &B)
{
  int m = A.nrows();
  int n = A.ncols();
  
  if (B.nrows() != m ||  B.ncols() != n ) {
    throw("operator+: matrix dimensions do not match");
    return o_TMatrix<T>();
  }
  else
  {
    o_TMatrix<T> C(m,n);
    
    for (int i=0; i<m; i++)
    {
      for (int j=0; j<n; j++)
        C[i][j] = A[i][j] + B[i][j];
    }
    return C;
  }
}

/**
  matrix matrix substraction
*/
template <class T>
o_TMatrix<T> operator-(const o_TMatrix<T> &A, const o_TMatrix<T> &B)
{
  int m = A.nrows();
  int n = A.ncols();
  
  if (B.nrows() != m ||  B.ncols() != n ) {
    throw("operator+: matrix dimensions do not match");
    return o_TMatrix<T>();
  }
  else
  {
    o_TMatrix<T> C(m,n);
    
    for (int i=0; i<m; i++)
    {
      for (int j=0; j<n; j++)
        C[i][j] = A[i][j] - B[i][j];
    }
    return C;
  }
}

/**
  matrix matrix addition
*/
template <class T>
o_TMatrix<T>&  operator+=(o_TMatrix<T> &A, const o_TMatrix<T> &B)
{
  int m = A.nrows();
  int n = A.ncols();
  
  if (B.nrows() == m ||  B.ncols() == n )
  {
    for (int i=0; i<m; i++)
    {
      for (int j=0; j<n; j++)
        A[i][j] += B[i][j];
    }
  }
  return A;
}

/**
  matrix matrix substraction 
*/
template <class T>
o_TMatrix<T>&  operator-=(o_TMatrix<T> &A, const o_TMatrix<T> &B)
{
  int m = A.nrows();
  int n = A.ncols();
  
  if (B.nrows() == m ||  B.ncols() == n )
  {
    for (int i=0; i<m; i++)
    {
      for (int j=0; j<n; j++)
        A[i][j] -= B[i][j];
    }
  }
  return A;
}

/**
  matrix matrix multiply 
*/
template <class T>
o_TMatrix<T> operator*(const o_TMatrix<T> &A, const o_TMatrix<T> &B)
{
  if (A.ncols() != B.nrows()) {
    throw("operator*: matrix dimensions do not match");
    return o_TMatrix<T>();
  }
  int M = A.nrows();
  int N = A.ncols();
  int K = B.ncols();
  
  o_TMatrix<T> C(M,K);
  
  T sum = 0;
  
  for (int i=0; i<M; i++)
    for (int j=0; j<K; j++)
    {
      sum = 0;
      
      for (int k=0; k<N; k++)
        sum += A[i][k] * B [k][j];
      
      C[i][j] = sum;
    }
  return C;
}

/**
  matrix vector multiply 
*/
template <class T>
o_TVector<T> operator*(const o_TMatrix<T> &A, const o_TVector<T> &B)
{
  if (A.ncols() != B.size()) {
    throw("operator*: matrix and vector dimensions do not match");
    return o_TVector<T>();
  }
  int M = A.nrows();
  int N = A.ncols();
  
  o_TVector<T> C(M);
  
  for (int i=0; i<M; i++) {
    T sum = T(0);
    for (int k=0; k<N; k++) {
      sum += A[i][k] * B [k];
    }
    C[i] = sum;
  }
  return C;
}

/**
  matrix scalar multiply
*/
template <class T>
o_TMatrix<T> operator*(const o_TMatrix<T> &A, const T &B)
{
  int m = A.nrows();
  int n = A.ncols();
  
  o_TMatrix<T> C(m,n);
  
  for (int i=0; i<m; i++)
  {
    for (int j=0; j<n; j++)
      C[i][j] = A[i][j] * B;
  }
  return C;
}

/**
  matrix scalar multiply
*/
template <class T>
o_TMatrix<T>&  operator*=(o_TMatrix<T> &A, const T &B)
{
  int m = A.nrows();
  int n = A.ncols();
  
  for (int i=0; i<m; i++)
  {
    for (int j=0; j<n; j++)
      A[i][j] *= B;
  }
  return A;
}

/**
  comparison of matrices
*/
template <class T>
bool compare(const o_TMatrix<T> &A, const o_TMatrix<T> &B, T error_threshhold=T(1e-8))
{
  if (A.ncols() != B.ncols() || A.nrows() != B.nrows())
    return false;
  
  int M = A.nrows();
  int N = A.ncols();
  
  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
      if ( A[i][j] > B[i][j] && A[i][j] - B[i][j] > error_threshhold ) return false;
      if ( B[i][j] > A[i][j] && B[i][j] - A[i][j] > error_threshhold ) return false;
    }
  }
  return true;
}

/**
  transpose matrix 
*/
template <class T>
o_TMatrix<T> transpose(const o_TMatrix<T> &A)
{
  int m = A.nrows();
  int n = A.ncols();
  
  o_TMatrix<T> C(n,m);
  
  for (int i=0; i<m; i++)
  {
    for (int j=0; j<n; j++)
      C[j][i] = A[i][j];
  }
  return C;
}

/**
  transpose vector (returning a matrix)
*/
template <class T>
o_TMatrix<T> transpose(const o_TVector<T> &A)
{
  int m = A.size();
  
  o_TMatrix<T> C(1,m);
  
  for (int i=0; i<m; i++)
  {
    C[0][i] = A[i];
  }
  return C;
}

/**
  crossproduct as matrix
*/
template <class T>
o_TMatrix<T> crosspro_mat(const o_TVector<T> &A)
{
  o_TMatrix<T> out(3,3);
  
  if (A.size() != 3 ) {
    throw("crosspro_mat: vector dimension not 3");
    return out;
  }
  
  out[0][0] =  T(0); out[0][1] = -A[2]; out[0][2] =  A[1];
  out[1][0] =  A[2]; out[1][1] =  T(0); out[1][2] = -A[0];
  out[2][0] = -A[1]; out[2][1] =  A[0]; out[2][2] =  T(0);
  
  return out;
}


/**
  submatrix
*/
template <class T>
o_TMatrix<T> submatrix(unsigned start_row, unsigned start_col, unsigned num_rows,
                       unsigned num_cols, const o_TMatrix<T> &A)
{
  
  o_TMatrix<T> C(num_rows,num_cols);
  
  for (unsigned i=0; i<num_rows; i++)
  {
    for (unsigned j=0; j<num_cols; j++){
      C[i][j] = A[start_row+i][start_col+j];
    }
  }
  return C;
}

/**
  get column
*/
template <class T>
o_TVector<T> getcolumn(unsigned col, const o_TMatrix<T> &A)
{
  int m = A.nrows();
  o_TVector<T> C(m);
  
  for (int i=0; i<m; i++)
  {
    C[i] = A[i][col];
  }
  return C;
}

/**
  set column
*/
template <class T>
void setcolumn(unsigned col, const o_TVector<T> A, o_TMatrix<T> &B )
{
  if (B.nrows() != A.size()) {
    throw("setcolumn: matrix and vector dimensions do not match");
    return;
  }
  int m = B.nrows();
  
  for (int i=0; i<m; i++)
  {
    B[i][col] = A[i];
  }
}

/**
  set column
*/
template <class T>
void setcolumn(unsigned col_in, const o_TMatrix<T> A, unsigned col_out, o_TMatrix<T> &B )
{
  if (B.nrows() != A.nrows()) {
    throw("setcolumn: matrix and matrix dimensions do not match");
    return;
  }
  int m = B.nrows();
  
  for (int i=0; i<m; i++)
  {
    B[i][col_out] = A[i][col_in];
  }
}


/**
  get row
*/
template <class T>
o_TVector<T> getrow(unsigned row, const o_TMatrix<T> &A)
{
  int n = A.ncols();
  o_TVector<T> C(n);
  
  for (int i=0; i<n; i++)
  {
    C[i] = A[row][i];
  }
  return C;
}

/**
  set column
*/
template <class T>
void setrow(unsigned row, const o_TVector<T> A, o_TMatrix<T> &B )
{
  if (B.ncols() != A.size()) {
    throw("setcolumn: matrix and vector dimensions do not match");
    return;
  }
  
  int m = B.ncols();
  
  for (int i=0; i<m; i++)
  {
    B[row][i] = A[i];
  }
}

/**
  frobenius norm of a matrix
*/

template< class T >
T frobenius_norm(const o_TMatrix<T> &in) {
  T result= T(0);
  for (int i0 = 0; i0 < in.nrows(); i0++) {
    for (int i1 = 0; i1 < in.ncols(); i1++) {
      result += (in[i0][i1] * in[i0][i1]);
    }
  }
  return sqrt(result);
}

/**
  matrix pretty output
*/
template <class T>
std::ostream& operator<<(std::ostream &s, const o_TMatrix<T> &A)
{
  int M=A.nrows();
  int N=A.ncols();
  s.setf(std::ios::fixed);
  for (int i=0; i<M; i++)
  {
    s << "[ ";
    for (int j=0; j<N; j++)
    {
      s << A[i][j] << " ";
    }
    
    s << "]\n";
  }
  s.unsetf(std::ios::fixed);
  return s;
}

/**
  matrix output
*/
template <class T>
std::ostream& write_ascii(std::ostream &s, const o_TMatrix<T> &A)
{
  int M=A.nrows();
  int N=A.ncols();
  streamsize old_prec = s.precision(16);
  s << M << " " << N << "\n";
  
  for (int i=0; i<M; i++)
  {
    for (int j=0; j<N; j++)
    {
      s << A[i][j] << " ";
    }
    s << "\n";
  }
  s.precision(old_prec);
  return s;
}

/**
  matrix input
*/
template <class T>
std::istream& read_ascii(std::istream &s, o_TMatrix<T> &A)
{
  
  int M, N;
  
  s >> M >> N;
  
  o_TMatrix<T> B(M,N);
  
  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++)
    {
      s >>  B[i][j];
    }
  }
  A = B;
  return s;
}


/**
  matrix binary output
*/
template <class T> bool write_binary(std::ostream &s, const o_TMatrix< T > &m) {
  unsigned row, col;
  row = m.nrows() ; col = m.ncols() ;
  if(!s.write((char*)&"nr_matrix",sizeof(char)*9)) return 0 ;
  if(!s.write((char*)&row,sizeof(unsigned))) return 0 ;
  if(!s.write((char*)&col,sizeof(unsigned))) return 0 ;
  for(unsigned r = 0; r < row; r++){
    for(unsigned c = 0; c < col; c++){
      T tmp = m[r][c];
      if(!s.write((char*)&tmp,sizeof(T))) return 0 ;
    }
  }
  return 1;
}

/**
  matrix binary input
*/
template <class T> 
bool read_binary(std::istream &s, o_TMatrix< T > &m) {
  
  unsigned row, col;
  char *type ;
  type = new char[9] ;
  if(!s.read(type,sizeof(char)*9)) return false;
  int suc = strncmp(type,"nr_matrix",9) ;
  if(suc) return 0 ;
  if(!s.read((char*)&row,sizeof(unsigned))) return false;
  if(!s.read((char*)&col,sizeof(unsigned))) return false;
  m.resize(row,col) ;
  for(unsigned r = 0; r < row; r++){
    for(unsigned c = 0; c < col; c++){
      if(!s.read((char*)&m[r][c],sizeof(T))) return false;
    }
  }
  delete []type ;
  return true;
}

/**
  identity matrix
*/
template <class T>
void identity(o_TMatrix<T> &A, int s)
{

  if (A.nrows() != s ||  A.ncols() != s ) {
    A.resize(s,s);
  }
  
  for (int i=0; i<s; i++){
    for (int j=0; j<s; j++){
      if(i==j) A[i][j] = T(1); else A[i][j] = T(0);
    }
  }
}

template< class T > T det(const o_TMatrix<T> &in)
{
  T result = T(0);
  
  if (in.size1() == 0) return result;
  if (in.size1() != in.size2())
  {
    throw(" determinant not defined for non square matrix");
    return T(0);
  }
  
  if (in.size1() == 1)
    return in [0] [0];
  else if (in.size1() == 2)
    return in [0] [0] * in [1] [1] - in [0] [1] * in [1] [0];
  else if (in.size1() == 3)
  {
    for (unsigned i0 = 0; i0 < 3; i0++)
    {
      result += in [0] [i0] * in [1] [(i0 + 1) % 3] * in [2] [(i0 + 2) % 3];
      result -= in [2] [i0] * in [1] [(i0 + 1) % 3] * in [0] [(i0 + 2) % 3];
    }
    
    return result;
  }
  else
  {
    o_TMatrix<T> decreased_m(in.size1() - 1, in.size1() - 1);
    
    for (unsigned i0 = 0; i0 < in.size1(); i0++)
    {
      for (unsigned i1 = 0; i1 < decreased_m.size1(); i1++)
        for (unsigned i2 = 0; i2 < decreased_m.size1(); i2++)
          decreased_m [i1] [i2] = in [i1 + 1] [(i2 >= i0 ? i2 + 1 : i2)];
      
       // speed ...
      if (in [0] [i0] != T(0))
        result += in [0] [i0] * (i0 % 2 ? -1 : 1) * det(decreased_m);
    }
    return result;
  }
}

template< class T > void copycolumn(unsigned source_col, const o_TMatrix<T> &source, unsigned dest_col, o_TMatrix<T> &dest)
{  
  for (unsigned i0 = 0; i0 < source.size1() && i0 < dest.size1(); i0++){
    dest [i0] [dest_col] = source [i0] [source_col];
  }
}

#endif
