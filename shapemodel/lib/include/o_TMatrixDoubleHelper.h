
//**************************************************************************************
// functions overview
//**************************************************************************************
//
// MRmatrix <double> helper functions 
//
// invert_gaussj         , matrix inversion with gaussj
// pseudo_inverse        , calculates the pseudo inverse

#ifndef	_o_TMatrixDoubleHelper_incl
#define	_o_TMatrixDoubleHelper_incl

#include <o_TVector.h>
#include <o_TVectorHelper.h>
#include <o_TMatrix.h>
#include <o_TMatrixHelper.h>
#include <nr3.h>
#include <gaussj.h>
#include <svd.h>

/**
  matrix invert gaussj
*/
template< class T >
NRmatrix< double > invert_gaussj(const NRmatrix< T > &A)
{
  int m = A.nrows();
  int n = A.ncols();
  
  if(m != n){
    throw("invert_gaussj: matrix is not squared");
  }
  
  NRmatrix < double > C (A);
  
  nr::gaussj(C);
  
  return C;
}


/** 
    Calculate the pseudo inverse of a matrix.
    Give the rank of the matrix,  if you know it.
    Otherwise set rank to 0 and give a tolerance threshhold for the singular values that will not be zeroed.
    If you have no idea of the treshhold set it to a negative value, 
    then the default tolerance threshhold is MAX(SIZE(pseudo)) * NORM(pseudo) * 2.2204e-16.
    The return value gives the number of zeroed sigular values.
*/

template< class T >
unsigned pseudo_inverse(NRmatrix < T > & pseudo, unsigned rank, T singular_value_threshhold)
{
  //cout.setf(ios::scientific, ios::floatfield); 
  nr::SVD mysvd( pseudo ); 
  NRmatrix< double > Wmat_plus(mysvd.w.size(), mysvd.w.size()); 
  
  unsigned ret = 0;
  
  if(rank != 0){
    for(int b = 0; b < mysvd.w.size(); b ++){ 
      for(int c = 0; c < mysvd.w.size(); c ++){ 
        if(b == c){
          if(c < (int)rank){ 
            Wmat_plus[c][c] = 1.0 / mysvd.w[c];
          }else{
            Wmat_plus[c][c] = 0.0;
            ret++;
          }
        }else{ 
          Wmat_plus[b][c] = 0.0; 
        }	 
      }  
    }
  }else{
    for(int b = 0; b < mysvd.w.size(); b ++){ 
      for(int c = 0; c < mysvd.w.size(); c ++){ 
        if(b == c){
          if(mysvd.w[c] >  mysvd.tsh){ 
            Wmat_plus[c][c] = 1.0 / mysvd.w[c];
          }else{
            Wmat_plus[c][c] = 0.0;
            ret++;
          }
        }else{ 
          Wmat_plus[b][c] = 0.0; 
        }	 
      }  
    }
  }

  pseudo = mysvd.v * Wmat_plus * transpose(mysvd.u);
  
  return ret;
}

#endif
