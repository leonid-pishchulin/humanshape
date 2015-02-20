/**
    This file is part of the implementation of the 3D human shape model as described in the paper:

    Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobald and Bernt Schiele
    Building Statistical Shape Spaces for 3D Human Modeling
    ArXiv, February 2015

    Please cite the paper if you are using this code in your work.
    
    Contributors: Arjun Jain, Juergen Gall, Leonid Pishchulin.
    This code also uses open source functionality for matrix operations by Thomas Brox.

    The code may be used free of charge for non-commercial and
    educational purposes, the only requirement is that this text is
    preserved within the derivative work. For any other purpose you
    must contact the authors for permission. This code may not be
    redistributed without permission from the authors.
*/

#include <math.h>
#include <stdlib.h>
#include "NMath.h"
#include <limits.h>

namespace NMath {

  const float Pi = 3.1415926536;

  // faculty
  int faculty(int n) {
    int aResult = 1;
    for (int i = 2; i <= n; i++)
      aResult *= i;
    return aResult;
  }

  // binCoeff
  int binCoeff(const int n, const int k) {
    if (k > (n >> 1)) return binCoeff(n,n-k);
    int aResult = 1;
    for (int i = n; i > (n-k); i--)
      aResult *= i;
    for (int j = 2; j <= k; j++)
      aResult /= j;
    return aResult;
  }

  // tangent
  float tangent(const float x1, const float y1, const float x2, const float y2) {
    float alpha;
    float xDiff = x2-x1;
    float yDiff = y2-y1;
    if (yDiff > 0) {
      if (xDiff == 0) alpha = 0.5*Pi;
      else if (xDiff > 0) alpha = atan(yDiff/xDiff);
      else alpha = Pi+atan(yDiff/xDiff);
    }
    else {
      if (xDiff == 0) alpha = -0.5*Pi;
      else if (xDiff > 0) alpha = atan(yDiff/xDiff);
      else alpha = -Pi+atan(yDiff/xDiff);
    }
    return alpha;
  }

  // absAngleDifference
  float absAngleDifference(const float aFirstAngle, const float aSecondAngle) {
    float aAlphaDiff = abs(aFirstAngle - aSecondAngle);
    if (aAlphaDiff > Pi) aAlphaDiff = 2*Pi-aAlphaDiff;
    return aAlphaDiff;
  }

  // angleDifference
  float angleDifference(const float aFirstAngle, const float aSecondAngle) {
    float aAlphaDiff = aFirstAngle - aSecondAngle;
    if (aAlphaDiff > Pi) aAlphaDiff = -2*Pi+aAlphaDiff;
    else if (aAlphaDiff < -Pi) aAlphaDiff = 2*Pi+aAlphaDiff;
    return aAlphaDiff;
  }

  // angleSum
  float angleSum(const float aFirstAngle, const float aSecondAngle) {
    float aSum = aFirstAngle + aSecondAngle;
    if (aSum > Pi) aSum = -2*Pi+aSum;
    else if (aSum < -Pi) aSum = 2*Pi+aSum;
    return aSum;
  }

  // round
  int round(const float aValue) {
    float temp1 = floor(aValue);
    float temp2 = ceil(aValue);
    if (aValue-temp1 < 0.5) return (int)temp1;
    else return (int)temp2;
  }

  // orderedRandom
  void orderedRandom(CVector<double>& aResult) {
    double tmp = 1.0;
    for(int i=aResult.size(); i>0; i--) {
      tmp *= pow( (double)random(), 1.0/(double)i );
      aResult[i-1] = tmp;
    }
  }



  // PATransformation
  // Cyclic Jacobi method for determining the eigenvalues and eigenvectors
  // of a symmetric matrix.
  // Ref.:  H.R. Schwarz: Numerische Mathematik. Teubner, Stuttgart, 1988.
  //        pp. 243-246.
  void PATransformation(const CMatrix<float>& aMatrix, CVector<float>& aEigenvalues, CMatrix<float>& aEigenvectors) {
    static const float eps = 0.0001;
    static const float delta = 0.000001;
    static const float eps2 = eps*eps;
    float sum,theta,t,c,r,s,g,h;
    // Initialization
    CMatrix<float> aCopy(aMatrix);
    int n = aEigenvalues.size();
    aEigenvectors = 0;
    for (int i = 0; i < n; i++)
      aEigenvectors(i,i) = 1;
    // Loop
    do {
      // check whether accuracy is reached
      sum = 0.0;
      for (int i = 1; i < n; i++)
        for (int j = 0; j <= i-1; j++)
          sum += aCopy(i,j)*aCopy(i,j);
      if (sum+sum > eps2) {
        for (int p = 0; p < n-1; p++)
          for (int q = p+1; q < n; q++)
            if (fabs(aCopy(q,p)) >= eps2) {
              theta = (aCopy(q,q) - aCopy(p,p)) / (2.0 * aCopy(q,p));
              t = 1.0;
              if (fabs(theta) > delta) t = 1.0 / (theta + theta/fabs(theta) * sqrt (theta*theta + 1.0));
              c = 1.0 / sqrt (1.0 + t*t);
              s = c*t;
              r = s / (1.0 + c);
              aCopy(p,p) -= t * aCopy(q,p);
              aCopy(q,q) += t * aCopy(q,p);
              aCopy(q,p) = 0;
              for (int j = 0; j <= p-1; j++) {
                g = aCopy(q,j) + r * aCopy(p,j);
                h = aCopy(p,j) - r * aCopy(q,j);
                aCopy(p,j) -= s*g;
                aCopy(q,j) += s*h;
              }
              for (int i = p+1; i <= q-1; i++) {
                g = aCopy(q,i) + r * aCopy(i,p);
                h = aCopy(i,p) - r * aCopy(q,i);
                aCopy(i,p) -= s * g;
                aCopy(q,i) += s * h;
              }
              for (int i = q+1; i < n; i++) {
                g = aCopy(i,q) + r * aCopy(i,p);
                h = aCopy(i,p) - r * aCopy(i,q);
                aCopy(i,p) -= s * g;
                aCopy(i,q) += s * h;
              }
              for (int i = 0; i < n; i++) {
                g = aEigenvectors(i,q) + r * aEigenvectors(i,p);
                h = aEigenvectors(i,p) - r * aEigenvectors(i,q);
                aEigenvectors(i,p) -= s * g;
                aEigenvectors(i,q) += s * h;
              }
            }
      }
    }
    // Return eigenvalues
    while (sum+sum > eps2);
    for (int i = 0; i < n; i++)
      aEigenvalues(i) = aCopy(i,i);
    // Order eigenvalues and eigenvectors
    for (int i = 0; i < n-1; i++) {
      int k = i;
      for (int j = i+1; j < n; j++)
        if (fabs(aEigenvalues(j)) > fabs(aEigenvalues(k))) k = j;
      if (k != i) {
        // Switch eigenvalue i and k
        float help = aEigenvalues(k);
        aEigenvalues(k) = aEigenvalues(i);
        aEigenvalues(i) = help;
        // Switch eigenvector i and k
        for (int j = 0; j < n; j++) {
          help = aEigenvectors(j,k);
          aEigenvectors(j,k) = aEigenvectors(j,i);
          aEigenvectors(j,i) = help;
        }
      }
    }
  }

  // PABackTransformation
  void PABacktransformation(const CMatrix<float>& aEigenvectors, const CVector<float>& aEigenvalues, CMatrix<float>& aMatrix) {
    for (int i = 0; i < aEigenvalues.size(); i++)
      for (int j = 0; j <= i; j++) {
         float sum = aEigenvalues(0) * aEigenvectors(i,0) * aEigenvectors(j,0);
         for (int k = 1; k < aEigenvalues.size(); k++)
           sum += aEigenvalues(k) * aEigenvectors(i,k) * aEigenvectors(j,k);
         aMatrix(i,j) = sum;
      }
    for (int i = 0; i < aEigenvalues.size(); i++)
      for (int j = i+1; j < aEigenvalues.size(); j++)
        aMatrix(i,j) = aMatrix(j,i);
  }

  // svd (übernommen von Bodo Rosenhahn)
  void svd(CMatrix<float>& U, CMatrix<float>& S, CMatrix<float>& V, bool aOrdering, int aIterations) {
    static float at,bt,ct;
    static float maxarg1,maxarg2;
    #define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ?  (ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))
    #define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?	(maxarg1) : (maxarg2))
    #define MIN(a,b) ((a) >(b) ? (b) : (a))
    #define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
    int flag,i,its,j,jj,k,l,nm;
	  float c,f,h,s,x,y,z;
	  float anorm=0.0,g=0.0,scale=0.0;
    int aXSize = U.xSize();
    int aYSize = U.ySize();
	  CVector<float> aBuffer(aXSize);
    for (i = 0; i < aXSize; i++) {
      l=i+1;
      aBuffer(i)=scale*g;
      g=s=scale=0.0;
      if (i < aYSize) {
        for (k = i; k < aYSize; k++)
          scale += fabs(U(i,k));
        if (scale) {
          for (k = i; k < aYSize; k++) {
            U(i,k) /= scale;
	   				s += U(i,k)*U(i,k);
          }
		  		f = U(i,i);
          g = -SIGN(sqrt(s),f);
		  		h = f*g-s;
          U(i,i) = f-g;
          for (j = l; j < aXSize; j++) {
  	   			for (s = 0.0, k = i; k < aYSize; k++)
              s += U(i,k)*U(j,k);
            f = s/h;
            for (k = i; k < aYSize; k++)
				     	U(j,k) += f*U(i,k);
          }
			    for ( k = i; k < aYSize; k++)
			     	U(i,k) *= scale;
        }
      }
   	  S(i,i) = scale*g;
      g=s=scale=0.0;
      if (i < aYSize && i != aXSize-1) {
     	  for (k = l; k < aXSize; k++)
          scale += fabs(U(k,i));
        if (scale != 0)	{
          for (k = l; k < aXSize; k++) {
            U(k,i) /= scale;
            s += U(k,i)*U(k,i);
          }
	   		  f = U(l,i);
          g = -SIGN(sqrt(s),f);
	   		  h = f*g-s;
          U(l,i) = f-g;
          for (k = l; k < aXSize; k++)
            aBuffer(k) = U(k,i)/h;
          for (j = l; j < aYSize; j++) {
            for (s = 0.0, k = l; k < aXSize; k++)
              s += U(k,j)*U(k,i);
	   			  for (k = l; k < aXSize; k++)
              U(k,j) += s*aBuffer(k);
          }
	   		  for (k = l; k < aXSize; k++)
            U(k,i) *= scale;
        }
      }
	    anorm = MAX(anorm,(fabs(S(i,i))+fabs(aBuffer(i))));
    }
   	for (i = aXSize-1; i >= 0; i--)	{
   		if (i < aXSize-1)	{
        if (g != 0)	{
          for (j = l; j < aXSize; j++)
            V(i,j) = U(j,i)/(U(l,i)*g);
          for (j = l; j < aXSize; j++) {
            for (s = 0.0, k = l; k < aXSize; k++)
              s += U(k,i)*V(j,k);
            for (k = l; k < aXSize; k++)
              V(j,k) += s*V(i,k);
          }
        }
        for (j = l; j < aXSize; j++)
          V(j,i) = V(i,j) = 0.0;
      }
  		V(i,i) = 1.0;
      g = aBuffer(i);
  		l = i;
    }
	  for (i = MIN(aYSize-1,aXSize-1); i >= 0; i--)	{
	  	l = i+1;
	  	g = S(i,i);
      for (j = l; j < aXSize; j++)
	  		U(j,i) = 0.0;
      if (g != 0) {
        g = 1.0/g;
        for (j = l; j < aXSize; j++) {
          for (s = 0.0, k = l; k < aYSize; k++)
            s += U(i,k)*U(j,k);
		   		f = (s/U(i,i))*g;
          for (k = i; k < aYSize; k++)
            U(j,k) += f*U(i,k);
        }
   			for (j = i; j < aYSize; j++)
          U(i,j) *= g;
      }
		  else {
		   	for (j = i; j < aYSize; j++)
			  	U(i,j) = 0.0;
      }
   		++U(i,i);
    }
   	for (k = aXSize-1; k >= 0; k--)	{
      for (its = 1; its <= aIterations; its++)	{
	   		flag = 1;
        for (l = k; l >= 0; l--) {
          nm = l - 1;
          if (fabs(aBuffer(l))+anorm == anorm)	{
            flag = 0; break;
          }
				  if (fabs(S(nm,nm))+anorm == anorm)	break;
        }
   			if (flag)	{
	  		  c = 0.0;
          s = 1.0;
          for (i = l; i <= k; i++) {
            f = s*aBuffer(i);
            aBuffer(i) = c*aBuffer(i);
            if (fabs(f)+anorm == anorm)	break;
            g = S(i,i);
				   	h = PYTHAG(f,g);
            S(i,i) = h;
            h = 1.0/h;
            c = g*h;
            s=-f*h;
            for (j = 0; j < aYSize; j++) {
              y = U(nm,j);
		   				z = U(i,j);
              U(nm,j) = y*c + z*s;
              U(i,j) = z*c - y*s;
            }
          }
        }
	  		z = S(k,k);
	   		if (l == k)	{
          if (z < 0.0) {
            S(k,k) = -z;
            for (j = 0; j < aXSize; j++)
              V(k,j) = -V(k,j);
          }
	   			break;
        }
		   	if (its == aIterations) std::cerr << "svd: No convergence in " << aIterations << " iterations" << std::endl;
			  x = S(l,l);
  			nm = k-1;
        y = S(nm,nm);
        g = aBuffer(nm);
        h = aBuffer(k);
		  	f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
        g = PYTHAG(f,1.0);
        f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
        c = s = 1.0;
   			for (j = l; j <= nm; j++)	{
	  			i = j+1;
          g = aBuffer(i);
          y = S(i,i);
          h = s*g;
          g = c*g;
          z = PYTHAG(f,h);
          aBuffer(j) = z;
          float invZ = 1.0/z;
          c = f*invZ;
          s = h*invZ;
          f = x*c+g*s;
          g = g*c-x*s;
          h = y*s;
          y *= c;
          for (jj = 0; jj < aXSize; jj++)	{
            x = V(j,jj);
		   			z = V(i,jj);
            V(j,jj) = x*c + z*s;
			   		V(i,jj) = z*c - x*s;
          }
  				z = PYTHAG(f,h);
          S(j,j) = z;
          if (z != 0)	{
            z = 1.0/z;
            c = f*z;
            s = h*z;
          }
  				f = (c*g)+(s*y);
          x = (c*y)-(s*g);
          for (jj = 0; jj < aYSize; jj++)	{
            y = U(j,jj);
		  			z = U(i,jj);
            U(j,jj) = y*c + z*s;
            U(i,jj) = z*c - y*s;
          }
        }
   			aBuffer(l) = 0.0;
        aBuffer(k) = f;
	   		S(k,k) = x;
      }
    }
    // Order singular values
    if (aOrdering) {
      for (int i = 0; i < aXSize-1; i++) {
        int k = i;
        for (int j = i+1; j < aXSize; j++)
          if (fabs(S(j,j)) > fabs(S(k,k))) k = j;
        if (k != i) {
          // Switch singular value i and k
          float help = S(k,k);
          S(k,k) = S(i,i);
          S(i,i) = help;
          // Switch columns i and k in U and V
          for (int j = 0; j < aYSize; j++) {
            help = U(k,j);
            U(k,j) = U(i,j);
            U(i,j) = help;
            help = V(k,j);
            V(k,j) = V(i,j);
            V(i,j) = help;
          }
        }
      }
    }
  }

  #undef PYTHAG
  #undef MAX
  #undef MIN
  #undef SIGN

  // svdBack
  void svdBack(CMatrix<float>& U, const CMatrix<float>& S, const CMatrix<float>& V) {
    for (int y = 0; y < U.ySize(); y++)
      for (int x = 0; x < U.xSize(); x++)
        U(x,y) = S(x,x)*U(x,y);
    U *= trans(V);
  }

  // Householder-Verfahren (nach Stoer), uebernommen von Bodo Rosenhahn
  // Bei dem Verfahren wird die Matrix A (hier:*this und die rechte Seite (b)
  // mit unitaeren Matrizen P multipliziert, so dass A in eine
  // obere Dreiecksmatrix umgewandelt wird.
  // Dabei ist P = I + beta * u * uH
  // Die Vektoren u werden bei jeder Transformation in den nicht
  // benoetigten unteren Spalten von A gesichert.

  void householder(CMatrix<float>& A, CVector<float>& b) {
	  int i,j,k;
	  float sigma,s,beta,sum;
	  CVector<float> d(A.xSize());
    for (j = 0; j < A.xSize(); ++j) {
      sigma = 0;
	    for (i = j; i < A.ySize(); ++i)
    	  sigma += A(j,i)*A(j,i);
      if (sigma == 0) {
	      std::cerr << "NMath::householder(): matrix is singular!" << std::endl;
 	      break;
      }
	    // Choose sign to avoid elimination
	    s = d(j) = A(j,j)<0 ? sqrt(sigma) : -sqrt(sigma);
	    beta = 1.0/(s*A(j,j)-sigma);
      A(j,j) -= s;
      // Transform submatrix of A with P
	    for (k = j+1; k < A.xSize(); ++k)	{
	      sum = 0.0;
	      for (i = j; i < A.ySize(); ++i)
		      sum += (A(j,i)*A(k,i));
	      sum *= beta;
	      for (i = j; i < A.ySize(); ++i)
          A(k,i) += A(j,i)*sum;
      }
      // Transform right hand side of linear system with P
	    sum = 0.0;
	    for (i = j; i < A.ySize(); ++i)
	      sum += A(j,i)*b(i);
	    sum *= beta;
	    for (i = j; i < A.ySize(); ++i)
	      b(i) += A(j,i)*sum;
    }
    for (i = 0; i < A.xSize(); ++i)
	    A(i,i) = d(i);
  }

  // leastSquares
  CVector<float> leastSquares(CMatrix<float>& A, CVector<float>& b) {
    CVector<float> aResult(A.xSize());
    householder(A,b);
    for (int i = A.xSize()-1; i >= 0; i--) {
      float s = 0;
	    for (int k = i+1; k < A.xSize(); k++)
	      s += A(k,i)*aResult(k);
      aResult(i) = (b(i)-s)/A(i,i);
    }
    return aResult;
  }

  // invRegularized
  void invRegularized(CMatrix<float>& A, int aReg) {
    if (A.xSize() != A.ySize()) throw ENonquadraticMatrix(A.xSize(),A.ySize());
    CVector<float> eVals(A.xSize());
    CMatrix<float> eVecs(A.xSize(),A.ySize());
    PATransformation(A,eVals,eVecs);
    for (int i = 0 ; i < A.xSize(); i++)
      if (eVals(i) < aReg) eVals(i) = 1.0/aReg;
      else eVals(i) = 1.0/eVals(i);
    PABacktransformation(eVecs,eVals,A);
  }

  /********************************************
  Cholesky decomposition. Upper triang of x
  overwritten by decomposition. Nonsquare
  and non-definite matrices return false
  ********************************************/

  bool cholesky(CMatrix<double>& A) {
    //bool noerr = true; 
    int i, j, k;
    double g, h;
    if ( A.ySize() != A.xSize() ) {
      return false; // not square matrix
    }
    for ( j = 0 ; j < A.xSize(); j++ ) {
      g = A(j,j) ;
      for ( k = 0 ; k < j ; k++ )
        g -= A(j,k) * A(j,k) ;
      //std::cout << "g: " << g << std::endl;
      if ( g <= 0.0 ) {
	std::cerr << "NPDM: " << j << " " << g << std::endl;
	//noerr = false;
	//g = 0.000001; 
	return false; // not positive definite matrix
      }
      A(j,j) = g = sqrt( g ) ;
      for ( i = j + 1 ; i < A.xSize() ; i++ ) {
	h = A(i,j) ;
	for ( k = 0 ; k < j ; k++ )
	  h -= A(i,k) * A(j,k) ;
	A(i,j) = h / g ;
      } // for i
    } // for j

    // delete lower triang 
    for ( j = 1 ; j < A.ySize(); j++ ) 
      for ( i = 0 ; i < j ; i++ )
	A(i,j) = 0;
    return true;
  } // cholesky

  // Invert triangle matrix (upper triangle) 
  void invTriangleMatrix(CMatrix<double>& A) {
    for(int i=0; i<A.xSize(); ++i) {
      A(i,i) = 1.0/A(i,i);
      for(int j=i+1; j<A.xSize(); ++j) {
	double sum = 0.0;
	for(int k=i;k<j;++k) {
	  sum -= A(j,k)*A(k,i);
	}
	A(j,i) = sum/A(j,j);
      }
    } 
  }

  // Solve cholesky Ax = b => x = invL'*invL*b
  void solveCholesky(const CMatrix<double>& A, const CVector<double>& b, CVector<double>& x) {
    for(int i = 0; i<A.xSize(); ++i) {
      x[i] = 0.0;
      for(int j = 0; j<=i; ++j)
	x[i] += A(i,j)*b[j];
    }
    //std::cout << x << std::endl; 
   
    for(int i = 0; i<A.ySize(); ++i) {
      double sum = 0.0;
      for(int j = i; j<A.xSize(); ++j)
	sum += A(j,i)*x[j];
      x[i] = sum;
    }
    //std::cout << x << std::endl; 
    
  }
 
  
  // For EuclideanDistanceTransform
  void dt(CVector<float>& f, CVector<float>& d, int n) {
    d.setSize(n);
    int *v = new int[n];
    float *z = new float[n+1];
    int k = 0;
    v[0] = 0;
    z[0] = -10e20;
    z[1] = 10e20;
    for (int q = 1; q <= n-1; q++) {
      float s  = ((f[q]+q*q)-(f(v[k])+v[k]*v[k]))/(2*(q-v[k]));
      while (s <= z[k]) {
	k--;
	s  = ((f[q]+q*q)-(f(v[k])+v[k]*v[k]))/(2*(q-v[k]));
      }
      k++;
      v[k] = q;
      z[k] = s;
      z[k+1] = 10e20;
    }
    k = 0;
    for (int q = 0; q <= n-1; q++) {
      while (z[k+1] < q)
	k++;
      int help = q-v[k];
      d(q) = help*help + f(v[k]);
    }
    delete[] v;
    delete[] z;
  }

  // euclideanDistanceTransform
  void euclideanDistanceTransform(CMatrix<float>& aMatrix, CMatrix<bool>* aExclude) {
    int aXSize = aMatrix.xSize();
    int aYSize = aMatrix.ySize();
    CMatrix<float> aComputed(aXSize,aYSize,10e20);
    // Find zero level line
    if (aExclude != 0) for (int y = 0; y < aYSize-1; y++)
      for (int x = 0; x < aXSize-1; x++) {
	if (aMatrix(x,y)*aMatrix(x+1,y) < 0 && (*aExclude)(x,y) == false && (*aExclude)(x+1,y) == false) {
	  aComputed(x,y) = 0.0; aComputed(x+1,y) = 0;
	}
	if (aMatrix(x,y)*aMatrix(x,y+1) < 0 && (*aExclude)(x,y) == false && (*aExclude)(x,y+1) == false) {
	  aComputed(x,y) = 0.0; aComputed(x,y+1) = 0;
	}
      }
    else for (int y = 0; y < aYSize-1; y++)
      for (int x = 0; x < aXSize-1; x++) {
	if (aMatrix(x,y)*aMatrix(x+1,y) < 0) {
	  aComputed(x,y) = 0.0; aComputed(x+1,y) = 0;
	}
	if (aMatrix(x,y)*aMatrix(x,y+1) < 0) {
	  aComputed(x,y) = 0.0; aComputed(x,y+1) = 0;
	}
      }
    CVector<float> f(mmax(aXSize,aYSize));
    // Transform along columns
    for (int x = 0; x < aXSize; x++) {
      for (int y = 0; y < aYSize; y++)
	f(y) = aComputed(x,y);
      CVector<float> d;
      dt(f,d,aYSize);
      for (int y = 0; y < aYSize; y++)
	aComputed(x,y) = d(y);
    }
    // Transform along rows
    for (int y = 0; y < aYSize; y++) {
      int aOffset = y*aXSize;
      for (int x = 0; x < aXSize; x++)
	f(x) = aComputed.data()[x+aOffset];
      CVector<float> d;
      dt(f,d,aXSize);
      for (int x = 0; x < aXSize; x++)
	aComputed.data()[x+aOffset] = d(x);
    }
    int aSize = aMatrix.size();
    for (int i = 0; i < aSize; i++) {
      if (aMatrix.data()[i] >= 0) aMatrix.data()[i] = 1.0+sqrt(aComputed.data()[i]);
      else aMatrix.data()[i] = -1.0-sqrt(aComputed.data()[i]);
    }
  }
  
  void euclideanDistanceTransform(CTensor<float>& aTensor, CMatrix<bool>* aExclude) {
    int aXSize = aTensor.xSize();
    int aYSize = aTensor.ySize();
    for (int k = 0; k < aTensor.zSize(); k++) {
      CMatrix<float> aComputed(aXSize,aYSize,10e20);
      // Find zero level line
      if (aExclude != 0) for (int y = 0; y < aYSize-1; y++)
	for (int x = 0; x < aXSize-1; x++) {
	  if (aTensor(x,y,k)*aTensor(x+1,y,k) < 0 && (*aExclude)(x,y) == false && (*aExclude)(x+1,y) == false) {
	    aComputed(x,y) = 0.0; aComputed(x+1,y) = 0;
	  }
	  if (aTensor(x,y,k)*aTensor(x,y+1,k) < 0 && (*aExclude)(x,y) == false && (*aExclude)(x,y+1) == false) {
	    aComputed(x,y) = 0.0; aComputed(x,y+1) = 0;
	  }
	}
      else for (int y = 0; y < aYSize-1; y++)
	for (int x = 0; x < aXSize-1; x++) {
	  if (aTensor(x,y,k)*aTensor(x+1,y,k) < 0) {
	    aComputed(x,y) = 0.0; aComputed(x+1,y) = 0;
	  }
	  if (aTensor(x,y,k)*aTensor(x,y+1,k) < 0) {
	    aComputed(x,y) = 0.0; aComputed(x,y+1) = 0;
	  }
	}
      CVector<float> f(mmax(aXSize,aYSize));
      // Transform along columns
      for (int x = 0; x < aXSize; x++) {
	for (int y = 0; y < aYSize; y++)
	  f(y) = aComputed(x,y);
	CVector<float> d;
	dt(f,d,aYSize);
	for (int y = 0; y < aYSize; y++)
	  aComputed(x,y) = d(y);
      }
      // Transform along rows
      for (int y = 0; y < aYSize; y++) {
	int aOffset = y*aXSize;
	for (int x = 0; x < aXSize; x++)
	  f(x) = aComputed.data()[x+aOffset];
	CVector<float> d;
	dt(f,d,aXSize);
	for (int x = 0; x < aXSize; x++)
	  aComputed.data()[x+aOffset] = d(x);
      }
      int aSize = aTensor.size();
      int i2 = k*aSize;
      for (int i = 0; i < aSize; i++, i2++) {
	if (aTensor.data()[i2] >= 0) aTensor.data()[i2] = 1.0+sqrt(aComputed.data()[i]);
	else aTensor.data()[i2] = -1.0-sqrt(aComputed.data()[i]);
      }
    }
  }
  
  // Chamfer distance transformation
  bool chamferDistanceTransform(CMatrix<float>& aMatrix) {
    //CMatrix<int> test(aMatrix.xSize(),aMatrix.ySize());

    int d1 = 3; int d2 = 4;
    int maxV = INT_MAX-d2;
	//// liuyebin
	int size = aMatrix.size();
	int* v = (int*)malloc(size*(sizeof(int)));
  //  int v[aMatrix.size()];
    bool isEmpty = true;
    for(int i=0; i<aMatrix.size(); ++i) {
      if(aMatrix.data()[i]>0) {
	v[i] = 0;
	isEmpty = false;
      } else
	v[i] = maxV;
    }
    if(isEmpty) return false;
 
    //test = CMatrix<int>(v,aMatrix.xSize(),aMatrix.ySize());
    //std::cout << test << std::endl;

    int xMax = aMatrix.xSize()-1; 
    int off, off1;

    // left to right
    for(int x=1; x<aMatrix.xSize(); ++x) {
      int s = v[x-1]+d1;
      if(s < v[x]) v[x] = s;
    }
    
    for(int y=1; y<aMatrix.ySize(); ++y){
      off = (y-1)*aMatrix.xSize();
      off1 = y*aMatrix.xSize();
      int s = v[off]+d1; 
      if(s<v[off1]) v[off1] = s;
      s = v[off+1]+d2;
      if(s<v[off1]) v[off1] = s;
      
      for(int x=1; x<xMax; ++x) {
	++off;
	++off1;
	// if clauses speeds up the computation only when 
	// many of the values are zero else it slows down
	if(v[off1]>0) {
	  s = v[off-1]+d2; 
	  if(s<v[off1]) v[off1] = s;
	  s = v[off]+d1; 
	  if(s<v[off1]) v[off1] = s;
	  s = v[off+1]+d2;
	  if(s<v[off1]) v[off1] = s;
	  s = v[off1-1]+d1;
	  if(s<v[off1]) v[off1] = s;
	}
      }
      
      ++off;
      ++off1;
    
      s = v[off-1]+d2; 
      if(s<v[off1]) v[off1] = s;
      s = v[off]+d1; 
      if(s<v[off1]) v[off1] = s;
      s = v[off1-1]+d1;
      if(s<v[off1]) v[off1] = s;
     
    }

    // right to left
    off1 = (aMatrix.ySize()-1)*aMatrix.xSize()+xMax;
    for(int x=xMax-1; x>=0; --x) {
      --off1;
      int s = v[off1+1]+d1;
      if(s < v[off1]) v[off1] = s;
    }

    for(int y=aMatrix.ySize()-2;y>=0; --y){
      off = (y+1)*aMatrix.xSize()+xMax;
      off1 = y*aMatrix.xSize()+xMax;
      int s = v[off]+d1; 
      if(s<v[off1]) v[off1] = s;
      s = v[off-1]+d2;
      if(s<v[off1]) v[off1] = s;
   
      for(int x=xMax-1; x>0; --x) {
	--off;
	--off1;
	// if clauses speeds up the computation only when 
	// many of the values are zero else it slows down
	if(v[off1]>0) {
	  s = v[off-1]+d2; 
	  if(s<v[off1]) v[off1] = s;
	  s = v[off]+d1; 
	  if(s<v[off1]) v[off1] = s;
	  s = v[off+1]+d2;
	  if(s<v[off1]) v[off1] = s;
	  s = v[off1+1]+d1;
	  if(s<v[off1]) v[off1] = s;
	}
      }

      --off;
      --off1;
      s = v[off+1]+d2; 
      if(s<v[off1]) v[off1] = s;
      s = v[off]+d1; 
      if(s<v[off1]) v[off1] = s;
      s = v[off1+1]+d1;
      if(s<v[off1]) v[off1] = s;
     
    }

    for(int i=0; i<aMatrix.size(); ++i)
      aMatrix.data()[i] = (float)v[i]/3.0f;

	free(v);
    return true;
  }

  // restrictToBoundary
  void getBoundary(CMatrix<float>& aRegion, CMatrix<float>& aBorder) {
    for(int i = 1; i<aRegion.xSize()-1; ++i)
      for(int j = 1; j<aRegion.ySize()-1; ++j) 
	if( aRegion(i,j)>0 && ( aRegion(i-1,j)<=0 || aRegion(i,j-1)<=0 || aRegion(i+1,j)<=0 || aRegion(i,j+1)<=0) ) 
	  aBorder(i,j) = 1;
  }

  // restrictToBoundary skip image border
  void getBoundary(CMatrix<float>& aRegion, CMatrix<float>& aBorder, int b) {
    for(int i = b; i<aRegion.xSize()-b; ++i)
      for(int j = b; j<aRegion.ySize()-b; ++j) 
	if( aRegion(i,j)>0 && ( aRegion(i-1,j)<=0 || aRegion(i,j-1)<=0 || aRegion(i+1,j)<=0 || aRegion(i,j+1)<=0) ) 
	  aBorder(i,j) = 1;
  }
}

