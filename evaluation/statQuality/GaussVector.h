/**
    This file is part of the evaluation of the 3D human shape model as described in the paper:

    Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobalt and Bernt Schiele
    Building Statistical Shape Spaces for 3D Human Modeling
    ArXiv, March 2015

    Please cite the paper if you are using this code in your work.
    
    Contributors: Stefanie Wuhrer, Leonid Pishchulin.

    The code may be used free of charge for non-commercial and
    educational purposes, the only requirement is that this text is
    preserved within the derivative work. For any other purpose you
    must contact the authors for permission. This code may not be
    redistributed without permission from the authors.
*/

#ifndef GAUSS_VECTOR
#define GAUSS_VECTOR

#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <vector>
#include "matrix.h"
//Mex 

#include <math.h>
#include "mex.h"
#include "blas.h"
//#include "stddef.h"
#include "lapack.h"    
//#include  "blascompat32.h"

//includes for CLAPACK
// namespace clapack
// {
// 	extern "C"
//     {
//         
// 		#include "blaswrap.h"
// 		#include "f2c.h"
// 		extern int dgemm_(char *transa, char *transb, integer *m, integer *
// 			n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
// 			doublereal *b, integer *ldb, doublereal *beta, doublereal *c, 
// 			integer *ldc);
// 		extern int dgemv_(char *trans, integer *m, integer *n, doublereal *
// 			alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
// 			doublereal *beta, doublereal *y, integer *incy);
// 		extern int dsyevx_(char *jobz, char *range, char *uplo, integer *n, 
// 			doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
// 			il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
// 			doublereal *z__, integer *ldz, doublereal *work, integer *lwork, 
// 			integer *iwork, integer *ifail, integer *info);
// 		extern double dlamch_(char *cmach);
// 		extern int dgetrf_(integer *m, integer *n, doublereal *a, integer *
// 			lda, integer *ipiv, integer *info);
// 		extern int dgetri_(integer *n, doublereal *a, integer *lda, integer 
// 			*ipiv, doublereal *work, integer *lwork, integer *info);
// 		extern int dgesdd_(char *jobz, integer *m, integer *n, doublereal *
// 			a, integer *lda, doublereal *s, doublereal *u, integer *ldu, 
// 			doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
// 			integer *iwork, integer *info);
// 		extern int dgels_(char *trans, integer *m, integer *n, integer *
// 			nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
// 			doublereal *work, integer *lwork, integer *info);
// 	}
// }


class GaussVector
{
public:
	//constructor:
	GaussVector();
	GaussVector(int dimension, double * mean, double * covMat); 
	~GaussVector();

	//generate a random vector of a Gaussian distribution:
    //Monica-Changed qualification error    
    //void GaussVector::generateRandomVector(double *& result);
    void generateRandomVector(double *& result);
	//if returnUnits = true, unitVec returns unit vectors. Otherwise, you must provide unit vectors!!!
	void generateRandomVector(double *& result, bool returnUnits, double *& unitVec);	
	//generate numberSamples vectors of a Gaussian distribution. If returnUnits = true, unitVec returns unit vectors.
	//Otherwise, you must provide unit vectors!!!
	void generateRandomVectors(int numberSamples, double *& results, bool returnUnits, double *& unitVecs);

	//Method to sample on a surface of equal likelihood (an ellipsoid)
	void generateRandomVector(double distFromMean, double *& result);

	//methods to return computed quantities:
	int getDimension(){return dimension;}
	double * getMean(){return mean;}
	double * getCovarianceMat(){return covarianceMat;}
	double * getEigenvalues(){return eigenvalues;}
	double * getEigenvectors(){return eigenvectors;}
	double * getSqrtEigenvalueMatrix(){return sqrtEigenvalueMatrix;}
	double * getMultMat(){return multMatrix;}
	double getDeterminantOfCov();

	//method to compute and print numberSamples random vectors to a file (filename) that is compatible with Maple 7 
	//for plotting:
	void printSampleFile(int numberSamples, char * filename);
	void printSampleFile(int numberSamples, char * filename, double *& samples);
	//print only selected dimensions:
	void printSampleFile(int numberSamples, char * filename, int numDimensions, int * dimensions);
	void printSampleFile(int numberSamples, char * filename, double *& samples, int numDimensions, int * dimensions);

private:
	//compute the eigenvalues and eigenvectors of the covariance matrix using the CLAPACK library
	bool computeEigen(double * covMatrix);
	//compute a one-dimensional unit normal distribution:
	double getUnitNormalDistribution();
	//member variables:
	int dimension;
	double * mean; 
	double * covarianceMat;
	double * eigenvalues;
	double * eigenvectors;
	double * sqrtEigenvalueMatrix;
	double * multMatrix;
};

#endif // GAUSS_VECTOR
