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

#ifndef GPA
#define GPA

#define ROTATE_THRESHOLD 0.00001
#define MAX_WHILE 100

#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <vector>



//Mex 
 extern "C" {
#include <math.h>
#include "mex.h"
#include "blas.h"
#include "stddef.h"
#include "lapack.h"    

}

// includes for CLAPACK
// namespace clapack
// {
//         extern "C"
//         {
//                 #include "blaswrap.h"
//                 #include "f2c.h"
//                 multiplications:
//                 int dsymm_(char *side, char *uplo, integer *m, integer *n, 
//                         doublereal *alpha, doublereal *a, integer *lda, doublereal *b, 
//                         integer *ldb, doublereal *beta, doublereal *c__, integer *ldc);
//                 int dgemm_(char *transa, char *transb, integer *m, integer *
//                         n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
//                         doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, 
//                         integer *ldc);
//                 int dgemv_(char *trans, integer *m, integer *n, doublereal *
//                         alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
//                         doublereal *beta, doublereal *y, integer *incy);
//                 compute eigenvectors and eigenvalues:
//                 int dsyevx_(char *jobz, char *range, char *uplo, integer *n, 
//                         doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
//                         il, integer *iu, doublereal *abstol, integer *m, doublereal *w, 
//                         doublereal *z__, integer *ldz, doublereal *work, integer *lwork, 
//                         integer *iwork, integer *ifail, integer *info);
//                 double dlamch_(char *cmach);
//                 compute inverse of a matrix:
//                 int dgetrf_(integer *m, integer *n, doublereal *a, integer *
//                         lda, integer *ipiv, integer *info);
//                 int dgetri_(integer *n, doublereal *a, integer *lda, integer 
//                         *ipiv, doublereal *work, integer *lwork, integer *info);
//                 int dgesv_(integer *n, integer *nrhs, doublereal *a, integer 
//                         *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);
//                 int dgesvx_(char *fact, char *trans, integer *n, integer *
//                         nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
//                         integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
//                         doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
//                         rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
//                         iwork, integer *info);
//                 int dgels_(char *trans, integer *m, integer *n, integer *
//                         nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
//                         doublereal *work, integer *lwork, integer *info);
//                 int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
//                         doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
//                         ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
//                         integer *info);
//         }
// }

using namespace std;

//Set of numberVertices vertices containing typically x,y, and z coordinates each
struct Shape
{
        int numberVertices;
        vector < vector < double > > vertexCoordinates;
        //vector < vector < double > > vertexCoordinates;
};


class GeneralizedProcrustes {

public:
        GeneralizedProcrustes();

        //As input the method expects a vector of shapes.
        GeneralizedProcrustes(vector < Shape >& shapes, int numberOfShapes, int minNumNodes);
        ~GeneralizedProcrustes();

        void computeGPAlignment( vector < Shape >& outShapes);
        //compute the average of the aligned shapes:
        void computeAverage(Shape & average);
        //evaluate the error energy:
        double evaluateEnergy();

private:
        int numberOfShapes;
        long int minNumNodes;
        vector < Shape > shapes;
        double ** partialAverages;

        void computePartialAverages();
        void centerAllShapes();
        bool computeRotationAlignment(Shape & cont, double * partialAverage);
};

#endif
