/**
    This file is part of the evaluation of the 3D human shape model as described in the paper:

    Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobald and Bernt Schiele
    Building Statistical Shape Spaces for 3D Human Modeling
    ArXiv, February 2015

    Please cite the paper if you are using this code in your work.
    
    Contributors: Stefanie Wuhrer, Leonid Pishchulin.

    The code may be used free of charge for non-commercial and
    educational purposes, the only requirement is that this text is
    preserved within the derivative work. For any other purpose you
    must contact the authors for permission. This code may not be
    redistributed without permission from the authors.
*/

#include "GaussVector.h"

#if !defined(_WIN32)
//#define dgemm dgemm_
//#define dsyevx dsyevx_
//#define  dgemv dgemv_ 
//#define  dlamch dlamch_
#endif



GaussVector::GaussVector()
{
	dimension = -1;
}

GaussVector::GaussVector(int dimension, double * mean, double * covMat)
{
	this->dimension = dimension;
	this->mean = new double[dimension];
	for(int i = 0; i < dimension; i++)
		this->mean[i] = mean[i];
	this->covarianceMat = new double [dimension * dimension];
	for(int i = 0; i < dimension * dimension; i++)
		this->covarianceMat[i] = covMat[i];
	this->eigenvectors = new double[dimension * dimension];
	this->eigenvalues = new double[dimension];
	this->sqrtEigenvalueMatrix = new double[dimension * dimension];
	this->multMatrix = new double[dimension * dimension];
	computeEigen(covMat);
}

GaussVector::~GaussVector()
{
	if(dimension != -1)
	{
		delete [] mean;
		delete [] covarianceMat;
		delete [] eigenvectors;
		delete [] eigenvalues;
		delete [] sqrtEigenvalueMatrix;
		delete [] multMatrix;
	}
}

//compute the eigenvalues and eigenvectors of the covariance matrix using the CLAPACK library
bool GaussVector::computeEigen(double * covMatrix)
{
	

    //find eigenvectors and eigenvalues:
	//initialize the input values to the CLAPACK function dsyevx
	
    //Monica type change long int by  ptrdiff_t

    char jobz = 'V';
	char range = 'A';
	char uplo = 'L';
	ptrdiff_t nN = (ptrdiff_t) dimension;

	double vl = 0;
	ptrdiff_t il = 0;
	char s = 'S';
    
 // 	double abstol = 2 * clapack::dlamch_(&s);
 	double abstol = 2 * dlamch_ (&s);
    
    ptrdiff_t m;
	ptrdiff_t lwork = 8 * dimension;
	double *work = new double[lwork];
	ptrdiff_t *iwork = new ptrdiff_t[5 * dimension];
	ptrdiff_t *ifail = new ptrdiff_t[dimension];
	ptrdiff_t info;
	//test for validity of parameters:
	if((work == NULL) || (iwork == NULL) || (ifail == NULL))
		return false;

// 	clapack::dsyevx_(&jobz, &range, &uplo, &nN, covarianceMat, &nN, &vl, &vl, &il, &il,
// 		&abstol, &m, eigenvalues, eigenvectors, &nN, work, &lwork, iwork, ifail, &info);

    dsyevx_(&jobz, &range, &uplo, &nN, covarianceMat, &nN, &vl, &vl, &il, &il,
		&abstol, &m, eigenvalues, eigenvectors, &nN, work, &lwork, iwork, ifail, &info);

    
    
	//check for correctness:
	if(info != 0)
		return false;

	//create the sqrt eigenvector matrix (symmetric diagonal matrix):
	int offset;
	for(int i = 0; i < dimension; i++)
	{
		offset = i * dimension;
		for(int j = 0; j < dimension; j++)
		{
			//set diagonal elements to sqrt(eigenvalues[i])
			if(i == j)
				sqrtEigenvalueMatrix[offset + j] = sqrt(eigenvalues[i]);
			//otherwise, set the matrix to 0:
			else
				sqrtEigenvalueMatrix[offset + j] = 0;
		}
	}

	//multiply the eigenvector matrix with the sqrt eigenvalue matrix and store the result in multMatrix:
	//initialize the input values for the function dgemm:
	char transa = 'N';
	char transb = 'N';
	double alpha = 1.0;
	double beta = 0.0;

	//clapack::dgemm_(&transa, &transb, &nN, &nN, &nN, &alpha, eigenvectors, &nN, 
		//sqrtEigenvalueMatrix, &nN, &beta, multMatrix, &nN);
        
        dgemm_(&transa, &transb, &nN, &nN, &nN, &alpha, eigenvectors, &nN, 
		sqrtEigenvalueMatrix, &nN, &beta, multMatrix, &nN);
        
	//free the allocated space:
	delete [] work;
	delete [] iwork;
	delete [] ifail;

	//Repopulate the covariance matrix correctly:
	for(int i = 0; i < dimension * dimension; i++)
		covarianceMat[i] = covMatrix[i];

	return true;
}

void GaussVector::generateRandomVector(double distFromMean, double *& result)
{
	int i, j;
	double parameter;

	double * unitRandomVector = new double[dimension];

	for(i = 0; i < dimension; i++) unitRandomVector[i] = distFromMean * sqrt(covarianceMat[i*dimension+i]);

	//Sample a point on an axes-aligned ellipsoid of the desired dimensions centered at the origin
	for(i = 0; i < dimension-1; i++) 
	{
		parameter = ((double)rand()/((double)(RAND_MAX)+(double)(1)))*((double)(2.0*3.14159));
		for(j = 0; j < i+1; j++) unitRandomVector[j] *= sin(parameter);
		unitRandomVector[i+1] *= cos(parameter);
	}

	//Transform to non-axes aligned:
	char trans = 'T';
	ptrdiff_t  nN = (ptrdiff_t ) dimension;
	double alpha = 1.0;
	ptrdiff_t incx = 1;
	double beta = 0.0;

// 	clapack::dgemv_(&trans, &nN, &nN, &alpha, eigenvectors, &nN, unitRandomVector, &incx, 
// 		&beta, result, &incx);
    
    dgemv_(&trans, &nN, &nN, &alpha, eigenvectors, &nN, unitRandomVector, &incx,&beta, result, &incx);
    

	//Then, adjust the mean by adding the mean vector to the result!
	for(int i = 0; i < dimension; i++)
		result[i] = result[i] + mean[i];

	delete [] unitRandomVector;
}

//generate a random vector of a Gaussian distribution:
void GaussVector::generateRandomVector(double *& result)
{
	//obtain a random vector with mean 0 (vector 0) and covariance matrix I
	double * unitRandomVector = new double[dimension];
	for(int i = 0; i < dimension; i++)
		unitRandomVector[i] = getUnitNormalDistribution();

	//compute the true Gaussian vector by multiplying the unit random vector with 
	//multMatrix = eigenvectors * sqrtEigenvalueMatrix:

	//initialize the input values for the function dgemv:
	char trans = 'N';
	ptrdiff_t nN = (ptrdiff_t) dimension;
	double alpha = 1.0;
	ptrdiff_t incx = 1;
	double beta = 0.0;

//	clapack::dgemv_(&trans, &nN, &nN, &alpha, multMatrix, &nN, unitRandomVector, &incx, 
//		&beta, result, &incx);
    dgemv_(&trans, &nN, &nN, &alpha, multMatrix, &nN, unitRandomVector, &incx, 
		&beta, result, &incx);

    
	//Then, adjust the mean by adding the mean vector to the result!
	for(int i = 0; i < dimension; i++)
		result[i] = result[i] + mean[i];

	//free the allocated memory:
	delete [] unitRandomVector;
}

//generate a random vector of a Gaussian distribution:
void GaussVector::generateRandomVector(double *& result, bool returnUnits, double *& unitVec)
{
	int i;
	//obtain a random vector with mean 0 (vector 0) and covariance matrix I
	if(returnUnits)
	{
		for(i = 0; i < dimension; i++)
			unitVec[i] = getUnitNormalDistribution();
	}

	//compute the true Gaussian vector by multiplying the unit random vector with 
	//multMatrix = eigenvectors * sqrtEigenvalueMatrix:

	//initialize the input values for the function dgemv:
	char trans = 'N';
	ptrdiff_t nN = (ptrdiff_t) dimension;
	double alpha = 1.0;
	ptrdiff_t incx = 1;
	double beta = 0.0;

	//clapack::dgemv_(&trans, &nN, &nN, &alpha, multMatrix, &nN, unitVec, &incx, 
//	&beta, result, &incx);

dgemv(&trans, &nN, &nN, &alpha, multMatrix, &nN, unitVec, &incx, 
		&beta, result, &incx);

    
//Then, adjust the mean by adding the mean vector to the result!
	for(i = 0; i < dimension; i++)
		result[i] = result[i] + mean[i];
}

//generate numberSamples vectors of a Gaussian distribution:
void GaussVector::generateRandomVectors(int numberSamples, double *& results, bool returnUnits, double *& unitVecs)
{
	int offset;
	double * resultVec = new double [dimension];
	double * unitVec = new double [dimension];
	//results and unitVecs need to be of dimension numberSamples x dimension:
	for(int i = 0; i < numberSamples; i++)
	{
		offset = i*dimension;
		if(!returnUnits)
		{
			for(int j = 0; j < dimension; j++)
				unitVec[j] = unitVecs[offset+j];
		}
		generateRandomVector(resultVec, returnUnits, unitVec);
		for(int j = 0; j < dimension; j++)
		{
			results[offset+j] = resultVec[j];
			if(returnUnits)
				unitVecs[offset+j] = unitVec[j];
		}
	}
	delete [] resultVec;
	delete [] unitVec;
}

//print numberSamples random vectors to a file (filename) that is compatible with Maple 7 for plotting
void GaussVector::printSampleFile(int numberSamples, char * filename)
{
	FILE * fp = fopen(filename,"w");
	double * result = new double [dimension];
	fprintf(fp, "[");
	for(int i = 0; i < numberSamples; i++)
	{
		generateRandomVector(result);
		fprintf(fp, "[");
		for(int j = 0; j < dimension; j++)
		{
			fprintf(fp, "%f", result[j]);
			if(j != dimension - 1)
				fprintf(fp, ",");
			else
				fprintf(fp, "]");
		}
		if(i != numberSamples -1)
			fprintf(fp, ", ");
		else
			fprintf(fp, "]");
	}
	fclose(fp);
	delete [] result;
}

void GaussVector::printSampleFile(int numberSamples, char * filename, double *& samples)
{
	FILE * fp = fopen(filename,"w");
	fprintf(fp, "[");
	for(int i = 0; i < numberSamples; i++)
	{
		fprintf(fp, "[");
		for(int j = 0; j < dimension; j++)
		{
			fprintf(fp, "%f", samples[i * dimension + j]);
			if(j != dimension - 1)
				fprintf(fp, ",");
			else
				fprintf(fp, "]");
		}
		if(i != numberSamples -1)
			fprintf(fp, ", ");
		else
			fprintf(fp, "]");
	}
	fclose(fp);
}

//print only selected dimensions:
void GaussVector::printSampleFile(int numberSamples, char * filename, int numDimensions, int * dimensions)
{
	FILE * fp = fopen(filename,"w");
	double * result = new double [dimension];
	fprintf(fp, "[");
	for(int i = 0; i < numberSamples; i++)
	{
		generateRandomVector(result);
		fprintf(fp, "[");
		for(int j = 0; j < numDimensions; j++)
		{
			fprintf(fp, "%f", result[dimensions[j]]);
			if(j != numDimensions - 1)
				fprintf(fp, ",");
			else
				fprintf(fp, "]");
		}
		if(i != numberSamples -1)
			fprintf(fp, ", ");
		else
			fprintf(fp, "]");
	}
	fclose(fp);
	delete [] result;
}

void GaussVector::printSampleFile(int numberSamples, char * filename, double *& samples, int numDimensions, int * dimensions)
{
	FILE * fp = fopen(filename,"w");
	fprintf(fp, "[");
	for(int i = 0; i < numberSamples; i++)
	{
		fprintf(fp, "[");
		for(int j = 0; j < numDimensions; j++)
		{
			fprintf(fp, "%f", samples[i * dimension + dimensions[j]]);
			if(j != numDimensions - 1)
				fprintf(fp, ",");
			else
				fprintf(fp, "]");
		}
		if(i != numberSamples -1)
			fprintf(fp, ", ");
		else
			fprintf(fp, "]");
	}
	fclose(fp);
}

//compute a one-dimensional unit normal distribution:
double GaussVector::getUnitNormalDistribution()
{
	double result = 0;
	//By central limit theorem use 48 = 12 * 4 trials of uniform random variable:
	for(int i = 0; i < 48; i++)
	{
		//add a random number in (0, 1]
		result += ((double)rand()/((double)(RAND_MAX)+(double)(1)));
	}
	//the variable result is now of a Gaussian distribution (24, 4).
	result = result / 2;
	//the variable result is now of a Gaussian distribution (12, 1).
	result = result - 12;
	//the variable result is now of a Gaussian distribution (0, 1), as wanted.
	return result;
}

double GaussVector::getDeterminantOfCov()
{
	double determinant = 1;
	for(int i = 0; i < dimension; i++) 
		determinant = determinant * eigenvalues[i];
	return determinant;
}


