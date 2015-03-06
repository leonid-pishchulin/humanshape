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

#include "Mle.h"

Mle::Mle()
{
}

Mle::~Mle()
{
}

void Mle::learnExponentialDistribution(int dimension, int sampleSize, double *& samples, double *& learnedLambda)
{
	int offset;
	for(int j = 0; j < dimension; j++)
		learnedLambda[j] = 0;
	for(int i = 0; i < sampleSize; i++)
	{
		offset = i*dimension;
		for(int j = 0; j < dimension; j++)
			learnedLambda[j] += samples[offset+j];
	}
	for(int j = 0; j < dimension; j++)
		learnedLambda[j] = sampleSize / learnedLambda[j];
}

void Mle::learnGaussDistribution(int dimension, int sampleSize, double * &samples, double *& learnedMean, 
		double *& learnedCovariance)      
        
{
  	//first, we learn the parameter mean:
	int offset;
   
    for(int j = 0; j < dimension; j++)
		learnedMean[j] = 0;
    
	for(int i = 0; i < sampleSize; i++)
	{
		offset = i*dimension;
		for(int j = 0; j < dimension; j++)
			learnedMean[j] +=samples[offset+j];
	}
	for(int j = 0; j < dimension; j++)
        learnedMean[j] = learnedMean[j]/ sampleSize;
 
    //second, we learn the parameter covariance:
	for(int j = 0; j < dimension * dimension; j++)
		learnedCovariance[j] = 0;
 
    double * differenceVec = new double [dimension];
	for(int i = 0; i < sampleSize; i++)
	{
		offset = i*dimension;
		for(int j = 0; j < dimension; j++)
			differenceVec[j] = samples[offset+j] - learnedMean[j];
		for(int j = 0; j < dimension; j++)
		{
			for(int k = 0; k < dimension; k++)
				learnedCovariance[j*dimension + k] += differenceVec[j] * differenceVec[k];
		}
	}
	for(int j = 0; j < dimension * dimension; j++)
		learnedCovariance[j] = learnedCovariance[j]/ sampleSize;
	delete [] differenceVec;
    
    
}
