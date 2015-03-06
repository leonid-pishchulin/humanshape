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

#ifndef UNSUP_LEARN
#define UNSUP_LEARN

#include "GaussVector.h"
#include "Mle.h"
#include <stdlib.h>
#include <search.h>
#include <list>

struct SComponent
{
	double importance;
	int index;
	static int componentCompare(const void * comp1, const void * comp2)
	{
		SComponent * c1 = (SComponent *)comp1;
		SComponent * c2 = (SComponent *)comp2;
		if(c1->importance > c2->importance) return -1;
		else if(c1->importance == c2->importance) return 0;
		else return 1;
	}
};

class UnsupervisedLearning
{
public:
	UnsupervisedLearning();
	UnsupervisedLearning(int dimension, int numberClasses, int numberSamples);
	~UnsupervisedLearning();

	void kMeansClustering(double * samples, double *& means);
	void nearestNeighbor(double * samples);

	void getClustering(double** samples, int** number);

	//Dimensionality reduction:
	//return reduced samples as well as reductionMatrix and eigenvalues:
    void performPCA(int inputDim, int outputDim, double * samples, double * outputSamples, double * reductionMatrix, 
		double * returnEigenvalues, double * returnMean);
    
	void performPCASpaceRestricted(int inputDim, int outputDim, double * samples, double * reductionMatrix, 
		double * returnEigenvalues, double * returnMean);

	//Same as above but using SVD instad of calculating the eigenvalues and eigenvectors of the covariance matrix.
	void performPCASpaceRestrictedSVD(int inputDim, int outputDim, double * samples, double * reductionMatrix, 
		double * returnEigenvalues, double * returnMean);

	//do a selective PCA for classification; choose the components according to how well they separate the data:
	//if outputDim == -1 then we choose all "well-separated" components where the distance between the means is 
	//more than twice the maximum variance
	//NOTE reductionMatrix and returnEigenvalues have full dimensionality (inputDim x inputDim)
	std::list<int> performClassificationPCA(int inputDim, int outputDim, int numSamplesClass1, int numSamplesClass2, double *
		inputSamples, double * outputSamples, double * reductionMatrix, double * returnEigenvalues, double * returnMean);

	//use the output from the above method to reduce the dimensionality of testing data:
	void reduceDimensionality(int numTestSamples, int inputDim, double * samples, double * outputSamples, 
		double * reductionMatrix, double * eigenvalues, double * mean, std::list<int> relevantComp);

	//Same as above but without PCA alignment:
	std::list<int> performClassificationReduction(int inputDim, int numSamplesClass1, int numSamplesClass2, 
		double * inputSamples, double * outputSamples);
	void reduceDimensionality(int numTestSamples, int inputDim, double * samples, double * outputSamples, 
		std::list<int> relevantComp);

private:
	int dimension;
	int numberClasses;
	int numberSamples;
	//inefficient in terms of storage, but simple and sufficient for small examples:
	double * samplesByClass;
	int * numberClass;
};

#endif //MLE
