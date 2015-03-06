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

#include "UnsupervisedLearning.h"
#include <iostream>
#include <exception>
UnsupervisedLearning::UnsupervisedLearning()
{
	samplesByClass = NULL;
}

UnsupervisedLearning::UnsupervisedLearning(int dimension, int numberClasses, int numberSamples)
{
	this->dimension = dimension;
	this->numberClasses = numberClasses;
	this->numberSamples = numberSamples;
	samplesByClass = new double[numberSamples * numberClasses * dimension];
	numberClass = new int[numberClasses];
}

UnsupervisedLearning::~UnsupervisedLearning()
{
	if(samplesByClass != NULL)
	{
		delete [] samplesByClass;
		delete [] numberClass;
		samplesByClass = NULL;
	}
}

void UnsupervisedLearning::getClustering(double** samples, int** number)
{
	*samples = samplesByClass;
	*number = numberClass;
}

void UnsupervisedLearning::kMeansClustering(double * samples, double *& means)
{
	int offset1, offset2, offset3, minDistIdx, i, j, k;
	double minDist;
	double *distance = new double[numberClasses];
	bool meanChanged = true;
	double * learnedCovariance = new double [dimension * dimension];
	double * oldMean = new double[numberClasses * dimension];
	Mle mle;

	while(meanChanged)
	{
		meanChanged = false;
		for(k = 0; k < numberClasses; k++)
			numberClass[k] = 0;
		//classify all of the samples:
		for(i = 0; i < numberSamples; i++)
		{
			offset1 = i * dimension;
			for(k = 0; k < numberClasses; k++)
			{
				distance[k] = 0;
				offset2 = k * dimension;
				//compute the distance of the i-th sample to the k-th mean:
				for(j = 0; j < dimension; j++)
					distance[k] += pow((means[offset2+j] - samples[offset1+j]),2);
				if(k == 0 || distance[k] < minDist)
				{
					minDist = distance[k];
					minDistIdx = k;
				}
			}
			//store the classified sample:
			offset3 = minDistIdx * numberSamples * dimension + numberClass[minDistIdx] * dimension;
			for(j = 0; j < dimension; j++)
				samplesByClass[offset3 + j] = samples[offset1 + j];
			numberClass[minDistIdx]++;
		}

		//store the old means and compute the new means:
		for(k = 0; k < numberClasses * dimension; k++)
			oldMean[k] = means[k];
		for(k = 0; k < numberClasses; k++)
		{
			double * samplePointer = &samplesByClass[k * numberSamples * dimension];
			double * meanPointer = &means[k*dimension];
			mle.learnGaussDistribution(dimension, numberClass[k], samplePointer, meanPointer, learnedCovariance);
		}
		//test if the means changed. If yes, set meanChanged = true:
		for(k = 0; k < numberClasses * dimension; k++)
		{
			if(fabs(oldMean[k] - means[k]) > 0.001)
			{
				meanChanged = true;
				break;
			}
		}
	}

	delete [] distance;
	delete [] learnedCovariance;
	delete [] oldMean;
}

void UnsupervisedLearning::nearestNeighbor(double * samples)
{
	int i, j, k, numberClusters, idx1, idx2, offset;
	double distance, minDist;
	int * clusterId = new int[numberSamples];
	numberClusters = numberSamples;
	for(i = 0; i < numberSamples; i++)
		clusterId[i] = i;

	while(numberClusters != numberClasses)
	{
		minDist = -1;
		idx1 = -1;
		idx2 = -1;
		for(i = 0; i < numberSamples; i++)
		{
			for(j = i+1; j < numberSamples; j++)
			{
				//compute the distance between the two samples:
				distance = 0;
				for(k = 0; k < dimension; k++)
					distance += pow((samples[i*dimension+k] - samples[j*dimension+k]), 2);
				if((minDist==-1 && clusterId[i]!=clusterId[j]) || (distance < minDist && clusterId[i]!=clusterId[j]))
				{
					minDist = distance;
					idx1 = std::min(clusterId[i], clusterId[j]);
					idx2 = std::max(clusterId[i], clusterId[j]) ;
				}
			}
		}
		//merge the clusters idx1 and idx2:
		for(i = 0; i < numberSamples; i++)
		{
			if(clusterId[i] == idx2)
				clusterId[i] = idx1;
		}
		numberClusters--;
	}

	//store the clusters in samplesByClass
	int * classIds = new int[numberClasses];
	int numberViewed = 0;
	for(i = 0; i < numberSamples; i++)
	{
		offset = -1;
		for(j = 0; j < numberViewed; j++)
		{
			if(classIds[j] == clusterId[i])
				offset = j * numberSamples * dimension + numberClass[j] * dimension;
		}
		if(offset == -1)
		{
			offset = numberViewed * numberSamples * dimension + numberClass[numberViewed] * dimension;
			classIds[numberViewed] = clusterId[i];
			numberViewed++;
		}
		for(j = 0; j < dimension; j++)
			samplesByClass[offset + j] = samples[i*dimension + j];
		numberClass[clusterId[i]]++;
	}
	
	delete [] clusterId;
}

//return reduces samples. Also scale the samples by their eigenvalue:
void UnsupervisedLearning::performPCA(int inputDim, int outputDim, double * samples, double * outputSamples, 
									  double * aMatrix, double * returnEigenvalues, double * returnMean)
{

	int i, j;
	Mle mle;
	double * learnedMean = new double[inputDim];
    
    double * learnedCovariance = new double[inputDim * inputDim];// breaks here 

    mle.learnGaussDistribution(inputDim, numberSamples, samples, learnedMean, learnedCovariance);
	GaussVector vec(inputDim, learnedMean, learnedCovariance);
   	
//get the eigenvalues and eigenvectors:
	double * eigenvalues = vec.getEigenvalues();
	double * eigenvectors = vec.getEigenvectors();
    
    //Monica
    double *reductionMatrix= new double [inputDim*outputDim];
    
    ;
	//take the last outputDim eigenvectors:
	for(i = 0; i < outputDim; i++)
	{
		for(j = 0; j < inputDim; j++)
			reductionMatrix[j*outputDim + i] = eigenvectors[(inputDim-i-1)*inputDim + j];
	}
     
     
     
	//reduce the dimensionality of each sample:
	char trans = 'N';
	ptrdiff_t nN = (ptrdiff_t) inputDim;
	ptrdiff_t nM = ( ptrdiff_t) outputDim;
	double alpha = 1.0;
	ptrdiff_t incx = 1;
	double beta = 0.0;
	double * sample = new double[inputDim];
	double * result = new double[outputDim];

    try{
        for(i = 0; i < numberSamples; i++)
        {
            for(j = 0; j < inputDim; j++)
                sample[j] = samples[i*inputDim + j] - learnedMean[j];
            //clapack::dgemv_(&trans, &nM, &nN, &alpha, reductionMatrix, &nM, sample, &incx, 
            //	&beta, result, &incx);

            dgemv(&trans, &nM, &nN, &alpha, reductionMatrix, &nM, sample, &incx, 
                &beta, result, &incx);

            for(j = 0; j < outputDim; j++)
            {
                outputSamples[i*outputDim + j] = result[j];
                //Note: before, the order of the eigenvalues was ASCENDING. For return, it is DESCENDING.
                if(eigenvalues[inputDim - j - 1] < 0) eigenvalues[inputDim - j - 1] = - eigenvalues[inputDim - j - 1];
                //ALREADY NORMALIZED BY CLAPACK => DO NOT NORMALIZE AGAIN!
                returnEigenvalues[j] = eigenvalues[inputDim - j - 1];
            }
        }
    
    }catch(std::exception& e){printf(e.what());}
     for(i = 0; i < inputDim; i++)
    returnMean[i] = learnedMean[i];

    
    
    
	delete [] learnedMean;
	delete [] learnedCovariance;
	delete [] sample;
//	delete [] result;
}

void  UnsupervisedLearning::performPCASpaceRestricted(int inputDim, int outputDim, double * samples, double * reductionMatrix, 
		double * returnEigenvalues, double * returnMean)
{
	int i, j;
	Mle mle;
	double * learnedMean = new double[inputDim];
	double * learnedCovariance = new double[inputDim * inputDim];
	mle.learnGaussDistribution(inputDim, numberSamples, samples, learnedMean, learnedCovariance);

	//Compute the eigen information directly to save space:
	//double * eigenvectors = new double[inputDim * inputDim];
	double * eigenvectors = new double[inputDim * outputDim];
    double * eigenvalues = new double[inputDim];

	char jobz = 'V';
	char range = 'A';
	char uplo = 'L';
	ptrdiff_t nN = (ptrdiff_t) inputDim;

	double vl = 0;
	ptrdiff_t il = 0;
	char s = 'S';
	double abstol = 2 * dlamch(&s);
    //double abstol = 2 * clapack::dlamch_(&s);
	ptrdiff_t m;
	ptrdiff_t lwork = 8 * inputDim;
	double *work = new double[lwork];
	ptrdiff_t *iwork = new ptrdiff_t[5 * inputDim];
	ptrdiff_t *ifail = new ptrdiff_t[inputDim];
	ptrdiff_t info;
	//test for validity of parameters:
	if((work == NULL) || (iwork == NULL) || (ifail == NULL))
		return;

	//clapack::dsyevx_(&jobz, &range, &uplo, &nN, learnedCovariance, &nN, &vl, &vl, &il, &il,
	//	&abstol, &m, eigenvalues, eigenvectors, &nN, work, &lwork, iwork, ifail, &info);
    
   dsyevx(&jobz, &range, &uplo, &nN, learnedCovariance, &nN, &vl, &vl, &il, &il,
		&abstol, &m, eigenvalues, eigenvectors, &nN, work, &lwork, iwork, ifail, &info);

	//check for correctness:
	if(info != 0)
		return;

	//take the last outputDim eigenvectors:
	for(i = 0; i < outputDim; i++)
	{
		for(j = 0; j < inputDim; j++)
			reductionMatrix[j*outputDim + i] = eigenvectors[(inputDim-i-1)*inputDim + j];
	}

	for(j = 0; j < outputDim; j++)
	{
		//Note: before, the order of the eigenvalues was ASCENDING. For return, it is DESCENDING.
		if(eigenvalues[inputDim - j - 1] < 0) eigenvalues[inputDim - j - 1] = - eigenvalues[inputDim - j - 1];
		returnEigenvalues[j] = eigenvalues[inputDim - j - 1];
	}

	for(i = 0; i < inputDim; i++)
		returnMean[i] = learnedMean[i];

	delete [] learnedMean;
	delete [] learnedCovariance;
	delete [] eigenvalues;
	delete [] eigenvectors;
}

void  UnsupervisedLearning::performPCASpaceRestrictedSVD(int inputDim, int outputDim, double * samples, double * reductionMatrix, 
		double * returnEigenvalues, double * returnMean)
{
	double * learnedMean = new double[inputDim];

	for(int i = 0; i < inputDim; ++i)
		learnedMean[i] = 0;

 
for(int i = 0; i < numberSamples; ++i)
	{
		int offset = i*inputDim;
		for(int j = 0; j < inputDim; ++j)
			learnedMean[j] += samples[offset+j];
	}
 	

	for(int i = 0; i < inputDim; ++i)
		learnedMean[i] /= static_cast<double>(numberSamples);
 
	const double denominator = sqrt(static_cast<double>(numberSamples));

	double * centeredSample = new double[inputDim*numberSamples];
	for(int i = 0; i < numberSamples; ++i)
	{
		int offset = i*inputDim;
		for(int j = 0; j < inputDim; ++j)
		{        
    
			centeredSample[j+offset] = (samples[j+offset] - learnedMean[j])/denominator;
		}
	}

	char jobz = 'S';
	ptrdiff_t m = inputDim;
	ptrdiff_t n = numberSamples;

	ptrdiff_t  lda = inputDim;

	double * s = new double [std::min(inputDim,numberSamples)];
	double * u = new double[inputDim*inputDim];
    //double * u = new double[inputDim*numberSamples];
	ptrdiff_t ldu = inputDim;
	double * vt = new double[numberSamples*numberSamples];

	ptrdiff_t ldvt = numberSamples;

        
	ptrdiff_t lwork = 3*std::min(m,n)*std::min(m,n)+std::max(std::max(m,n),4*std::min(m,n)*std::min(m,n)+4*std::min(m,n));
	double * work = new double[lwork];
	
	ptrdiff_t * iwork = new ptrdiff_t[8*std::min(m,n)];
	ptrdiff_t info = 0;

//if JOBZ = 'S', U contains the first min(M,N) columns of U
//          (the left singular vectors, stored columnwise)

    dgesdd(&jobz, &m, &n, centeredSample, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info);

	if(info != 0)
	{
		delete [] learnedMean;
		delete [] centeredSample;
		delete [] s;
		delete [] u;
		delete [] vt;
		delete [] iwork;
		return;
	}


    for(int i = 0; i < outputDim; ++i)
    {
        returnEigenvalues[i] = s[i]*s[i];
        
        const int offset = i*inputDim;
        for(int j = 0; j < inputDim; ++j){
            
            
   	 		 reductionMatrix[j*outputDim+i] = u[offset+j];
           //  std::cout<<"vt[offset+j]"<<vt[offset+j]<<std::endl;
        
        }
	}
    
    
    
   
	for(int i = 0; i < inputDim; ++i)
        returnMean[i] = learnedMean[i];
       

	//std::cout<<"first"<<reductionMatrix[0]<<std::endl;
    
	delete [] learnedMean;
	delete [] centeredSample;
	delete [] s;
	delete [] u;
	delete [] vt;
	delete [] iwork;
}

std::list<int> UnsupervisedLearning::performClassificationPCA(int inputDim, int outputDim, int numSamplesClass1, int numSamplesClass2, double *
		inputSamples, double * outputSamples, double * reductionMatrix, double * returnEigenvalues, double * returnMean)
{
	int i, j, offset, componentCounter, stopIndex;
	double key;
	numberSamples = numSamplesClass1 + numSamplesClass2;
	std::list<int>returnList;
	SComponent * helpList = new SComponent[inputDim];
	Mle mle;
	double * samplesClass1 = new double[inputDim * numSamplesClass1];
	double * samplesClass2 = new double[inputDim * numSamplesClass2];
	double * learnedMean1 = new double[inputDim];
	double * learnedMean2 = new double[inputDim];
	double * learnedCovariance1 = new double[inputDim * inputDim];
	double * learnedCovariance2 = new double[inputDim * inputDim];
	double * outputSamplesPCA = new double[inputDim * numberSamples];

	//perform PCA on the pooled data set:
	performPCA(inputDim, inputDim, inputSamples, outputSamplesPCA, reductionMatrix, returnEigenvalues, returnMean);

	//learn the means and covariance matrices for the two aligned data sets separately:
	for(i = 0; i < numSamplesClass1; i++)
	{
		offset = i*inputDim;
		for(j = 0; j < inputDim; j++)
		{
			samplesClass1[offset+j] = outputSamplesPCA[offset+j];
		}
	}
	mle.learnGaussDistribution(inputDim, numSamplesClass1, samplesClass1, learnedMean1, learnedCovariance1);
	for(i = 0; i < numSamplesClass2; i++)
	{
		offset = (i+numSamplesClass1)*inputDim;
		for(j = 0; j < inputDim; j++)
			samplesClass2[i*inputDim+j] = outputSamplesPCA[offset+j];
	}
	mle.learnGaussDistribution(inputDim, numSamplesClass2, samplesClass2, learnedMean2, learnedCovariance2);

	//keep the components with larger mean distance than twice the covariance in sorted order:
	componentCounter = 0;
	for(i = 0; i < inputDim; i++)
	{
		key = fabs(learnedMean1[i] - learnedMean2[i]) - 2*std::max(sqrt(learnedCovariance1[i*inputDim+i]),
			sqrt(learnedCovariance2[i*inputDim+i]));
		if(key > 0)
		{
			helpList[componentCounter].importance = key;
			helpList[componentCounter].index = i;
			componentCounter++;
		}
	}
	qsort(helpList, inputDim, sizeof(SComponent), SComponent::componentCompare);
	//copy the indices into the output array:
	if(outputDim == -1) stopIndex = componentCounter;
	else stopIndex =outputDim;
	for(i = 0; i < stopIndex; i++)
	{
		returnList.push_back(helpList[i].index);
	};
	//populate the outputSample matrix:
	for(i = 0; i < numberSamples; i++)
	{
		std::list<int>::iterator it = returnList.begin();
		for(j = 0; j < stopIndex; j++)
		{
			outputSamples[i*stopIndex+j] = outputSamplesPCA[i*inputDim + *it];
			it++;
		}
	}

	delete [] samplesClass1;
	delete [] samplesClass2;
	delete [] learnedMean1;
	delete [] learnedMean2;
	delete [] learnedCovariance1;
	delete [] learnedCovariance2;
	delete [] helpList;
	delete [] outputSamplesPCA;

	return returnList;
}

void UnsupervisedLearning::reduceDimensionality(int numTestSamples, int inputDim, double * samples, double * outputSamples, 
		double * reductionMatrix, double * eigenvalues, double * mean, std::list<int> relevantComp)
{
	int i, j, count;
	char trans = 'T';
	 ptrdiff_t nN = ( ptrdiff_t) inputDim;
	double alpha = 1.0;
	ptrdiff_t incx = 1;
	double beta = 0.0;
	double * unalignedPoint = new double[inputDim];
	double * alignedPoint = new double[inputDim];
	for(i = 0; i < numTestSamples; i++)
	{
		//align the point along the learned principal axes:
		for(j = 0; j < inputDim; j++)
			unalignedPoint[j] = samples[i*inputDim+j] - mean[j];
		dgemv(&trans, &nN, &nN, &alpha, reductionMatrix, &nN, unalignedPoint, &incx, 
			&beta, alignedPoint, &incx);
		
		//finally store the relevant dimensions in outputSamples:
		std::list<int>::iterator it;
		count = 0;
		for(it = relevantComp.begin(); it != relevantComp.end(); it++, count++)
			outputSamples[i*relevantComp.size()+count] = alignedPoint[*it];
	}
	delete [] unalignedPoint;
	delete [] alignedPoint;
}

std::list<int> UnsupervisedLearning::performClassificationReduction(int inputDim, int numSamplesClass1, int numSamplesClass2, 
		double * inputSamples, double * outputSamples)
{
	int i, j, offset, componentCounter;
	double key;
	numberSamples = numSamplesClass1 + numSamplesClass2;
	std::list<int>returnList;
	SComponent * helpList = new SComponent[inputDim];
	Mle mle;
	double * samplesClass1 = new double[inputDim * numSamplesClass1];
	double * samplesClass2 = new double[inputDim * numSamplesClass2];
	double * learnedMean1 = new double[inputDim];
	double * learnedMean2 = new double[inputDim];
	double * learnedCovariance1 = new double[inputDim * inputDim];
	double * learnedCovariance2 = new double[inputDim * inputDim];

	//learn the means and covariance matrices for the two aligned data sets separately:
	for(i = 0; i < numSamplesClass1; i++)
	{
		offset = i*inputDim;
		for(j = 0; j < inputDim; j++)
		{
			samplesClass1[offset+j] = inputSamples[offset+j];
		}
	}
	mle.learnGaussDistribution(inputDim, numSamplesClass1, samplesClass1, learnedMean1, learnedCovariance1);
	for(i = 0; i < numSamplesClass2; i++)
	{
		offset = (i+numSamplesClass1)*inputDim;
		for(j = 0; j < inputDim; j++)
			samplesClass2[i*inputDim+j] = inputSamples[offset+j];
	}
	mle.learnGaussDistribution(inputDim, numSamplesClass2, samplesClass2, learnedMean2, learnedCovariance2);

	//keep the components with larger mean distance than twice the covariance in sorted order:
	componentCounter = 0;
	for(i = 0; i < inputDim; i++)
	{
		key = fabs(learnedMean1[i] - learnedMean2[i]) - 2*std::max(sqrt(learnedCovariance1[i*inputDim+i]),
			sqrt(learnedCovariance2[i*inputDim+i]));
		if(key > 0)
		{
			helpList[componentCounter].importance = key;
			helpList[componentCounter].index = i;
			componentCounter++;
		}
	}
	qsort(helpList, inputDim, sizeof(SComponent), SComponent::componentCompare);
std::cout<<componentCounter<<std::endl;
	//copy the indices into the output array:
	for(i = 0; i < componentCounter; i++)
	{
//		if(helpList[i].importance >= 0.001)
//		if(i < 12)
			returnList.push_back(helpList[i].index);
	};
	componentCounter = returnList.size();
std::cout<<componentCounter<<std::endl;
	//populate the outputSample matrix:
	for(i = 0; i < numberSamples; i++)
	{
		std::list<int>::iterator it = returnList.begin();
		for(j = 0; j < componentCounter; j++)
		{
			outputSamples[i*componentCounter+j] = inputSamples[i*inputDim + *it];
			it++;
		}
	}

	delete [] samplesClass1;
	delete [] samplesClass2;
	delete [] learnedMean1;
	delete [] learnedMean2;
	delete [] learnedCovariance1;
	delete [] learnedCovariance2;
	delete [] helpList;

	return returnList;
}

void UnsupervisedLearning::reduceDimensionality(int numTestSamples, int inputDim, double * samples, double * outputSamples, 
		std::list<int> relevantComp)
{
	int i, j, count;
	double * unalignedPoint = new double[inputDim];
	for(i = 0; i < numTestSamples; i++)
	{
		for(j = 0; j < inputDim; j++)
			unalignedPoint[j] = samples[i*inputDim+j];
		//finally store the relevant dimensions in outputSamples:
		std::list<int>::iterator it;
		count = 0;
		for(it = relevantComp.begin(); it != relevantComp.end(); it++, count++)
			outputSamples[i*relevantComp.size()+count] = unalignedPoint[*it];
	}
	delete [] unalignedPoint;
}
