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

#include "ErrorEvaluation.h"
#include <string.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include "time.h"
#include <algorithm> 
using namespace std;
 

ErrorEvaluation::ErrorEvaluation()
{
    dataMatrix = NULL;
}   

ErrorEvaluation::ErrorEvaluation(double *& dataMatrix, int dimension, int numMeshes)
{
    this->numMeshes = numMeshes;
    this->dimension = dimension;
    this->dataMatrix = dataMatrix;
}

ErrorEvaluation::ErrorEvaluation(double *& dataMatrix,double *& newSample, int dimension, int numMeshes)
{
    this->numMeshes  = numMeshes;
    this->dimension  = dimension;
    this->dataMatrix = dataMatrix;
    this->newSample  = newSample;
    
}

ErrorEvaluation::ErrorEvaluation(double *& dataMatrix,double *& gtMatrix, double *& newSample, int dimension, int numMeshes)
{
    this->numMeshes  = numMeshes;
    this->dimension  = dimension;
    this->dataMatrix = dataMatrix;
    this->gtMatrix = gtMatrix;
    this->newSample  = newSample;
    
}

ErrorEvaluation::ErrorEvaluation(double *& dataMatrix,int numMeshes)
{
   this->dataMatrix = dataMatrix;
   this->numMeshes=numMeshes;
   
}

ErrorEvaluation::~ErrorEvaluation()
{
    dataMatrix = NULL;
}

/*Overload Compactness*/
void ErrorEvaluation::compactness(int startNumPrincipalComponents, int endNumPrincipalComponents, std::vector<double> & result, double * ptrEvalues, double * ptrEvectors, double *  ptrMean, double *newSamplePCA, int pcaDim)
{
    result.clear();
    printf("Compactness()\n");
    if((dataMatrix == NULL) || (endNumPrincipalComponents > dimension) || (startNumPrincipalComponents < 0)) return;
    
    if(startNumPrincipalComponents > endNumPrincipalComponents)
    {
        int swap = endNumPrincipalComponents;
        endNumPrincipalComponents = startNumPrincipalComponents;
        startNumPrincipalComponents = swap;
    }
    
    int i, numPrincipalComponents;
    patternRecognitionPCA pca;    
    
    pca.learnPCA(dimension, numMeshes-1, numMeshes, dataMatrix);
    
    double * eigenvalues = pca.getEigenvalues();
    double * reducedMatrix =pca.getReductionMatrix();
    
    double * reducedSample = new double[endNumPrincipalComponents];
    double * reconstructedSample = new double[dimension];
    double * sampleToProject = &(dataMatrix[0]);
    
    double * mean=pca.getMean();
    double * normalizedEigenvalues = new double[numMeshes-1];
    double sum = 0;
    
    for(i = 0; i < numMeshes-1; i++) 
    {
        sum += eigenvalues[i];
        ptrEvalues[i] =eigenvalues[i];
    }
    
    for(i = 0; i < numMeshes-1; i++) normalizedEigenvalues[i]=eigenvalues[i]/sum; //
    
    for(numPrincipalComponents = startNumPrincipalComponents; numPrincipalComponents <= endNumPrincipalComponents; numPrincipalComponents++)
    {
      double returnVal = 0;
      for(i = 0; i < numPrincipalComponents; i++) returnVal += normalizedEigenvalues[i];
      returnVal = returnVal * (sqrt(2.0/(double)numMeshes));
      result.push_back(returnVal);
    }
    
    for(int j = 0; j < dimension; ++j)
      ptrMean[j]=mean[j];
    
    float uppBound=3*sqrt(eigenvalues[0]);
    float lowBound=-3*sqrt(eigenvalues[0]);
    float range = uppBound - lowBound;
    double randomSample;
    
    for(int m=0; m<pcaDim; m++)
      reducedSample[m]=(randomSample*range) +lowBound;
    
    // Saving eigenvectors after PCA
    for(i = 0; i < dimension; i++)
      for(int j = 0; j < (numMeshes-1); j++)
	ptrEvectors[i*(numMeshes-1) + j] = reducedMatrix[i*(numMeshes-1)+ j];
    
}

void ErrorEvaluation::generalization(int startNumPrincipalComponents, int endNumPrincipalComponents, std::vector<double> & distancesOut, int idSample, int numTrainingMeshes, bool bUseGT)
{
  if((dataMatrix == NULL) || (endNumPrincipalComponents > dimension) || (startNumPrincipalComponents < 0)) return;
  
  if(startNumPrincipalComponents > endNumPrincipalComponents)
    {
      int swap = endNumPrincipalComponents;
      endNumPrincipalComponents = startNumPrincipalComponents;
      startNumPrincipalComponents = swap;
    }
  
  int  j, k;
  int l,m,n;
  int nPoints = dimension/3;

  if (numTrainingMeshes > numMeshes - 1)
    numTrainingMeshes = numMeshes - 1;
  
  //Do leave-1-out tests and save error result in a vector:
  double  distances; //new double[numMeshes*dimension+ dimension];
  double * subsample = new double[(numMeshes-1) * dimension];
  double * reducedSample = new double[endNumPrincipalComponents];
  double * reconstructedSample = new double[dimension];
  double auxsum=0.0;
  double sum=0.0;
  
  //Counting Time
  clock_t beginT, endT;
  
  Mle mle1;
  double * learnedMean1 = new double[1];
  double * learnedCovariance1 = new double[1];
    
  if((subsample == NULL) || (reducedSample == NULL) || (reconstructedSample == NULL))
    {
      printf("There is not enough memory to complete this computation\n");
      return;
    }
  
  /*This loop was removed from the original code to allow for parallelization*/
  //for(i = 0; i < numMeshes; i++)
  // {
  
  // live-one-out 
  for(j = 0; j < numMeshes; j++)
    {
      if(j < idSample)
	{
	  for(k = 0; k < dimension; k++)
	    subsample[j*dimension+k] = dataMatrix[j*dimension+k];
	}
      else if(j > idSample)
	{
	  for(k = 0; k < dimension; k++)
	    subsample[(j-1)*dimension+k] = dataMatrix[j*dimension+k];
	}
    }
  
  if (numTrainingMeshes < numMeshes-1){
    double * subsampleSubsample = new double[numTrainingMeshes * dimension];
    std::vector<int> sampleidxs(numMeshes-1);
    for(j = 0; j < numMeshes-1; j++)
      sampleidxs[j] = j;
    
    srand(42u);
    std::random_shuffle(sampleidxs.begin(), sampleidxs.end());
    for(j = 0; j < numTrainingMeshes; j++){
      int jj = sampleidxs[j];
      for(k = 0; k < dimension; k++)
	subsampleSubsample[j*dimension+k] = subsample[jj*dimension+k];
    }
    subsample = subsampleSubsample;
    numMeshes = numTrainingMeshes+1;
  }
  
  patternRecognitionPCA pca;
  pca.learnPCA(dimension, endNumPrincipalComponents, numMeshes-1, subsample);
  double * sampleToProject = &(dataMatrix[dimension*idSample]);
  
  for(j = startNumPrincipalComponents; j <= endNumPrincipalComponents; j++)
    {
      
      cout << "pcaComp " << j << "/" << endNumPrincipalComponents << endl;

      pca.projectToPCA(sampleToProject, reducedSample, j);
      pca.reconstructFromPCA(reducedSample, reconstructedSample, j);//not clear
      
      // compute the distance between the original and reconstructed sample
      distances = 0;
      if (~bUseGT)
	for(k = 0; k < nPoints; k++)
	  {
	    auxsum = pow(dataMatrix[idSample*dimension             + k] - reconstructedSample[k],             2.0) +
	      pow(dataMatrix[idSample*dimension +   nPoints + k] - reconstructedSample[k +   nPoints], 2.0) +
	      pow(dataMatrix[idSample*dimension + 2*nPoints + k] - reconstructedSample[k + 2*nPoints], 2.0);
	    
	    auxsum = sqrt(auxsum);
	    distances += auxsum;
	  }
      else
	for(k = 0; k < nPoints; k++)
	  {
	    auxsum = pow(gtMatrix[idSample*dimension             + k] - reconstructedSample[k],             2.0) +
	      pow(gtMatrix[idSample*dimension +   nPoints + k] - reconstructedSample[k +   nPoints], 2.0) +
	      pow(gtMatrix[idSample*dimension + 2*nPoints + k] - reconstructedSample[k + 2*nPoints], 2.0);
	    
	    auxsum = sqrt(auxsum);
	    distances += auxsum;
	  }
      distances = distances/(double)nPoints;
      distancesOut.push_back(distances);
    }
  
  delete [] subsample;
  delete [] reducedSample;
  delete [] reconstructedSample;
  
 }

/*Compute mean/std for specificity*/
void ErrorEvaluation::learnDistribution(int startNumPrincipalComponents, int endNumPrincipalComponents, std::vector<double> &mean, std::vector<double> &covariance){
    
  mean.clear();
  covariance.clear();
  Mle mle;
  int j;
  double * learnedMean = new double[1];
  double * learnedCovariance = new double[1];
  
  for(j = startNumPrincipalComponents; j <= endNumPrincipalComponents; j++)
    {
      double * distancesForMle = &(dataMatrix[(j-startNumPrincipalComponents)*numMeshes]);
      mle.learnGaussDistribution(1, numMeshes, dataMatrix, learnedMean, learnedCovariance);    
      
      // Note to Leonid: here, you need to report once the mean (learnedMean[0]), and once the covariance (learnedCovariance[0]).
      double c=learnedCovariance[0];
      double me=learnedMean[0];
      covariance.push_back(c);
      mean.push_back(me);
    }    
}

void ErrorEvaluation::specificity(int startNumPrincipalComponents, int endNumPrincipalComponents, std::vector<double> & covariance, std::vector<double> & meanVector, int sampleSize, int numPrincipalComponents,  int numTrainingMeshes, bool bUseGT)
{
  
  printf("ErrorEvaluation::specificity()\n");
  covariance.clear();
  meanVector.clear();
  
  if((dataMatrix == NULL) || (endNumPrincipalComponents > dimension) || (startNumPrincipalComponents < 0)) return;
  
  if(startNumPrincipalComponents > endNumPrincipalComponents)
    {
      int swap = endNumPrincipalComponents;
      endNumPrincipalComponents = startNumPrincipalComponents;
      startNumPrincipalComponents = swap;
    }
  
  int i, j, k;//, numPrincipalComponents;
  double distance, returnVal = 0;
  int l,m,n;
  int nPoints=dimension/3;
  double auxsum=0.0;
  
  double * currentRandomVec;
  double * reducedData = new double[numMeshes * endNumPrincipalComponents];
  double * learnedMean = new double[endNumPrincipalComponents];
  double * learnedCovariance = new double[endNumPrincipalComponents*endNumPrincipalComponents];
  double * randomVec = new double[endNumPrincipalComponents * sampleSize];
  double * randomShape = new double[dimension];
  double * distances = new double[sampleSize];
  double * dataMatrixTrain = new double[numTrainingMeshes * dimension];
  
  if (numTrainingMeshes > numMeshes)
    numTrainingMeshes = numMeshes;

  if((reducedData == NULL) || (learnedMean == NULL) || (learnedCovariance == NULL) || (randomVec == NULL)
     || (randomShape == NULL) || (distances == NULL))
    {
      printf("There is not enough memory to complete this computation\n");
      return;
    }
  
  if (numTrainingMeshes < numMeshes){

    std::vector<int> sampleidxs(numMeshes);
    for(j = 0; j < numMeshes; j++)
      sampleidxs[j] = j;
    
    srand(42u);
    std::random_shuffle(sampleidxs.begin(), sampleidxs.end());
    for(j = 0; j < numTrainingMeshes; j++){
      int jj = sampleidxs[j];
      //std::cout << jj << std::endl;
      for(k = 0; k < dimension; k++)
	dataMatrixTrain[j*dimension+k] = dataMatrix[jj*dimension+k];
    }
  }
  else
    dataMatrixTrain = dataMatrix;
  
  //Sample randomly in PCA space and compare the reconstructed shapes to the given samples
  patternRecognitionPCA pca;
  pca.learnPCA(dimension, endNumPrincipalComponents, numTrainingMeshes, dataMatrixTrain, reducedData);
    
  Mle mle1;
  mle1.learnGaussDistribution(endNumPrincipalComponents, numTrainingMeshes, reducedData, learnedMean, learnedCovariance);
  
  GaussVector gauss(endNumPrincipalComponents, learnedMean, learnedCovariance);
  
  //Generate and store all random vectors. This way, we only generate sampleSize random vectors.
  srand(42u);
  for(i = 0; i < sampleSize; i++)
    {
      currentRandomVec = &(randomVec[i*endNumPrincipalComponents]);
      gauss.generateRandomVector(currentRandomVec);
    }

  /*This loop was removed from the original code to allow parallelization*/  
  // for(numPrincipalComponents = startNumPrincipalComponents; numPrincipalComponents <= endNumPrincipalComponents; numPrincipalComponents++)
   // {
       
  //sampleSize number of samples to calculate PCA
  for(i = 0; i < sampleSize; i++)
    {
      cout << "sample " << i << "/" << sampleSize << endl;
      currentRandomVec = &(randomVec[i*endNumPrincipalComponents]);
      pca.reconstructFromPCA(currentRandomVec, randomShape, numPrincipalComponents);
      
      for(j = 0; j < numMeshes; j++)
	{
	  distance = 0;
	  
	  if (~bUseGT)
	    for(k = 0; k < nPoints; k++)
	      {
		auxsum = pow(dataMatrix[j*dimension             + k] - randomShape[k],             2.0)+
		  pow(dataMatrix[j*dimension +   nPoints + k] - randomShape[k +   nPoints], 2.0)+
		  pow(dataMatrix[j*dimension + 2*nPoints + k] - randomShape[k + 2*nPoints], 2.0);
		auxsum=sqrt(auxsum);
		distance += auxsum;
	      }
	  else
	    for(k = 0; k < nPoints; k++)
	      {
		auxsum = pow(gtMatrix[j*dimension             + k] - randomShape[k],             2.0)+
		  pow(gtMatrix[j*dimension +   nPoints + k] - randomShape[k +   nPoints], 2.0)+
		  pow(gtMatrix[j*dimension + 2*nPoints + k] - randomShape[k + 2*nPoints], 2.0);
		auxsum=sqrt(auxsum);
		distance += auxsum;
	      }
	  distance = distance/(double)nPoints;

	  if((j == 0) || (distance < distances[i])) distances[i] = distance;
	}
    }
  
  //Learn the distribution of the distances:
  Mle mle2;
  mle2.learnGaussDistribution(1, sampleSize, distances, learnedMean, learnedCovariance);
  covariance.push_back(learnedCovariance[0]);
  meanVector.push_back(learnedMean[0]); 
  
  delete [] reducedData;
  delete [] learnedMean;
  delete [] learnedCovariance;
  delete [] randomVec;
  delete [] randomShape;
  delete [] distances;
  
  if (numTrainingMeshes < numMeshes)
    delete [] dataMatrixTrain;

}


/*This function returns the vector as an array for Matlab*/

mxArray * getMexArray(const std::vector<double>& v){
    
   mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
    std::copy(v.begin(), v.end(), mxGetPr(mx));
    return mx;
}




/* Mex function called from Matlab*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
  if (nrhs < 6) {
        mexErrMsgTxt("Error evaluation requires seven input arguments.");
    }
  
  /*input data*/
  double * dataMatrix;
  double * gtMatrix;
  double * dataNewSample;
  int numMeshes;
  int dimension;
  int startNumPrincipalComponents;
  int endNumPrincipalComponents;
  int sampleSize;
  int evalMode;
  int pcaDim;
  int id;
  int numTrainingMeshes;
  
  /*Setting input data*/
  evalMode=(int)mxGetScalar(prhs[0]);
  dataMatrix=mxGetPr(prhs[1]); //matrix with all 3d points
  dataNewSample=mxGetPr(prhs[2]);
  dimension=(int)mxGetScalar(prhs[3]); // # points from a sample in 1D array
  numMeshes=(int)mxGetScalar(prhs[4]); //  # of sampels (scans)
  
  startNumPrincipalComponents=(int)mxGetScalar(prhs[5]);
  endNumPrincipalComponents=(int)mxGetScalar(prhs[6]);
  sampleSize=(int)mxGetScalar(prhs[7]);
  pcaDim=(int)mxGetScalar(prhs[8]);
  
  /* Switch case for compactness, generalization & specificity   */  
  switch(evalMode){
    
  case 0:{ /* Compactness*/
    /*Get input Data*/
    mexPrintf("Compactness()\n");
    mexEvalString("pause(0.001);");
    
    /*output data*/
    std::vector<double> cumulativeStd;
    cumulativeStd.reserve(numMeshes);
    
    /*For saving  eigenvalues and eigenvectors*/
    mwSize m;
    m=(mwSize)(numMeshes-1);//(numMeshes);
    mxArray * evalues = NULL;
    evalues = mxCreateDoubleMatrix(1,m,mxREAL);
    double  * ptrEvalues=mxGetPr(evalues);
    
    mwSize n;
    int x=(dimension)*(numMeshes-1);
    n=(mwSize)(x);
    mxArray * evectors = NULL;
    evectors = mxCreateDoubleMatrix(1,n,mxREAL);
    double  * ptrEvectors=mxGetPr(evectors);
    
    mwSize l;
    l=(mwSize)(dimension);
    mxArray * mean = NULL;
    mean = mxCreateDoubleMatrix(1,l,mxREAL);
    double  * ptrMean=mxGetPr(mean);
    
    mwSize ns;
    int xx=(dimension)*(numMeshes-1);
    ns=(mwSize)(xx);
    mxArray * newSample = NULL;
    newSample = mxCreateDoubleMatrix(1,ns,mxREAL);
    double  * ptrNewSample=mxGetPr(newSample);
    
    ErrorEvaluation *evaluation = new ErrorEvaluation(dataMatrix,dataNewSample, dimension,numMeshes);
    
    evaluation->compactness(startNumPrincipalComponents, endNumPrincipalComponents,cumulativeStd, ptrEvalues, ptrEvectors, ptrMean,ptrNewSample, pcaDim);
    
    plhs[0] = getMexArray(cumulativeStd);
    plhs[1]=evalues; // eigenvalues
    plhs[2]=evectors; //reduction matrix
    plhs[3]=mean; //reduction matrix
    
    std::vector<double>().swap(cumulativeStd);
  }
    break;
    
  case 1:{ /*Generalization*/
    
    /*Get input Data*/
    mexPrintf("Generalization()\n");
    mexEvalString("pause(0.001);");
    id=(int)mxGetScalar(prhs[9]);
    numTrainingMeshes=(int)mxGetScalar(prhs[10]);
    
    /*output Data*/
    std::vector<double> distancesGeneralization;
    distancesGeneralization.reserve(numMeshes);
    
    ErrorEvaluation *evaluation = new ErrorEvaluation(dataMatrix,dataNewSample, dimension,numMeshes);
    evaluation->generalization(startNumPrincipalComponents,endNumPrincipalComponents,distancesGeneralization,id,numTrainingMeshes);
    plhs[0] = getMexArray(distancesGeneralization);
    std::vector<double>().swap(distancesGeneralization);
  }
    break;
    
  case 2:{/*Specificity*/
    
    /*Get input Data*/
    mexPrintf("Specificity()\n");
    mexEvalString("pause(0.001);");
    id=(int)mxGetScalar(prhs[9]);
    numTrainingMeshes=(int)mxGetScalar(prhs[10]);
    
    /*output Data*/
    std::vector<double> Mean;
    Mean.reserve(numMeshes);
    std::vector<double> Covariance;
    Covariance.reserve(numMeshes);
    
    /*For saving  eigenvalues and eigenvectors*/
    mwSize m;
    m=(mwSize)(dimension);//(numMeshes);
    mxArray * evalues = NULL;
    evalues = mxCreateDoubleMatrix(1,m,mxREAL);
    double  * ptrEvalues=mxGetPr(evalues);
    
    ErrorEvaluation *evaluation = new ErrorEvaluation(dataMatrix,dataNewSample, dimension,numMeshes);
    evaluation->specificity(startNumPrincipalComponents,endNumPrincipalComponents,Mean,Covariance,sampleSize,id,numTrainingMeshes);
    plhs[0] = getMexArray(Mean);
    plhs[1] = getMexArray(Covariance);
    std::vector<double>().swap(Mean);
    std::vector<double>().swap(Covariance);
  }
    break;
    
    //Compute mean/std for generalization
  case 3:
    {
      mexPrintf("Compute mean/std\n");
      mexEvalString("pause(0.001);");
      
      /*ouput data*/
      std::vector<double> Covariance;
      Covariance.reserve(1);
      std::vector<double> Mean;
      Mean.reserve(1);
      
      ErrorEvaluation *evaluation = new ErrorEvaluation(dataMatrix,  numMeshes);
      evaluation->learnDistribution(startNumPrincipalComponents,  endNumPrincipalComponents, Mean,Covariance);
      plhs[0] = getMexArray(Mean);
      plhs[1] = getMexArray(Covariance);
      std::vector<double>().swap(Covariance);
      std::vector<double>().swap(Mean);
    }
    break;
    
  case 4:{ /*GeneralizationGT - compute distance to the scan */
    
    /*Get input Data*/
    mexPrintf("GeneralizationGT()\n");
    mexEvalString("pause(0.001);");
    id=(int)mxGetScalar(prhs[9]);
    numTrainingMeshes=(int)mxGetScalar(prhs[10]);
    gtMatrix=mxGetPr(prhs[11]); //matrix with all 3d scan points

    /*output Data*/
    std::vector<double> distancesGeneralization;
    distancesGeneralization.reserve(numMeshes);
    
    ErrorEvaluation *evaluation = new ErrorEvaluation(dataMatrix,gtMatrix,dataNewSample, dimension,numMeshes);
    evaluation->generalization(startNumPrincipalComponents,endNumPrincipalComponents,distancesGeneralization,id,numTrainingMeshes,true);
    plhs[0] = getMexArray(distancesGeneralization);
    std::vector<double>().swap(distancesGeneralization);
  }
    break;

  case 5:{/*SpecificityGT*/
    
    /*Get input Data*/
    mexPrintf("SpecificityGT()\n");
    mexEvalString("pause(0.001);");
    id=(int)mxGetScalar(prhs[9]);
    numTrainingMeshes=(int)mxGetScalar(prhs[10]);
    gtMatrix=mxGetPr(prhs[11]); //matrix with all 3d scan points
    
    /*output Data*/
    std::vector<double> Mean;
    Mean.reserve(numMeshes);
    std::vector<double> Covariance;
    Covariance.reserve(numMeshes);
    
    /*For saving  eigenvalues and eigenvectors*/
    mwSize m;
    m=(mwSize)(dimension);//(numMeshes);
    mxArray * evalues = NULL;
    evalues = mxCreateDoubleMatrix(1,m,mxREAL);
    double  * ptrEvalues=mxGetPr(evalues);
    
    ErrorEvaluation *evaluation = new ErrorEvaluation(dataMatrix,gtMatrix,dataNewSample, dimension,numMeshes);
    evaluation->specificity(startNumPrincipalComponents,endNumPrincipalComponents,Mean,Covariance,sampleSize,id,numTrainingMeshes,true);
    plhs[0] = getMexArray(Mean);
    plhs[1] = getMexArray(Covariance);
    std::vector<double>().swap(Mean);
    std::vector<double>().swap(Covariance);
  }
    break;

    
  default:
    break;
    
  }
            
}
