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
#include "patternRecognitionPCA.h"

 extern "C" {
#include <math.h>
}

patternRecognitionPCA::patternRecognitionPCA()
{
	reductionMatrix = NULL;
	mean = NULL;
	eigenvalues = NULL;
}

patternRecognitionPCA::~patternRecognitionPCA()
{
	if(reductionMatrix != NULL)
	{
		delete [] reductionMatrix;
		reductionMatrix = NULL;
	}
	if(mean != NULL)
	{
		delete [] mean;
		mean = NULL;
	}
	if(eigenvalues != NULL)
	{
		delete [] eigenvalues;
		eigenvalues = NULL;
	}
}

void patternRecognitionPCA::learnPCA(int inputDim, int outputDim, int numberSamples, double * &samples)
{
	int i, j;
	
	this->inputDim = inputDim;
	this->outputDim = outputDim;

	UnsupervisedLearning learn(inputDim, 1, numberSamples);
	eigenvalues = new double[outputDim];
	reductionMatrix = new double[inputDim * outputDim];
	mean = new double[inputDim];

	learn.performPCASpaceRestrictedSVD(inputDim, outputDim, samples, reductionMatrix, eigenvalues, mean);
        
 
}

// void patternRecognitionPCA::learnPCA(int inputDim, int outputDim, int numberSamples, double * samples, double * reducedSamples)
// {
// 	int i, j;
//     
//     mexWarnMsgTxt("Start patternRecognitionPCA\n");
// 
// 	this->inputDim = inputDim;
// 	this->outputDim = outputDim;
// 	UnsupervisedLearning learn(inputDim, 1, numberSamples);
// 	eigenvalues = new double[outputDim];
// 	reductionMatrix = new double[inputDim * outputDim];
// 	mean = new double[inputDim];
//     mexWarnMsgTxt("patternRecognitionPCA\n");
// 	learn.performPCA(inputDim, outputDim, samples, reducedSamples, reductionMatrix, eigenvalues, mean);
// 
//     
// }



void patternRecognitionPCA::learnPCA(int inputDim, int outputDim, int numberSamples, double *& samples, double *& reducedSamples)
{
	int i, j;

	this->inputDim = inputDim;
	this->outputDim = outputDim;

	UnsupervisedLearning learn(inputDim, 1, numberSamples);
	eigenvalues = new double[outputDim];
	reductionMatrix = new double[inputDim * outputDim];
	mean = new double[inputDim];

	//learn.performPCA(inputDim, outputDim, samples, reducedSamples, reductionMatrix, eigenvalues, mean);
	learn.performPCASpaceRestrictedSVD(inputDim, outputDim, samples, reductionMatrix, eigenvalues, mean);

	//Compute the reduced samples:
	double * sample = new double[inputDim];
	double * reducedSample = new double[outputDim];

	for(i = 0; i < numberSamples; i++)
	{
		for(j = 0; j < inputDim; j++) 
			sample[j] = samples[i*inputDim+j];

		projectToPCA(sample, reducedSample);

		for(j = 0; j < outputDim; j++)
			reducedSamples[i*outputDim+j] = reducedSample[j];
	}

	delete [] sample;
	delete [] reducedSample;
}






void patternRecognitionPCA::projectToPCA(double *& newSample, double *& reducedSample)
{
	int i;

	char transN = 'N';
	ptrdiff_t n = (ptrdiff_t) inputDim;
	ptrdiff_t m = (ptrdiff_t) outputDim;
	double alpha = 1.0;
	ptrdiff_t incx = 1;
	double beta = 0.0;

	double * modifiedSample = new double[inputDim];

	for(i = 0; i < inputDim; i++) modifiedSample[i] = newSample[i] - mean[i];

	//Use the reduction matrix to add a sample to the PCA space:
	dgemv(&transN, &m, &n, &alpha, reductionMatrix, &m, modifiedSample, &incx, &beta, reducedSample, &incx);
	
	delete [] modifiedSample;
}

void patternRecognitionPCA::projectToPCA(double *& newSample, double *& reducedSample, int pcaDim)
{
	if(pcaDim > outputDim) return;
    
	else if(pcaDim == outputDim) return projectToPCA(newSample, reducedSample);
	else
	{
		int i, j;

		double * modifiedSample = new double[inputDim];

		char transN = 'N';
		ptrdiff_t n = (ptrdiff_t) inputDim;
		ptrdiff_t m = (ptrdiff_t) pcaDim;
		double alpha = 1.0;
		ptrdiff_t incx = 1;
		double beta = 0.0;

		double * smallReductionMatrix = new double[inputDim * pcaDim];
		for(i = 0; i < inputDim; i++)
		{
			for(j = 0; j < pcaDim; j++) smallReductionMatrix[i*pcaDim + j] = reductionMatrix[i*outputDim + j];
		}

		for(i = 0; i < inputDim; i++) modifiedSample[i] = newSample[i] - mean[i];

		//Use the reduction matrix to add a sample to the PCA space:
        //y := alpha*smallReductionMatrix*inc + beta*y.
		dgemv(&transN, &m, &n, &alpha, smallReductionMatrix, &m, modifiedSample, &incx, &beta, reducedSample, 
			&incx);

		delete [] modifiedSample;

		delete [] smallReductionMatrix;
	}
}	

void patternRecognitionPCA::reconstructFromPCA(double *& newRedSample, double *& outputSample)
{
	int i;

	char transY = 'T';
	ptrdiff_t inc = 1;
	ptrdiff_t m = (ptrdiff_t) outputDim;
	ptrdiff_t n = (ptrdiff_t) inputDim;
	double alpha = 1.0;
	double beta = 0.0;

	//Use the transposed reduction matrix here:
    // MONICA-y := alpha*reductionMatrix'*newRedSample + beta*y
	dgemv(&transY, &m, &n, &alpha, reductionMatrix, &m, newRedSample, &inc, &beta, outputSample, &inc);

        
	for(i = 0; i < inputDim; i++) outputSample[i] = mean[i] + outputSample[i];
}

void patternRecognitionPCA::reconstructFromPCA(double *& newRedSample, double *& outputSample, int pcaDim)
{
	if(pcaDim > outputDim) return;
	else if(pcaDim == outputDim) return reconstructFromPCA(newRedSample, outputSample);
	else
	{
		int i, j;

		char transY = 'T';
		ptrdiff_t inc = 1;
		ptrdiff_t m = (ptrdiff_t) pcaDim;
		ptrdiff_t n = (ptrdiff_t) inputDim;
		double alpha = 1.0;
		double beta = 0.0;

		double * smallReductionMatrix = new double[inputDim * pcaDim];
		for(i = 0; i < inputDim; i++)
		{
			for(j = 0; j < pcaDim; j++) smallReductionMatrix[i*pcaDim + j] = reductionMatrix[i*outputDim + j];
		}

        
		//Use the transposed reduction matrix here:
        //size ( 1 + ( inputDim- 1 ) )
		dgemv(&transY, &m, &n, &alpha, smallReductionMatrix, &m, newRedSample, &inc, &beta, outputSample, &inc);

		for(i = 0; i < inputDim; i++) outputSample[i] = mean[i] + outputSample[i];
       // printf(" smallReductionMatrix%f\n", smallReductionMatrix[19347]);
       //  printf(" newRedSample%f\n", newRedSample[pcaDim]);
		delete [] smallReductionMatrix;
	}
    
       
}

void patternRecognitionPCA::setAllInfo(int inputDim, int outputDim, double *& reductionMatrix, double *& mean)
{
    
	int i;
	this->inputDim = inputDim;
	this->outputDim = outputDim;
	this->reductionMatrix = new double[inputDim * outputDim];
	for(i = 0; i < inputDim*outputDim; i++) this->reductionMatrix[i] = reductionMatrix[i];
	this->mean = new double[inputDim];
	for(i = 0; i < inputDim; i++) this->mean[i] = mean[i];
}

bool patternRecognitionPCA::exportPCA(char * filename)
{
	int i;

	FILE * fp = fopen(filename, "w");
	if((fp == NULL) || (mean == NULL) || (reductionMatrix == NULL)) return false;

	fprintf(fp, "%%File created by patternRecognitionPCA\n");
	fprintf(fp, "%%Line 1: input dimension, output dimension\n");
	fprintf(fp, "%d %d\n", inputDim, outputDim);
	fprintf(fp, "%%Line 2: mean\n");
	for(i = 0; i < inputDim; i++) fprintf(fp, "%f ", mean[i]);
	fprintf(fp, "\n");
	fprintf(fp, "%%Line 3: normalized eigenvector matrix in column-major order\n");
	for(i = 0; i < inputDim*outputDim; i++) fprintf(fp, "%f ", reductionMatrix[i]);
	fprintf(fp, "\n");
	fprintf(fp, "%%Line 4: eigenvalues\n");
	for(i = 0; i < outputDim; i++) fprintf(fp, "%f ", eigenvalues[i]);
	fprintf(fp, "\n");

	fclose(fp);

	return true;
}

bool patternRecognitionPCA::importPCA(char * filename)
{
	int i, intRead1, intRead2;
	float floatRead;

	FILE * fp = fopen(filename, "r");
	if(fp == NULL) return false;

	readComment(fp);
	fscanf(fp, "%d %d ", &intRead1, &intRead2);
	inputDim = intRead1;
	outputDim = intRead2;

	if(mean != NULL) delete [] mean;
	mean = new double[inputDim];
	readComment(fp);
	for(i = 0; i < inputDim; i++) 
	{
		fscanf(fp, "%f ", &floatRead);
		mean[i] = floatRead;
	}

	if(reductionMatrix != NULL) delete [] reductionMatrix;
	reductionMatrix = new double[inputDim*outputDim];
	readComment(fp);
	for(i = 0; i < inputDim*outputDim; i++) 
	{
		fscanf(fp, "%f ", &floatRead);
		reductionMatrix[i] = floatRead;
	}

	if(eigenvalues != NULL) delete [] eigenvalues;
	eigenvalues = new double[outputDim];
	readComment(fp);
	for(i = 0; i < outputDim; i++)
	{
		fscanf(fp, "%f ", &floatRead);
		eigenvalues[i] = floatRead;
	}

	fclose(fp);

	return true;
}

void patternRecognitionPCA::readComment(FILE *& fp)
{
	int c;
	bool comment;
	char buffer[255];

	do
	{
		c = getc (fp);
		if (c == '%')
		{
			ungetc (c, fp);
			fgets (buffer, 255, fp);
			comment = true;
		}
		else
		{
			ungetc (c, fp);
			comment = false;
		}
	}
	while(comment);
}
