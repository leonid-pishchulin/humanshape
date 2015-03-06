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

#ifndef PCA_CLASS
#define PCA_CLASS

#include "UnsupervisedLearning.h"

class patternRecognitionPCA
{
public:
	patternRecognitionPCA();
	~patternRecognitionPCA();

	//Method learns the PCA space given a set of input samples:
	void learnPCA(int inputDim, int outputDim, int numberSamples, double * &samples);
	//Like before, but gives back the reduced samples:
	void learnPCA(int inputDim, int outputDim, int numberSamples, double *& samples, double  *& reducedSamples);
   // void learnPCA(int inputDim, int outputDim, int numberSamples, double * samples, double * reducedSamples);


    
	//Method takes a new sample and projects it to the PCA space: returns reducedSample
	//void projectToPCA(double * newSample, double * reducedSample);
	//Method takes PCA coordinates and transforms them back to original space: returns outputSample
	//void reconstructFromPCA(double * newRedSample, double * outputSample);

    void projectToPCA(double *& newSample, double *& reducedSample);
	void reconstructFromPCA(double *& newRedSample, double *& outputSample);

   
    
	//Same methods as above, but only a subset of the PCA space is used:
	void projectToPCA(double *& newSample, double *& reducedSample, int pcaDim);
	void reconstructFromPCA(double *& newRedSample, double *& outputSample, int pcaDim);

	//Methods to set the information:
	void setAllInfo(int inputDim, int outputDim, double *& reductionMatrix, double *& mean);

	//Methods for import and export:

	//Give the filename (including path) where the PCA should be exported to. Method expects path to be valid.
	//Returns true if everything went well and false otherwise.
	bool exportPCA(char * filename);
	//Give the filename (including path) where the PCA should be imported from. Method expects path to be valid.
	//Returns true if everything went well and false otherwise.
	bool importPCA(char * filename);
    
	double * getEigenvalues(){return eigenvalues;}
	int getOutputDim(){return outputDim;}
	
    double * getReductionMatrix(){return reductionMatrix;}
    double * getMean(){return mean;}
    double * getEigenvectors(){ return eigenvectors;}

private:

	void readComment(FILE *& fp);

	int inputDim, outputDim;
	//reductionMatrix is the normalized eigenvalue matrix:
	double * reductionMatrix;
	double * mean;
	double * eigenvalues;
    double * eigenvectors;
};

#endif
