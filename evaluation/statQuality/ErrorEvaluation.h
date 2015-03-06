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

#ifndef CLASS_ERROR_EVAL
#define CLASS_ERROR_EVAL

#include "patternRecognitionPCA.h"
#include "GeneralizedProcrustes.h"


 extern "C" {
#include "mex.h"
}


/*
Class that evaluates the quality of a statistical model (arbitrary representation).
Code based on error evaluations by Styner et al. (but adapted for model evaluation).

Written by Stefanie Wuhrer, Feb. 2010.
*/

class ErrorEvaluation
{
public:
  ErrorEvaluation();
  ErrorEvaluation(double *& dataMatrix, double *& newSample,int dimension, int numMeshes);
  ErrorEvaluation(double *& dataMatrix, double *& gtMatrix, double *& newSample,int dimension, int numMeshes);
  ErrorEvaluation(double *& dataMatrix, int dimension, int numMeshes);
  ErrorEvaluation(double *& dataMatrix, int numMeshes);
    
  ~ErrorEvaluation();
  
  void compactness(int startNumPrincipalComponents, int endNumPrincipalComponents, std::vector<double> & result);
  void compactness(int startNumPrincipalComponents, int endNumPrincipalComponents, std::vector<double> & result,double * eValues, double *eVectors,double *mean, double *newSample, int pcaDim);
  
  void generalization(int startNumPrincipalComponents, int endNumPrincipalComponents, std::vector<double> & result);
  
  void generalization(int startNumPrincipalComponents, int endNumPrincipalComponents, std::vector<double> & distances, int proc, int numTrainingMeshes, bool bUseGT = false);
  
  void learnDistribution(int startNumPrincipalComponents, int endNumPrincipalComponents, std::vector<double> &mean, std::vector<double> &covariance);
  void specificity(int startNumPrincipalComponents, int endNumPrincipalComponents, std::vector<double> & mean,std::vector<double> & covariance, int sampleSize,int proc, int numTrainingMeshes, bool bUseGT = false);
  
  void learnGauss(int dimension, int sampleSize, double * &samples, double *& learnedMean, double *& learnedCovariance);
  
private:
    double *newSample;
    double * dataMatrix;
    double * gtMatrix;
    int numMeshes, dimension;
};

#endif
