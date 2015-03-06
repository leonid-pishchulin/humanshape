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

#ifndef MLE
#define MLE

#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <vector>

class Mle
{
public:
	Mle();
	~ Mle ();

	void learnExponentialDistribution(int dimension, int sampleSize, double *& samples, double *& learnedLambda);
	void learnGaussDistribution(int dimension, int sampleSize, double *& samples, double *& learnedMean, 
		double *& learnedCovariance);
};

#endif //MLE
