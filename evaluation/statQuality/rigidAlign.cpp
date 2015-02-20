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

#include "GeneralizedProcrustes.h"
#include <string.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include "time.h"
using namespace std;

/*This function returns the vector as an array for Matlab*/
mxArray * getMexArray(const std::vector<double>& v){
    
   mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
    std::copy(v.begin(), v.end(), mxGetPr(mx));
    return mx;
}


/* Mex function called from Matlab*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
  /*input data*/
  double * dataMatrix;
  int numMeshes;
  int dimension;
    
  /*Setting input data*/
  dataMatrix=mxGetPr(prhs[0]); //matrix with all 3d points
  dimension=(int)mxGetScalar(prhs[1]); // # points from a sample in 1D array
  numMeshes=(int)mxGetScalar(prhs[2]); //  # of sampels (scans)

  int numVertices=dimension/3;
  
  /* output data */
  plhs[0] = mxCreateDoubleMatrix(dimension,numMeshes,mxREAL);
  double *data = mxGetPr(plhs[0]);
    
  vector < Shape > shapes(numMeshes);
  mwSize n;
  
  //Stack 3D points
  for(int i = 0; i < numMeshes; i++){
    
    int startCurrMesh = i*dimension;
    shapes[i].vertexCoordinates.resize(numVertices);
    
    for(int j = 0; j < numVertices; j++){     
      shapes[i].vertexCoordinates[j].resize(3);
      shapes[i].vertexCoordinates[j][0] = dataMatrix[startCurrMesh + j                ];
      shapes[i].vertexCoordinates[j][1] = dataMatrix[startCurrMesh + j +   numVertices];
      shapes[i].vertexCoordinates[j][2] = dataMatrix[startCurrMesh + j + 2*numVertices];
      shapes[i].numberVertices = numVertices;
    }
  }
  
  mexPrintf("ComputeGPAlignment()\n");
  mexEvalString("pause(0.001);");
  GeneralizedProcrustes *alignment= new GeneralizedProcrustes(shapes, numMeshes, numVertices);
  alignment->computeGPAlignment(shapes);
  
  for(int i = 0; i < numMeshes; i++){
    int startCurrMesh = i*dimension;
    for(int j = 0; j < numVertices; j++){
      data[startCurrMesh + j                ] = shapes[i].vertexCoordinates[j][0];
      data[startCurrMesh + j +   numVertices] = shapes[i].vertexCoordinates[j][1];
      data[startCurrMesh + j + 2*numVertices] = shapes[i].vertexCoordinates[j][2];
    }
  }
  alignment->~GeneralizedProcrustes();
}
