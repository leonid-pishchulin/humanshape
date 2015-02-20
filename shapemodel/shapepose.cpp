/**
    This file is part of the implementation of the 3D human shape model as described in the paper:

    Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobald and Bernt Schiele
    Building Statistical Shape Spaces for 3D Human Modeling
    ArXiv, February 2015

    Please cite the paper if you are using this code in your work.
    
    Contributors: Arjun Jain, Juergen Gall, Leonid Pishchulin.
    This code also uses open source functionality for matrix operations by Thomas Brox.

    The code may be used free of charge for non-commercial and
    educational purposes, the only requirement is that this text is
    preserved within the derivative work. For any other purpose you
    must contact the authors for permission. This code may not be
    redistributed without permission from the authors.
*/

#include <vector>
#include <exception>

#include "CMatrix.h"
#include "Show.h"

#include "NMath.h"
#include "NRBM.h"
#include "CTMesh.h"
#include "onlyDefines.h"
#include "paramMap.h"

#ifdef USEGPU
#include "GL/glew.h"
#include "GL/glut.h"
#endif

#include "time.h"
#include "stdlib.h"
#include "mex.h"

using namespace std;

#define SMOOTHMODEL true

int getpose(const string inputDir, const double* motionParamsIn, const double* shapeParamsIn, double *eigenVectorsIn,
	    int numEigenVectors, double* pointsOut, double* jointsOut){

  // ************** Read Input ************** //
  string s = inputDir;
  if( (char)(*--s.end()) != '/' )
    s.append("/");
  NShow::mInputDir = s;
  
#ifndef EXEC_FAST
  cout << "Input dir: " << NShow::mInputDir << endl;
#endif

#ifndef EXEC_FAST
  cout << "LOAD MESH...\n";
#endif 
  std::string aModelFile = NShow::mInputDir + "model.dat";
  // read motion params from the precomputed 3D poses
  const int numMotionParams = 31;
  CVector<double> mParams(numMotionParams);
  for(int i = 0; i < numMotionParams;i++)
    mParams(i) = motionParamsIn[i];
#ifndef EXEC_FAST
  cout << "3D pose parameters:" << endl;
  cout << mParams << endl;
#endif
  
#ifndef EXEC_FAST
  cout << "semantic shape parameters: ";
#endif
  CVector<double> shapeParams(numEigenVectors);  
  
  for(unsigned int i0 = 0; i0 < numEigenVectors; i0++){
    shapeParams[i0] = shapeParamsIn[i0];
  }

  CVector<float> shapeParamsFloat(numEigenVectors);
  for(int i=0;i<numEigenVectors;i++)
    shapeParamsFloat(i) = float(shapeParams(i));
  
  // Read object model
  CMatrix<float> mRBM(4,4);
  NRBM::RVT2RBM(&mParams, mRBM);	
  CMesh initMesh;
  
  initMesh.readModel(aModelFile.c_str(), SMOOTHMODEL);
  initMesh.updateJntPos();
  
  initMesh.centerModel();
  
  // use eingevectors passed as an argument
  initMesh.readShapeSpaceEigens(eigenVectorsIn,numEigenVectors,initMesh.GetPointSize());
  
  // reshape the model
  initMesh.shapeChangesToMesh(shapeParamsFloat);

  // update joints
  initMesh.updateJntPos();
  
  CVector<CMatrix<float> > M(initMesh.joints()+1);
  CVector<float> TW(initMesh.joints()+6);
  float tmp[100];
  
  for(int j=6; j<mParams.size();++j) {
    tmp[j] = (float)mParams(j);
    TW(j) = (float)mParams(j);
  }

  initMesh.angleToMatrix(mRBM, TW, M);

  // rotate joints
  initMesh.rigidMotion(M, TW, true, true);

  // Fill in resulting joints array
  int nJoints = initMesh.joints();
  for (int i = 0; i < nJoints; i++){
    CJoint joint = initMesh.joint(i+1);
    jointsOut[i            ] = i + 1;
    jointsOut[i + nJoints  ] = joint.getDirection()[0];
    jointsOut[i + nJoints*2] = joint.getDirection()[1];
    jointsOut[i + nJoints*3] = joint.getDirection()[2];
    jointsOut[i + nJoints*4] = joint.getPoint()[0];
    jointsOut[i + nJoints*5] = joint.getPoint()[1];
    jointsOut[i + nJoints*6] = joint.getPoint()[2];
    jointsOut[i + nJoints*7] = double(joint.mParent);
  }
  
  // Fill in resulting points array
  int nPoints = initMesh.GetPointSize();
  for (int i = 0; i < nPoints; i++){
    float x,y,z;
    initMesh.GetPoint(i,x,y,z);
    pointsOut[i            ] = x;
    pointsOut[i +   nPoints] = y;
    pointsOut[i + 2*nPoints] = z;
  }
  return 0;
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {
  int m, n;
  
  double *pointsOut; /* 6449x3*/
  double *jointsOut; /* 24x8: jid directions_XYZ positions_XYZ jparent_id*/
  double *motionParamsIn;  /* 20x1*/
  double *shapeParamsIn;  /* 31x1*/
  double *eigenVectorsIn; /* nEigenVec x 6449 x 3*/
  
  if (nrhs < 4){
    mexErrMsgTxt("Not enough input arguments");
    return;
  }
  
  /* Find the dimensions of the data */
  m = 6449; /*mxGetM(prhs[0]);*/
  n = 3; /*mxGetN(prhs[0]);*/
  
  /* Create an mxArray for the output data */
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);

  /* Create a pointer to the output data */
  pointsOut = mxGetPr(plhs[0]);
  
  m = 25; n = 8;
  
  plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
  
  /* Create a pointer to the output data */
  jointsOut = mxGetPr(plhs[1]);

  /* Retrieve the input data */
  motionParamsIn = mxGetPr(prhs[0]);
  shapeParamsIn = mxGetPr(prhs[1]);
  eigenVectorsIn= mxGetPr(prhs[2]);
  int numEigenVectors = mxGetM(prhs[2]);
  
  int modelDirNameLen = (mxGetM(prhs[3]) * mxGetN(prhs[3])) + 1;
  char *modelDirName = (char*) malloc (modelDirNameLen);
  int status = mxGetString(prhs[3], modelDirName, modelDirNameLen); 
  if (status != 0){
    mexErrMsgTxt("modelDirName is wrong\n");
  }

  getpose(modelDirName, motionParamsIn, shapeParamsIn, eigenVectorsIn, numEigenVectors, pointsOut, jointsOut);
} 
