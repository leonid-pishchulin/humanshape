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

#include "paramMap.h"
#include "onlyDefines.h"

#include <utility>
#include <iosfwd>
#include "CMatrix.h"
// #include <xstring>
#include <string.h>

#include "o_Vector.h"
#include "o_TMatrix.h"
#include "o_TMatrixDoubleHelper.h"

paramMap::paramMap(std::string inputDir)
{
  mInputDir = inputDir;

  mActiveParams = NO_SEMANTIC_PARAMS;

  stateSemanticParams[Height]= true;
  stateSemanticParams[Weight]= true;
  stateSemanticParams[BreastSize]= true;
  stateSemanticParams[WaistSize]= true;
  stateSemanticParams[HipSize]= true;
  stateSemanticParams[LegLength]= true;

  nameSemanticParams[Height]= "Height";
  nameSemanticParams[Weight]= "Weight";
  nameSemanticParams[BreastSize]= "Breast Size";
  nameSemanticParams[WaistSize]= "Waist Size";
  nameSemanticParams[HipSize]= "Hip Size";
  nameSemanticParams[LegLength]= "Leg Length";
  
  mSemanticParamRange[Height] = std::make_pair(-70,70);
  mSemanticParamRange[Weight] = std::make_pair(-100,50);
  mSemanticParamRange[BreastSize] = std::make_pair(-900,500);
  mSemanticParamRange[WaistSize] = std::make_pair(-500,500);
  mSemanticParamRange[HipSize] = std::make_pair(-500,500);
  mSemanticParamRange[LegLength] = std::make_pair(-600,600);
  
  loadSemdata();
  loadP_mlab();
  computeMap();	
}

std::string& paramMap::getNameforParam( unsigned int i0 )
{
	return nameSemanticParams[(semanticParams)i0];
}

void paramMap::loadSemdata()
{
  std::string path = mInputDir + "data/semdata_small.dat";		
    
	std::ifstream loadFile;
	loadFile.open(path.c_str());
	semdata.setSize(228,NO_SEMANTIC_PARAMS);

	unsigned int row = 0, col = 0;
	unsigned totalLines = 1;
	int noPts = 0;
	while(loadFile)
	  {
		float dumm;
		loadFile>>dumm;
		semdata(row,col) = dumm;
		col++;
		if(col == NO_SEMANTIC_PARAMS)
		{
			col = 0;
			row++;
			if(row == 228) break;
		}
	}
	loadFile.close();
}

void paramMap::loadP_mlab()
{	
  std::string path = mInputDir + "data/P_mlab.dat";		
  std::ifstream loadFile;
	loadFile.open(path.c_str());
	P_mlab.setSize(20,228);
	
	unsigned int row = 0, col = 0;
	unsigned totalLines = 1;
	int noPts = 0;
	while(loadFile)
	{
		float dumm;
		loadFile>>dumm;
		P_mlab(row,col) = dumm;
		col++;
		if(col == 228)
		{
			col = 0;
			row++;
			if(row == 20) break;
		}
	}
	loadFile.close();
}



void paramMap::computeMap()
{
  o_TMatrix<double> Finv (mActiveParams + 1, semdata.xSize());
  
  //set F based on  stateSemanticParams
	int row = 0;
	for(unsigned int i0 = 0; i0 < NO_SEMANTIC_PARAMS; i0++)
	{
		if(stateSemanticParams[(semanticParams)i0])
		{
			for(unsigned int i1 = 0; i1 < semdata.xSize(); i1++)			
			{
				float t = semdata(i1,i0);
				Finv[row][i1] = semdata(i1,i0);
			}
			row++;
		}
	}
	for(unsigned int i1 = 0; i1 < semdata.xSize(); i1++)			
	  Finv[row][i1] = 1.0;	
	
	pseudo_inverse(Finv,0,-1.0); //this changes the orientation?!
	mMap.setSize(NO_OF_EIGENVECTORS, mActiveParams + 1);
	mMap_T.setSize(mActiveParams + 1, NO_OF_EIGENVECTORS);
	CMatrix<double> Finv_cmat(semdata.xSize(), mActiveParams + 1);
	//copy TMatrix to CMatrix
	for(unsigned int i0 = 0; i0 < semdata.xSize(); i0++)
		for(unsigned int i1 = 0; i1 < mActiveParams + 1; i1++)
			Finv_cmat(i0,i1) = Finv[i0][i1];

	mMap = Finv_cmat * P_mlab ;	

	for(unsigned int i0 = 0; i0 < mMap.xSize(); i0++)
		for(unsigned int i1 = 0; i1 < mMap.ySize(); i1++)
			mMap_T(i1, i0) = mMap(i0, i1);

}

void paramMap::setSemanticParams(const std::vector<bool> &states)
{
	int active_params = 0;
	for (unsigned int i0 = 0; i0 < states.size(); i0++)
	{
		stateSemanticParams[(semanticParams)i0] = states[i0];
		if(states[i0] == true)
			active_params++;
	}
	mActiveParams =  active_params;
}

void paramMap::mapSemanticToEigenSpace( const std::map<semanticParams, float> &iSemanticParamValues, CVector<double> &shapeParams )
{
	CVector<double> f1; f1.setSize(mActiveParams + 1); f1 = 0; f1(mActiveParams) = 1.0;
	CVector<double> f0; f0.setSize(mActiveParams + 1); f0(mActiveParams) = 1.0;
	for(unsigned int i0 = 0, j0 = 0; i0 < NO_SEMANTIC_PARAMS; i0++)
		if(stateSemanticParams[(semanticParams)i0])
			f0[j0++] = (iSemanticParamValues.find((semanticParams)i0))->second;

	CVector<double> p0, p1; p0.setSize(NO_OF_EIGENVECTORS); p1.setSize(NO_OF_EIGENVECTORS);
	p0 = mMap_T * f0;
	p1 = mMap_T * f1; 
	
	shapeParams = p1 - p0;
}

int paramMap::getRange( semanticParams param, bool lowEnd )
{
	if(lowEnd == true)
		return mSemanticParamRange[param].first;
	else
		return mSemanticParamRange[param].second;
}
