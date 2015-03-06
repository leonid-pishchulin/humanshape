/**
    This file is part of the implementation of the 3D human shape model as described in the paper:

    Leonid Pishchulin, Stefanie Wuhrer, Thomas Helten, Christian Theobalt and Bernt Schiele
    Building Statistical Shape Spaces for 3D Human Modeling
    ArXiv, March 2015

    Please cite the paper if you are using this code in your work.
    
    Contributors: Arjun Jain, Juergen Gall, Leonid Pishchulin.
    This code also uses open source functionality for matrix operations by Thomas Brox.

    The code may be used free of charge for non-commercial and
    educational purposes, the only requirement is that this text is
    preserved within the derivative work. For any other purpose you
    must contact the authors for permission. This code may not be
    redistributed without permission from the authors.
*/

#ifndef __PARAM_MAP__
#define __PARAM_MAP__

#include <map>
//#include <xstring>
#include "CMatrix.h"
#include <nr3.h>
#include <svd.h>
#include <vector>

class paramMap
{
public:
	enum semanticParams {Height,Weight,BreastSize,WaistSize,HipSize,LegLength};

private:
	std::map<semanticParams, bool> stateSemanticParams;	int mActiveParams;
	std::map<semanticParams, std::string> nameSemanticParams;
	std::map<semanticParams, std::pair<int, int> > mSemanticParamRange;
	CMatrix<double> P_mlab, semdata; //note, F is a subset of semdata with addition '1' col
	CMatrix <double> mMap, mMap_T;

private:
	void loadSemdata();
	void loadP_mlab();
	std::string mInputDir;
public:
	paramMap(std::string inputDir);
	paramMap(){paramMap("");};
	void computeMap(); //using param map
	std::string& getNameforParam( unsigned int i0 );
	void setSemanticParams(const std::vector<bool> &states);
	void mapSemanticToEigenSpace( const std::map<paramMap::semanticParams, float> &mSemanticParamValues, CVector<double> &shapeParams );
	int getRange(paramMap::semanticParams param, bool lowEnd);
};

#endif //__PARAM_MAP__
