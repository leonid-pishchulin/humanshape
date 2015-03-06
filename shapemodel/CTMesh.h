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

#ifndef CMeshH
#define CMeshH

#include <fstream>
#include "CMatrix.h"
#include "CTensor.h"
#include <iostream>
#include "Show.h"

#define NO_OF_EIGENVECTORS 20
#define SKEL_SCALE_FACTOR 1.0

#define EXEC_FAST 1
//#include "RenderWidget.h"

class CJoint 
{
public:
	// constructor
	inline CJoint() {mParent = 0;mDirection.setSize(0);mPoint.setSize(0);mMoment.setSize(0);  };
	~CJoint(){};
	inline CJoint(CJoint& aCopyFrom) {*this = aCopyFrom;};
	CJoint(CVector<float>& aDirection, CVector<float>& aPoint, int aParent);
	// Performs a rigid motion M of the joint
	void rigidMotion(CMatrix<float>& M);
	// Constructs the motion matrix from the joint axis and the rotation angle
	void angleToMatrix(float aAngle, CMatrix<float>& M);
	// Access to joint's position and axis
	inline void set(CVector<float>& aDirection, CVector<float>& aPoint) {mDirection = aDirection; mPoint = aPoint; mMoment = aPoint/aDirection;};
	inline void setDirection(CVector<float>& aDirection) {mDirection = aDirection; mMoment = mPoint/aDirection;};
	inline CVector<float>& getDirection() {return mDirection;};
	inline void setPoint(CVector<float>& aPoint) {mPoint = aPoint; mMoment = aPoint/mDirection;};
	inline CVector<float>& getPoint() {return mPoint;};
	inline CVector<float>& getMoment() {return mMoment;};
	// Copy operator
	CJoint& operator=(CJoint& aCopyFrom);
	// Parent joint
	int mParent;
protected:
	// Defines joint's position and axis
	CVector<float> mDirection;
	CVector<float> mPoint;
	CVector<float> mMoment;
};

class CMeshMotion {
public:
	// constructor
	inline CMeshMotion() {};
	// Resets the model to zero motion
	void reset(int aJointNumber);
	// Shows current configuration
	void print();
	// Writes current configuration to file
	void writeToFile();
	// Copy operator
	CMeshMotion& operator=(const CMeshMotion& aCopyFrom);
	// Access to pose parameters
	inline float& operator()(int aIndex) const {return mPoseParameters(aIndex);};

	// Main body motion
	CMatrix<float> mRBM;
	// Vector with all pose parameters (including joint angles)
	CVector<float> mPoseParameters;
};

/* class CAppearanceModel { */
/* public: */
/*   inline CAppearanceModel() {}; */
/*   CVector<int> mLastFrameVisible; */
/*   CTensor<float> mHistogram; */
/* }; */


class CMesh {
public:
	void writeMeshDat( std::string fname );
	void printPoints( std::string fname );
	int shapeChangesToMesh( CVector<float> shapeParams );
	int updateJntPos();
	static void readShapeSpaceEigens( std::string fileName, int numEigenVectors);
	static void readShapeSpaceEigens(const double* eigenVectors, int numEigenVectors, int nPoints);
	// constructor
	inline CMesh();
	inline CMesh(const CMesh& aMesh);
	CMesh(int aPoints, int aPatches);
	// destructor
	~CMesh();
	// Reads the mesh from a file
	// liuyebin collect the readmodel functions together. The first is for initISA, and the second is for trackHybOpt
	bool readModel(const char* aFilename, bool smooth = false);
	bool readModelHybOpt(const char* aFilename, bool smooth = false);
	bool adaptOFF(const char* aFilename, float lambda);
	bool readOFF(const char* aFilename);
	void centerModel();
	// Writes the mesh to a file
	void writeModel(const char* aFilename);
	void writeAdaptModel(const char* aFilename, const CMesh* adM);
	void writeSkel(const char* aFilename);

	//void draw();

	// Performs rigid body motions M of the mesh points in the kinematic chain
	void rigidMotion(CVector<CMatrix<float> >& M,CVector<float>& X, bool smooth = false, bool force = false);
	void smoothMotionDQ(CVector<CMatrix<float> >& M,CVector<float>& X);
	// Reuses InitMesh to set up Smooth Pose: Global transformation
	void makeSmooth(CMesh* initMesh, bool dual = false);

	void angleToMatrix(const CMatrix<float>& aRBM, CVector<float>& aJAngles, CVector<CMatrix<float> >& M);
	void invAngleToMatrix(const CMatrix<float>& aRBM, CVector<float>& aJAngles, CVector<CMatrix<float> >& M);
	void twistToMatrix(CVector<float>& aTwist, CVector<CMatrix<float> >& M);

	// Fast projection of a 3-D point to the image plane
	inline void projectPoint(CMatrix<float>& P, float X, float Y, float Z, int& x, int& y);
	inline void projectPoint(CMatrix<float>& P, float X, float Y, float Z, float& x, float& y);

	// Indicates a new frame, causing copying the current motion into the motion history
	//void newFrame();

	// Copies aMesh
	void operator=(const CMesh& aMesh);

	void setJointDir(int aJointID, CVector<float> dir) {
		mJoint(aJointID).setDirection(dir);
	};

	// Returns the number of joints
	int joints() {return mJointNumber;};
	// Returns a joint
	CJoint& joint(int aJointID) {return mJoint(aJointID);};
	// turns whether a point is influenced by a certain joint
	bool influencedBy(int aJointIDOfPoint, int aJointID) { return mInfluencedBy(aJointIDOfPoint,aJointID);};
	bool isNeighbor(int i,int j) {return mNeighbor(mJointMap(i),mJointMap(j));};
	bool isEndJoint(int aJointID) {return mEndJoint[aJointID];}

	inline void GetPoint(int i,  float& x, float& y, float& z);
	inline void GetPoint(int i,  float& x, float& y, float& z, int& j);
	inline void GetPoint(int i,  float& x, float& y, float& z, float& j);
	inline void GetPatch(int i, int& x, int& y, int& z);
	inline void GetBounds(int J, int i, float& x, float& y, float& z);
	int GetBoundJID(int J) {return (int)mBounds(J,8,0);};
	float GetCenter(int i) {return mCenter[i];};
	int GetBodyPart(int jID) {return mJointMap[jID];};

	int GetMeshSize() {return mNumPatch;};
	int GetPointSize() {return mNumPoints;};
	int GetBoundSize() {return mBounds.xSize();};

	int GetJointID(int i) {return (int)mPoints[i](3);};

	bool IsCovered(int i) {return mCovered[i];};
	bool IsExtremity(int i) {return mExtremity[i];};

	void projectSurface(CMatrix<float>& aImage, CMatrix<float>& P, int xoffset = 0, int yoffset = 0);


	// Projects the mesh to the image plane given the projection matrix P
	void projectToImage(CMatrix<float>& aImage, CMatrix<float>& P, int aLineVal = 255);
	void projectToImage(CMatrix<float>& aImage, CMatrix<float>& P, int aLineVal, int aNodeVal);
	void projectToImage(CTensor<float>& aImage, CMatrix<float>& P, int aLineR = 255, int aLineG = 255, int aLineB = 255);
	void projectToImage(CTensor<float>& aImage, CMatrix<float>& P, int aLineR, int aLineG, int aLineB, int aNodeR, int aNodeG, int aNodeB);
	// Just project a joint index + neighbor ...
	void projectToImageJ(CMatrix<float>& aImage, CMatrix<float>& P, int aLineValue, int aJoint, int xoffset = 0, int yoffset = 0);

	// Projects all mesh points to the image plane
	void projectPointsToImage(CMatrix<float>& aImage, CMatrix<float>& P, int aNodeVal = 128);
	// Projects all mesh points to the image plane and returns the 3-D coordinates of these points
	void projectPointsToImage(CTensor<float>& a3DCoords, CMatrix<float>& P);

	// Initialization of appearance model
	//void initializeAppearanceModel(int aFeatureSize, int aViewCount);
	// Update of appearance model and occlusion detection
	//void updateAppearance(CTensor<float>& aData, CMatrix<float>& aOccluded, CMatrix<float>& P, int aViewNo, int aFrameNo, bool aOcclusionCheckOnly); 

	// Writes accumulated motion to standard output
	//inline void printAccumulatedMotion() {mAccumulatedMotion.print();};
	// Writes accumulated motion to file
	//inline void writeAccumulatedMotion() {mAccumulatedMotion.writeToFile();};
	// Resets accumulated motion to 0
	inline void resetAccumulation() {mAccumulatedMotion.reset(mJointNumber);};
	inline void resetCurrentMotion() {mCurrentMotion.reset(mJointNumber);};
	// Gives access to current motion
	inline CMeshMotion& currentMotion() {return mCurrentMotion;};
	// Gives access to pose history
	//inline CVector<CMeshMotion>& history() {return mHistory;};

protected:

	int mJointNumber;

	std::vector<CVector<float> >  mPoints;
	std::vector<CVector<int> >  mPatch;
	CTensor<float> mBounds;
	CVector<float> mCenter;
	CVector<int> mJointMap;
	CMatrix<bool> mNeighbor;
	CVector<bool> mEndJoint;

	CVector<bool> mCovered;
	CVector<bool> mExtremity;

	int mNumPoints;
	int mNumPatch;
	int mNumSmooth; // how many joints can influence any given point
	int mNoOfBodyParts;
	CVector<bool> mBoundJoints;

	CMatrix<bool> mInfluencedBy;
	// CMatrix<float> mRBM; // Overall Rigid Body Motion;

	//CMatrix<CAppearanceModel>* mAppearanceField;
	// CTensor<bool>* mOccluded;
	CMeshMotion mAccumulatedMotion;
	CMeshMotion mCurrentMotion;
	//CVector<CMeshMotion> mHistory;
	CVector<CJoint>  mJoint;
	static std::vector<CMatrix<double> >  eigenVectors;
	static CMatrix<float> weightMatrix;

	// Sure Fields before the Draping ...
	// CVector<float>* mPointsS;

	// true if aParentJoint is an ancestor of aJoint
	bool isParentOf(int aParentJointID, int aJointID);
	//~aj - small functions, can move elsewhere (think)
	bool minMLab( CMatrix<float> weightMatrix, int i0, int i1, CVector<float> &minEle );
	void componet_wise_mul_with_pnts( CVector<float> minEle, CMatrix<float> &tmp_wtMatrix );
	double sumTheVector( CVector<float> minEle );
	void sumTheMatrix( CMatrix<float> tmpMatrix, CVector<float> & t1 );
	bool copyMatColToVector( CMatrix<float> weightMatrix, int i0, CVector<float> &singleWts );
	void copyJointPos( int to, int from, CMatrix<float> &newJntPos );
	void findNewAxesOrient( CVector<double> &axisLeftArm, CVector<double> &axisRightArm, CMatrix<float> newJntPos );	
	void getWeightsinBuffer( char *buffer, int ptNo );
};

// constructor
inline CMesh::CMesh() {
	//mPoints=0;
	//mPointsS=0;
	//mPatch=0;
	//mAppearanceField = 0;
	//mOccluded = 0;
}

inline CMesh::CMesh(const CMesh& aMesh) {

	//mPoints=0;
	//mPointsS=0;
	//mPatch=0; 
	//mAppearanceField = 0;
	//mOccluded = 0;
	*this = aMesh;
}

// projectPoint
inline void CMesh::projectPoint(CMatrix<float>& P, float X, float Y, float Z, int& x, int& y) {

	float hx = P.data()[0]*X + P.data()[1]*Y + P.data()[2]*Z + P.data()[3];
	float hy = P.data()[4]*X + P.data()[5]*Y + P.data()[6]*Z + P.data()[7];
	float hz = P.data()[8]*X + P.data()[9]*Y + P.data()[10]*Z + P.data()[11];


	float invhz = 1.0/hz;
	x = (int)(hx*invhz+0.5);
	y = (int)(hy*invhz+0.5);
}

inline void CMesh::projectPoint(CMatrix<float>& P, float X, float Y, float Z, float& x, float& y) {

	float hx = P.data()[0]*X + P.data()[1]*Y + P.data()[2]*Z + P.data()[3];
	float hy = P.data()[4]*X + P.data()[5]*Y + P.data()[6]*Z + P.data()[7];
	float hz = P.data()[8]*X + P.data()[9]*Y + P.data()[10]*Z + P.data()[11];

	float invhz = 1.0/hz;
	x = hx*invhz;
	y = hy*invhz;
}

inline void CMesh::GetPoint(int i,  float& x, float& y, float& z) {
	x=mPoints[i](0);
	y=mPoints[i](1);
	z=mPoints[i](2);
}

inline void CMesh::GetPatch(int i, int& x, int& y, int& z) {
	x=mPatch[i](0);
	y=mPatch[i](1);
	z=mPatch[i](2);   
}


inline void CMesh::GetBounds(int J, int i,  float& x, float& y, float& z) {
	x=mBounds(J,i,0);
	y=mBounds(J,i,1);
	z=mBounds(J,i,2);
}

inline void CMesh::GetPoint(int i,  float& x, float& y, float& z, int& j) {
	x=mPoints[i](0);
	y=mPoints[i](1);
	z=mPoints[i](2);
	j=int(mPoints[i](3));
}

inline void CMesh::GetPoint(int i,  float& x, float& y, float& z, float& j) {
	x=mPoints[i](0);
	y=mPoints[i](1);
	z=mPoints[i](2);
	j=mPoints[i](3);
}


#endif

