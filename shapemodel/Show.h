#ifndef ShowH
#define ShowH

#include "CMatrix.h"
#include "CTensor.h"
#include <vector>



namespace NShow {

  // Current view no in image sequence
  extern int mViewNo;
  // Current image no in image sequence
  extern int mImageNo;
  // Current input images
  extern std::vector<CTensor<float> >* mInputImage;
  // Directory into which the results are written
  extern std::string mResultDir;
  // Directory where the input data can be found
  extern std::string mInputDir;
  // File name without ending
  extern std::string mFileName;
  // cout stream for Windows
  extern std::ostringstream mOut;
  // Logfile
  extern std::ofstream mLogFile;

  extern CMatrix<std::string> mInput;

  extern bool mForceContourColor;
  extern int mContourColorR;
  extern int mContourColorG;
  extern int mContourColorB;
  extern int mMeshColorR;
  extern int mMeshColorG;
  extern int mMeshColorB;
  extern int mMeshPointColorR;
  extern int mMeshPointColorG;
  extern int mMeshPointColorB;  

  // Allow window update
  void idle();
  // Show message in status bar
  void message(bool aClear = true);
  // Show image in popup window
  // void show(int aNo, CMatrix<float>& aImage, const char* aText = 0);
  // void show(int aNo, CTensor<float>& aImage, const char* aText = 0);
  // Write image to file in result directory
  // void writeToFile(CMatrix<float>& aImage, const char* aFileName, const char* type = "png");
  // void writeToFile(CTensor<float>& aImage, const char* aFileName, const char* type = "png");
  // void writeToFile(CMatrix<float>& aImage, const char* aFileName, int i, const char* type = "png");
  // void writeToFile(CTensor<float>& aImage, const char* aFileName, int i, const char* type = "png");
  // void writeToFile(std::vector<CMatrix<float> >& aImage, const char* aFileName, const char* type = "png");
  // void writeToFile(std::vector<CTensor<float> >& aImage, const char* aFileName, const char* type = "png");
  // void writeToFile(std::vector<CMatrix<float> >& aImage, const char* aFileName, int i, const char* type = "png");
  // void writeToFile(std::vector<CTensor<float> >& aImage, const char* aFileName, int i, const char* type = "png");
  
}

std::ostream& operator<<(std::ostream& aStream, const std::vector<double>& vec);
std::ostream& operator<<(std::ostream& aStream, const std::vector<float>& vec);
std::ostream& operator<<(std::ostream& aStream, const std::vector<int>& vec);
std::ostream& operator<<(std::ostream& aStream, const std::vector<bool>& vec);
// std::ostream& operator<<(std::ostream& aStream, const std::vector<GLubyte>& vec);

#endif

