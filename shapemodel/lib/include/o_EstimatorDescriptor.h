//
//
//	o_EstimatorDescriptor.h	header for o_EstimatorDescriptor
//
//	Manfred Ribitzki		23-JAN-1996

#ifndef _o_EstimatorDescriptor_incl
#define _o_EstimatorDescriptor_incl

#include "o_Property.h"

/*!
  \class o_EstimatorDescriptor o_EstimatorDescriptor.h
  \brief Property class used by the pose estimation classes

  If an instance of this class is found to be a local property of a surface
  and the estimation_enable_flag is not set,
  the pose estimation classes (\sa o_PoseBase and its descendants) will 
  skip this surface while setting up the parameter equations.
  */

class o_EstimatorDescriptor : public o_Property
{
private:
  int m_estimation_enable_flag;

public:
  //! Creates the EstimatorDescriptor.
  o_EstimatorDescriptor();

  ~o_EstimatorDescriptor();

  o_Property *MakeInstance();
  void    Copy(o_Property &org);

  const char *GetClassName() const;
  
  //! Indicates whether motion estimation is enabled or not.
  int GetEstimationEnableFlag() const;

  //! Sets the estimation_enable_flag.
  void SetEstimationEnableFlag (int);
};


#endif
