//
//
//	o_TransformationAccu.h	header for o_TransformationAccu
//
//	Oliver Grau			28-FEB-1997

#ifndef _o_TransformationAccu_incl
#define _o_TransformationAccu_incl

#include "o_Property.h"
#include "o_HMatrix.h"
#include "o_Surface.h"

/*!
  \class o_TransformationAccu o_TransformationAccu.h
  \brief Describes the accumulated transformation of a body or a surface

   The transformation is described by a homogenous matrix (\sa o_HMatrix ).
   */

class o_TransformationAccu : public o_Property
{
	private:
		o_HMatrix m_trans;
	public:
		/*! initialized with the transformation of the next parent
		object. Note that the surface might share common
		points with adjoining surfaces, which will be affected by this
		operation. It is up to you to take care, that this operation
		does not hurt. Never create transformation objects for
		surfaces, which are rigid connected or having the 
		estimation_enable_flag (\sa o_EstimatorDescriptor ) not
		set. */
		o_TransformationAccu();

		~o_TransformationAccu();

		const char *GetClassName() const;

		/*! Gets the transformation matrix. The matrix describes the 
		absolute transformation of the associated object. */
		o_HMatrix &GetTransform() ;

		/*! Applies transformation to this transformation object and
		to the transformation objects of all subsurfaces of the 
		associated object */
		void ApplyTransform (o_HMatrix &);

		o_Property *MakeInstance();
		void Copy(o_Property &org);
};

#endif
