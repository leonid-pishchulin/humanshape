//
//
//	o_ObjectTransformation.h	header for o_ObjectTransformation
//
//	Manfred Ribitzki		23-JAN-1996
//	Oliver Grau			28-FEB-1997

#ifndef _o_ObjectTransformation_incl
#define _o_ObjectTransformation_incl

#include "o_Property.h"
#include "o_HMatrix.h"
#include "o_Surface.h"


/*!
  \class o_ObjectTransformation o_ObjectTransformation.h
  \brief Describes the transformation of a body or a surface

  The transformation is described by a homogenous matrix (\sa o_HMatrix ).
  Instances of this class must be instantiated dynamically via operator new,
  as for all classes derived from o_Property . You must not add instances
  of this class to a body or a surface (or any other properlist) by yourself.
  The constructor will do this job for you. That means, if you want to
  add an instance of o_ObjectTransformation to my_body for example, just 
  write:

  new o_ObjectTransformation (my_body);
  */

class o_ObjectTransformation : public o_Property
{
private:
  //list of subsurfaces of the associated object
  const std::list<o_Surface*> *m_psurf_list;

  o_HMatrix m_trans;

  void ApplyObTransform (o_Surface &, const o_HMatrix &);

 public:
    o_ObjectTransformation();

    // old constructors for compatibility
    /*! Creates a transformation object for body and adds it to body
      as a property.*/
    o_ObjectTransformation(o_Body &);

    /*! Creates a transformation object for surface and adds it to
      surface as a property. The transformation object will be
      initialized with the transformation of the next parent
      object. Note that the surface might share common
      points with adjoining surfaces, which will be affected by this
      operation. It is up to you to take care, that this operation
      does not hurt. Never create transformation objects for
      surfaces, which are rigid connected or having the 
      estimation_enable_flag (\sa o_EstimatorDescriptor ) not set. */
    o_ObjectTransformation(o_Surface &);

    ~o_ObjectTransformation();

    const char *GetClassName() const;
    
    /*! Gets the transformation matrix. The matrix describes the 
      absolute transformation of the associated object. */
    const o_HMatrix &GetTransform() const;
    
    /*! Applies transformation to this transformation object and
      to the transformation objects of all subsurfaces of the
      associated object. */
    void ApplyTransform (const o_HMatrix &);

    o_Property *MakeInstance();
    void Copy(o_Property &org);

};

#endif
