//
//
//	o_TransformationList.h	header for o_TransformationList
//
//	Oliver Grau			28-FEB-1997

#ifndef _o_TransformationList_incl
#define _o_TransformationList_incl

#include "o_Property.h"
#include "o_TransformationAccu.h"
#include <vector>
using std::vector;

/*!
  \class o_TransformationList o_TransformationList.h
  \brief Describes the accumulated transformation of a body or a surface

   The transformation is described by a homogenous matrix (\sa o_HMatrix ).
   */

class o_TransformationList : public o_Property
{

	public:
		o_TransformationList();
		~o_TransformationList();

		const char *GetClassName() const;
		void	Clear();

		o_Property *MakeInstance();
		void Copy(o_Property &org);

    //! returns the accumulated transformation list
    const std::vector<o_TransformationAccu*> & GetTransformationAccuList() const { return translist; } 
    //! adds accumulated transformation to the end of the accumulated transformation list
    void AddTransformationAccu(o_TransformationAccu & itm);
    //!  substract accumulated transformation itm from the accumulated transformation list
    bool SubTransformationAccu( o_TransformationAccu & itm );

  private:
    std::vector<o_TransformationAccu*> translist;
};

#endif
