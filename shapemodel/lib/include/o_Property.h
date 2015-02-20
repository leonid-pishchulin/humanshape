#ifndef	_incl_property
#define	_incl_property
//
//      o_Property.h  - Header for oGeM property
//
//      Oliver Grau, 28-OCT-1993
//

#include	"oGeM.h"

class	o_DataBaseOut;
class	o_DataBaseIn;
class   o_Triangle;

/*!
  \class o_Property o_Property.h
  \brief meta class for properties

  A property is a soft attribute that can assigned to a scene object
  at the run time. A property object should assigned exactly to one
  scene object using the interface of o_PropertyList
  */

class	o_HMatrix;
class	o_PropertyList;

class	o_Property : public o_MetaClass     {
	public:
		o_Property() ;
		~o_Property();
		
		/*! get name of the property. The default (if not specified) is
		  GetClassName() */
		virtual	const char	*PropertyName() const;

		/*! this method is empty for this base class. Some methods like
		o_Point::ApplyTransform() call for all assigned properties
		objects of the point object this method. So for speciall
		derived classes this method can be overwritten */
		virtual void ApplyTransform(o_HMatrix &matrix); 

		/*! create a new instance of the objects class. Must be
		  specified for each child class */
		virtual o_Property *MakeInstance()=0;

		virtual void	Copy(o_Property &org);
		int     CheckForField( o_DataBaseIn * dbptr, const char *name );

		/*! indicates a value update. The ChangedFlag is evaluated
		by update functions and after a successfull update the
		flag is resetted. \sa o_StreamOut . The method
		should be called after any value update of the property. */
		void	SetChangedFlag();

		//! get value of the changed flag.
		int	GetChangedFlag() const;

		//! returns pointer to parent object or 0 if no parent present
		o_MetaClass	*GetParentRef() const;

		/*!  returns an infostring describing property-values. Possible
		modes are 0: long infostring, 1: short infostring */
		virtual const char *GetInfostring(int mode);

	protected:
		friend class o_StreamOut;
		friend class o_DataBaseOut;
		friend class o_Triangle;
		void	UnsetChangedFlag();

		friend class o_PropertyList;
		o_MetaClass	*parent;	// points to parent object
		void	SetParent(o_MetaClass *par);
		void    StoreInfoString(const char *str);
		const char *RetrieveInfoString(void);
	private:
		int	changed_flag;
	        char *infostring;
};

#endif
