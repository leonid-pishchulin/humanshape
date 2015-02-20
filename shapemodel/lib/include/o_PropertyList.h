
#ifndef	_o_PropertyList_incl
#  define	_o_PropertyList_incl
//
//	o_PropertyList.h	- abstract class for holding properties
//
//	Oliver Grau, 19-MAY-1994
//


#include	"oGeM.h"
#include <list>
using std::list;

//
//      Forward declaration of classes
//
class   o_Property;
class o_DataBaseOut;
class o_StreamOut;

/*!
  \class o_PropertyList o_PropertyList.h
  \brief mixin class for keeping properties

  o_PropertyList is an abstract class do not instantiate. Note that
  on deletetion	the destructor deletes all appended property objectes
*/

class	o_HMatrix;
class	o_MetaClass;

class	o_PropertyList : public o_MetaClass {
	public:
                /*! get surface property 'propertyname'.
		If no property with the name is found it returns NULL.*/
		o_Property *GetLocalProperty( char *propertyname );

		//! adds a property to the property list.
		void AddProperty( o_Property &property );

		//! strips a property from the property list.
		bool SubProperty( o_Property &property );

		/*! get a list of all properties of the surface - including the
		inherited. The list should be deleted with the delete operator */
    std::list<o_Property*> *GetPropertyList( char *propertyname );

		//! calls ApplyTransform for every included property object

		virtual void ApplyTransformToProperties(o_HMatrix &matrix);

		/*! checks if any o_Property in the list has its
		  changed flag set */
		int UpdatePropertylist();

		void    Copy( o_PropertyList &org );

	protected:
		o_PropertyList();
		virtual ~o_PropertyList( );
	private:
		friend class o_DataBaseOut;
		friend class o_StreamOut;
    std::list<o_Property*> proplist;
};


#endif
