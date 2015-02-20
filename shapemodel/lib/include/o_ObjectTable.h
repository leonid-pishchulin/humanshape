#ifndef _o_ObjectTable_incl
#  define       _o_ObjectTable_incl
//
//      o_ObjectTable.h     - Header for o_ObjectTable class
//
//      Oliver Grau, Nov.1996
//

#include        "oGeM.h"
#include        "o_HashTable.h"

/*!
  \class o_ObjectTable o_ObjectTable.h
  \brief hash table for objects and ident handling 
  */

class o_StreamOut;
class o_StreamIn;
class o_MetaClass;

class   o_ObjectTable {
	public:
		o_ObjectTable();
		~o_ObjectTable();
		UniqueIdentifier LastIdent();
		void	IncrementIdent();
	private:
		friend class o_StreamOut;
		friend class o_StreamIn;
		friend class o_MetaClass;

		void SetLastIdent(UniqueIdentifier newid);
		o_HashTableSmpl<UniqueIdentifier,o_MetaClass>	objtable;
		UniqueIdentifier	lastid;
};
 
#endif
  

