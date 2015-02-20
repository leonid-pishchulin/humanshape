#ifndef o_HashTable_incl
#define o_HashTable_incl


//
//
//	o_HashTable.h	header for o_HashTable
//
//	Manfred Ribitzki	06-APR-1995
//
//
// In this hashtable the reference to both the keys and the values is stored,
// but the data pointed to is evaluated for searching. The size of the memory 
// region that is evaluated is determined via the sizeof operator. If the key
// fits into a memory region of a long, you should use o_HashTableSmpl.

#include "dynhash.h"

/*!
  \class o_HashTable o_HashTable.h
  \brief Template class of a hashtable.

  This hashtable can not be used with keys
  derived from o_MetaClass, because these keys might change during life time
  (the internal reference counter might change for example). You can use
  o_HashTableSmpl with a pointer to the object as the key instead.
  */

template<class K, class V>
class o_HashTable
{
private:
  //do not copy

  o_HashTable (const o_HashTable<K, V> &);	//there are no bodys!
  o_HashTable &operator= (const o_HashTable<K, V> &);

protected:
  HashTable *m_pht;

public:
  /*! o_HastTable\<K, V\> creates a HashTable with data type K as the key and data type V as the value. Note that on destruction the entered elements
    are not deleted.*/
  o_HashTable();

  virtual ~o_HashTable();
  
  //void Clear();

  //! Removes the element with the given key. The element is not deleted.
  void RemoveElement (K &Key);

  /*! Returns the reference to the corresponding value if Key
    exists in the HashTable, otherwise NULL. */
  V *Lookup (K &Key);

  /*! Retrieves the references to the key and the value of the first
    element in the HashTable. If there is no element in the 
    HashTable the function returns 0. */
  int GetFirstElement (K **ppKey, V **ppVal);
  
  /*! Retrieves the next element of the HashTable. This function
    may only be called when preceded by a call to
    GetFirstElement() or GetNextElement() . Do not 
    call any update functions between subsequent calls of these
    functions. If there is no element left in the HashTable the
    function returns 0. */
  int GetNextElement (K **ppKey, V **ppVal);

  //! Returns the number of elements in the HashTable.
  unsigned long GetKeyCount ();

  /*! Enters the given element into the HashTable. If an element with
    the given key already exists in the HashTable it will be
    replaced. Note that the references to the given Key and Value
    are stored. You have to make sure that the entered elements are
    valid as long as they are stored in the HashTable.  A return
    value of 0 indicates an error. */
  int EnterElement (K &Key, V &Val);

  //! Dumps statistics to pf.
  void Dump (FILE *);

};


//
//
//	methods for class o_HashTable
//
//


template<class K, class V> 
inline o_HashTable<K, V>::o_HashTable()
{
  m_pht = ::CreateHashTable (RefMemoryValue, o_True, RefMemoryValue, o_True);
}


template<class K, class V>
inline o_HashTable<K, V>::~o_HashTable()
{
  if (m_pht)
    ::DestroyHashTable (m_pht);
}

/*template<class K, class V>
inline void o_HashTable<K, V>::Clear()
{
  if (m_pht)
    {
      ::DestroyHashTable (m_pht);
      m_pht = NULL;
    }
}*/


template<class K, class V>
inline unsigned long o_HashTable<K, V>::GetKeyCount ()
{
    return ::GetHashKeyCount (m_pht);
}

template<class K, class V>
inline V *o_HashTable<K, V>::Lookup (K &Key)
{
  Value k;
  RefMemVal rmv;
  HashElement *phe;

  rmv.Leng = sizeof (K);
  rmv.Data = (unsigned char *)&Key;
  k.RefMem = &rmv;

  phe = ::LookupHashElement (m_pht, k);
  if (phe)
    return (V *)phe->Data.RefMem->Data;
  else
    return NULL;
}



template<class K, class V>
inline void o_HashTable<K, V>::RemoveElement (K &Key)
{
  Value k;
  RefMemVal rmv;

  rmv.Leng = sizeof (K);
  rmv.Data = (unsigned char *) &Key;
  k.RefMem = &rmv;

  ::RemoveHashElement (m_pht, k);
}

template<class K, class V>
inline int o_HashTable<K, V>::GetFirstElement (K **ppKey, V **ppVal)
{
  HashElement *phe;
  phe = ::GetFirstHashElement (m_pht);
  if (phe)
    {
      *ppKey = (K *)phe->Key.RefMem->Data;
      *ppVal = (V *)phe->Data.RefMem->Data;
      return 1;		//o.k.
    }
  else
    {
      ppKey = NULL;
      ppVal = NULL;
      return 0;
    }
}

template<class K, class V>
inline int o_HashTable<K, V>::GetNextElement (K **ppKey, V **ppVal)
{
  HashElement *phe;
  phe = ::GetNextHashElement (m_pht);
  if (phe)
    {
      *ppKey = (K *)phe->Key.RefMem->Data;
      *ppVal = (V *)phe->Data.RefMem->Data;
      return 1;		//o.k.
    }
  else
    {
      ppKey = NULL;
      ppVal = NULL;
      return 0;
    }
}      

template<class K, class V>
inline int o_HashTable<K, V>::EnterElement (K &Key, V &Val)
{
  Value k, v;
  RefMemVal krmv, vrmv;

  krmv.Leng = sizeof (K);
  krmv.Data = (unsigned char *) &Key;
  vrmv.Leng = sizeof (V);
  vrmv.Data = (unsigned char *) &Val;
  k.RefMem = &krmv;
  v.RefMem = &vrmv;

  if (::EnterHashElement (m_pht, k, v))
    return 1;
  else
    return 0;
}

template<class K, class V>
inline void o_HashTable<K, V>::Dump (FILE *fp)
{
  fDumpHashTable (0, m_pht, fp);
}




// This is a special version of the hashtable where the key is small enough to 
// fit into one long. In this case the key itself is stored into the hashtable.

/*!
  \class o_HashTableSmpl o_HashTable.h
  \brief Hash table with keys of type unsigned long.

  Template class of a hashtable with a key which can be safely casted into an
  unsigned long. You have to make sure, that you are using a proper key.
  */

template<class K, class V>
class o_HashTableSmpl
{
private:
  //do not copy
  o_HashTableSmpl (const o_HashTable<K, V> &);		//there are no bodys
  o_HashTableSmpl &operator= (const o_HashTable<K, V> &);

protected: 
  HashTable *m_pht;

public:
        /*! o_HashTableSmpl\<K, V\> creates a HashTable with data type K
	  as the key and data type
	  V as the value. If you don´t use a proper key, i.e. it cannot
	  be converted into unsigned long, you will probably get a 
	  compiler error. Use o_HashTable instead. */
	o_HashTableSmpl() {
	  m_pht = ::CreateHashTable (SimpleValue, o_False, RefMemoryValue, o_True);
	}

	virtual ~o_HashTableSmpl() {
	  if (m_pht)
	    ::DestroyHashTable (m_pht);
	}
	
	//! Returns the number of elements in the HashTable.
	unsigned long GetKeyCount () {
	  return ::GetHashKeyCount (m_pht);
	}

	/*! Returns the reference to the corresponding value if Key
	  exists in the HashTable, otherwise NULL. */
	V *Lookup (K Key) {
	  Value k;
	  HashElement *phe;
  
	  k.Smpl = (unsigned long) Key;
	  phe = ::LookupHashElement (m_pht, k);
	  if (phe)
	    return (V *)phe->Data.RefMem->Data;
	  else
	    return NULL;
	}

	/*! Removes the element with the given key. The corresponding 
	  value itself is not deleted. */
	void RemoveElement (K Key) {
	  Value k;

	  k.Smpl = (unsigned long) Key;
	  ::RemoveHashElement (m_pht, k);
	}

	/*! Retrieves the key and the reference to the value of 
	  the first element in the HashTable. If there is no element 
	  in the HashTable the function returns 0. */
	int GetFirstElement (K *pKey, V **ppVal)
	{
	  int rc;
	  HashElement *phe;
	  phe = ::GetFirstHashElement (m_pht);
	  if (phe) {
	      *pKey = (K)phe->Key.Smpl;
	      *ppVal = (V *)phe->Data.RefMem->Data;
	      rc = 1;		//o.k.	       
	  } else {
	      *pKey = 0L;
	      *ppVal = NULL;
	      rc = 0;
	  }
	  return rc;
	}
	
	/*! Retrieves the next element of the HashTable. This function
	  may only be called when preceded by a call to
	  GetFirstElement() or GetNextElement() . Do not 
	  call any update functions between subsequent calls of these
	  functions. If there is no element left in the HashTable the
	  function returns 0. */
	int GetNextElement (K *pKey, V **ppVal)
	{
	  int rc;
	  HashElement *phe;
	  phe = ::GetNextHashElement (m_pht);
	  if (phe) {
	      *pKey = (K)phe->Key.Smpl;
	      *ppVal = (V *)phe->Data.RefMem->Data;
	      rc = 1;		//o.k.	       
	    } else {
	      *pKey = 0L;
	      *ppVal = NULL;
	      rc = 0;
	    }
	  return rc;
	}

	/*! Enters the given element into the HashTable. If an element with
	  the given key already exists in the HashTable it will be
	  replaced. Note that the 
	  reference to the given Value is stored. You have to
	  make sure that the entered values are valid as long as they
	  are stored in the HashTable.
	  A return value of 0 indicates an error. */
	int EnterElement (K Key, V &Val)
	{
	  Value k, v;
	  RefMemVal vrmv;

	  k.Smpl = (unsigned long) Key;
	  vrmv.Leng = sizeof (V);
	  vrmv.Data = (unsigned char *) &Val;
	  v.RefMem = &vrmv;
    
	  int rc = (::EnterHashElement (m_pht, k, v)) ? 1 : 0;
	  return rc;
	}

	//! Dumps statistics to pf.
	void Dump (FILE *fp)
	{
	  fDumpHashTable (0, m_pht, fp);
	}
      

};

#endif
