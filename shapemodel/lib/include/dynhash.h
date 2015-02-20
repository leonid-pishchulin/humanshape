/* This code is copyright by
 *     Elmar Bartel 1993,
 *     Institut fuer Informatik, TU Muenchen, Germany  
 *     bartel@informatik.tu-muenchen.de
 * You may use this code as you wish, as long as these
 * comment-lines are keept intact and in place.
 *
 * Its an implementation of dynamic hashing, so that the hashtable
 * grows and shrinks with the number of elements stored in the table.
 * The algorithms used in this package are straight forward
 * copied from:
 * Dynamic hashing, by Per-Ake Larson, in CACM April 1988 pp 446-457.
 * There are other implementations of these algorithms around, but
 * none provided a flexible and uniform handling of different data
 * types used for keys or data in the hashtable.
 *
 * This implementation should fit in a lot of applications.
 * Each element of the hash table (HashElement) is a tuple
 * of a key and associated data.
 * The key and the data values can be one out of four 
 * choices (HashTableValues):
 *
 *     SimpleValue     fits in one long
 *     String          Pointer to a usual C-String
 *     CompactMemVal   Pointer to a memory region
 *                     prefixed with a long-value
 *                     giving the length of the region
 *     MemVal          Pointer to a struct, containing
 *                     The length of data and the pointer
 *                     to the data
 *
 * On creation of the hash table the type of key
 * and and the type of data is given.
 * Also the handling of the hash elements, whether the insert
 * routine should copy the key or data or both.
 *
 * It is also possible to deliver all key/data pairs one after
 * another. See GetFirstHashElement and GetNextHashElement.
 *
 * The name of the procedure should explain themselves.
 * A small demo can be seen on a unix host, if there
 * is a /etc/hosts file. Compile hash.c with -DMAIN!
 */


/* M. Ribitzki 18.11.1994
 * New Value type added: RefMemVal
 * Bug fixed on iterating.
 */


#ifndef DYNHASH_INCLUDED
#define DYNHASH_INCLUDED

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef GDATA
#define DATA
#else
#define DATA extern
#endif

/* The Boolean type -- what else */
typedef enum {
       o_False,
       o_True
} Boolean;

/* This is returned from some routines to report
 * success or failure of their task
 */
typedef enum {
       OkStatus,
       FailedStatus
} o_Status;

/* Some simple types for all modules */
typedef unsigned long   uLong;
typedef unsigned short  uShort;
typedef unsigned char   uChar;
typedef char *          String;
#define NilString       (String)0
#define FreeString(s)   free(s)

/* Here we define the value types used throughout the source */
typedef enum {
       SimpleValue,
       StringValue,
       CompactMemoryValue,
       MemoryValue,
       RefMemoryValue,
       ValueTypeCnt
} ValueType;
DATA char *ValueTypeToString[ValueTypeCnt]
#ifdef GDATA
       = { "SimpleValue", "StringValue", "CompactMemoryValue", "MemoryValue",
           "RefMemoryValue" }
#endif
;

/* This is the counted bytestring.
 * It is allocated via one malloc call.
 */
typedef struct CompactMemVal {
       uLong   Leng;
       uChar   Data[1];
} CompactMemVal;
#define NilCompactMemVal        (CompactMemVal *)0
#define NewCompactMemVal(n)     (CompactMemVal *)malloc(sizeof(CompactMemVal)+n*sizeof(uChar))
#define FreeCompactMemVal(v)    free(v)

/* This is nearly the same as the CompactMemVal
 * but you will need two malloc-calls to duplicate
 * this. But can be more efficient, if no data in the
 * hash-table is needed.
 */
typedef struct MemVal {
       uLong   Leng;
       uChar   *Data;
} MemVal;

#define NilMemVal       (MemVal *)0
#define NewMemVal       (MemVal *)malloc(sizeof(MemVal))
#define FreeMemVal(v)   free(v)
#define DeleteMemVal(v) if (v!=NilMemVal) free(v->Data), free(v)


/* This looks like MemVal but instead of copying the data referenced by
 * Data into the HashTable the reference itself is copied.
 */
typedef struct RefMemVal {
       uLong   Leng;
       uChar   *Data;
} RefMemVal;

#define NilRefMemVal	(RefMemVal *)0
#define NewRefMemVal	(RefMemVal *)malloc(sizeof(RefMemVal))
#define FreeRefMemVal(v) free (v)
#define DeleteRefMemVal(v) if (v!=NilRefMemVal) free(v)


typedef union {
       void            *Ptr;           /* SimpleValue */
       uLong           Smpl;           /* SimpleValue */
       String          Str;            /* StringValue */
       CompactMemVal   *CompactMem;
       MemVal          *Mem;
       RefMemVal       *RefMem;
} Value;
#endif



typedef struct {
       Value   Key;
       Value   Data;
} HashElement;
#define NilHashElement  (HashElement *)0

struct HashTable;
#define NilHashTable    (struct HashTable *)0

struct HashTable *CreateHashTable     (ValueType KeyTyp, Boolean CopyKey,
                                       ValueType DataTyp, Boolean CopyData);
void             DestroyHashTable     (struct HashTable *);
HashElement      *EnterHashElement    (struct HashTable *, Value key, Value data);
void             RemoveHashElement    (struct HashTable *, Value key);
HashElement      *LookupHashElement   (struct HashTable *, Value key);
HashElement      *GetFirstHashElement (struct HashTable *);
HashElement      *GetNextHashElement  (struct HashTable *);
void             fDumpHashTable       (int ind, struct HashTable *, FILE *);
uLong		 GetHashKeyCount      (struct HashTable *);


#ifdef __cplusplus
}
#endif


