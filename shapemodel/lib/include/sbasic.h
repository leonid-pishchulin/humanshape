
/* SCCS-ID: $Id: sbasic.h 431 2007-09-26 00:26:07Z rhys $ (c) Oliver Grau */



#ifndef	_SBASIC_H
#  define	_SBASIC_H


#include <iostream>
#include <stdio.h>
#include <math.h>

//extern	"C" void	exit(int);

# include	"o_Basis.h"
# include "o_Notify.h"
# define        error_handler(l,m)      { \
			o_Notify(O_NFY_WARN,"%s line:%d :%s",__FILE__,__LINE__ , m); }

#ifndef nil
# define        nil     (0)
#endif
#ifndef true
# define        true    (1)
#endif
#ifndef false
# define        false   (0)
#endif

//------------------------------------------------------------------------------// Memory Management
//------------------------------------------------------------------------------
#ifdef	notdef
struct  memory_elem_type { memory_elem_type* next; };
typedef memory_elem_type* memory_elem_ptr;

extern memory_elem_ptr memory_free_list[];

extern memory_elem_ptr memory_allocate_block(int);
extern memory_elem_ptr allocate_bytes(int);
extern memory_elem_ptr allocate_words(int);

extern void            deallocate_bytes(void*,int);
extern void            deallocate_bytes_with_check(void*,int);
extern void            deallocate_words(void*,int);
extern void            memory_clear();
extern void            memory_kill();
extern void            print_statistics();

inline void deallocate_list(void* head,void* tail, int bytes)
{ memory_elem_ptr(tail)->next = memory_free_list[bytes];
  memory_free_list[bytes] = memory_elem_ptr(head);
}


#define OPERATOR_NEW(bytes)\
void* operator new(size_t)\
{ memory_elem_ptr* q = memory_free_list+bytes;\
  if (*q==0) *q = memory_allocate_block(bytes);\
  memory_elem_ptr p = *q;\
  *q = p->next;\
  return p;\
 }

#define OPERATOR_DEL(bytes)\
void  operator delete(void* p)\
{ memory_elem_ptr* q = memory_free_list+bytes;\
  memory_elem_ptr(p)->next = *q;\
  *q = memory_elem_ptr(p);\
 }

#define OPERATOR_DEL_WITH_CHECK(bytes)\
void  operator delete(void* p) { deallocate_bytes_with_check(p,bytes); }


#define LEDA_MEMORY(type)\
OPERATOR_NEW(sizeof(type))\
OPERATOR_DEL(sizeof(type))

#endif
#endif	/* _SBASIC_H */

