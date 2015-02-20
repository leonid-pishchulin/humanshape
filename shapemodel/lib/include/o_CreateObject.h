#ifndef __OGEM_CREATEOBJECT_H
#define __OGEM_CREATEOBJECT_H

//	Create Property Objects dynamicly - for restoration
//
//      O.Grau OCT-1996

#include	"o_Material.h"
#include	"o_TextureMap.h"
#include	"o_TextureBinding.h"


o_MetaClass * CreateMaterial() ;
o_MetaClass * CreateTextureMap();
o_MetaClass * CreateTextureBinding();

typedef	o_MetaClass* (*CreFuncPtr) ();

class	o_CreationLookUp {
	public:
		o_CreationLookUp(const char *clasn, CreFuncPtr );
		~o_CreationLookUp();
		char *classname;
		CreFuncPtr funcptr;
};

void	AddCreationEntry(char *clasn, CreFuncPtr fp);

void	InitCreationTable();

int	SearchObjectClass(const char *cn) ;
o_MetaClass* CreateObject(const char *cn) ;


#endif
