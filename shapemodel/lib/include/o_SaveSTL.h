//
// o_SaveSTL.h
//
// Copyright (C) 2006 University of Adelaide
//
// Permission is granted to any individual or institution to use, copy, modify,
// and distribute this software, provided that this complete copyright and
// permission notice is maintained, intact, in all copies and supporting
// documentation.
//
// University of Hannover provides this software "as is" without
// express or implied warranty.
//
//
// History:
//
// Created: T. Thormaehlen  April 2006      Initial design and implementation.
// Updated:

#ifndef _o_SaveSTL_h
#define _o_SaveSTL_h

#include	"fileio.h"
#include	"o_World.h"
#include	"o_SaveBody.h"


/*!
  \class o_SaveSTL o_SaveSTL.h
  \brief save Body to STL-file (Stereolithography)
  */

class	o_SaveSTL : public o_SaveBody {
  
  public:
   //! creates STL-file object
   o_SaveSTL( );

   int SaveWorld (o_World &w);
   int Save( o_Body &body );

   virtual void Save_All();
   virtual void Save_Pointlist();
   virtual void Save_Polygons();
  
};

#endif

