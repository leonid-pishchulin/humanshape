/*
 *
 * oGeM_OpenGLTexBind.h   - class-definition for oGeM/GLTexBinding-property
 *
 * Copyright 1997, University of Hannover,
 * Institut fuer Theoretische Nachrichtentechnik und Informationsverarbeitung
 *
 * written by O. Stahlhut, stahlhut@tnt.uni-hannover.de
 *
 */

#ifndef o_OpenGLTexBind_h
#define o_OpenGLTexBind_h

#include "o_Property.h"
#include "o_TextureBinding.h"
#include "o_Surface.h"
#include "o_Body.h"

/*!
  \class o_OpenGLTexBind o_OpenGLTexBind.h
  \brief property-class for generating and storing texturebindings for use with oGeM_OpenGLTexMap (OpenGL-compliant texturemaps)
  */

class o_OpenGLTexBind : public o_Property
{
 public:
  o_OpenGLTexBind() {}
  //! generates an o_OpenGLTexBind object from the oGeM-texturebinding tb
  o_OpenGLTexBind(o_TextureBinding *tb);

  ~o_OpenGLTexBind();

  int NoOfCoord() { return num_bindings; }
  void GetCoordinates(int no, double &x, double &y);
  double GetCoordinateX( int no );
  double GetCoordinateY( int no );
  void SetCoordinates(int no, double x, double y);
  void UpdateCoordinates(o_TextureBinding *tb);

  const char *GetClassName() const { return "o_OpenGLTexBind"; }
  virtual o_Property *MakeInstance() { return new o_OpenGLTexBind; }

 private:
  double *bindings;
  int num_bindings;
};

#endif
