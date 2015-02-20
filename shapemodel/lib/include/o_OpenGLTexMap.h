/*
 *
 *     oGeM_OpenGLTexMap.h   - class-definition for oGeM/GLTexMap-property
 *
 * Copyright 1997, University of Hannover,
 * Institut fuer Theoretische Nachrichtentechnik und Informationsverarbeitung
 *
 * written by O. Stahlhut, stahlhut@tnt.uni-hannover.de
 *
 */

#ifndef o_OpenGLTexMap_h
#define o_OpenGLTexMap_h

#include "o_Property.h"
#include "o_TextureMap.h"
#include "o_Surface.h"
#include "o_Body.h"

#define TEXF_OFF  0
#define TEXF_2    1
#define TEXF_4    2
#define TEXF_AUTO 3

/*!
  \class o_OpenGLTexMap o_OpenGLTexMap.h
  \brief property-class for generating and storing OpenGL-compliant texturemaps
  */

class o_OpenGLTexMap : public o_Property
{
 public:
  o_OpenGLTexMap();
  //! generates an o_OpenGLTexMap object from the oGeM-texturemap tm
  o_OpenGLTexMap(o_TextureMap *tm);

  //! allocates an OpenGL-texturemap with fixed sizes x, y
  o_OpenGLTexMap(int x, int y);

  ~o_OpenGLTexMap();

  //! returns a pointer to the picture data
  unsigned char * Picture() { return tm_p; }

  /*! converts an oGeM standard o_TextureMap to an OpenGL-compliant map 
    scale_x and scale_y are set according to the necessary scaling to
    get OpenGL texturemaps with 2^m,2^n pixel boundary lengths.
    */
  void Convert2GL(o_TextureMap *);

  /*! scales the x- and y-value of the texturebindings on surface sp
    by "scale_x" and "scale_y" [the values returned by Convert2GL()].
    If sub is TRUE, the bindings on all subsurfaces are changed, too */
  void CreateGLBindings(o_Surface *sp, float scale_x, float scale_y, int sub);
  /*!
    \overload
    */
  void CreateGLBindings(o_Body *bp, float scale_x, float scale_y);

  //! returns the x-scaling factor of the OpenGL map. \sa Convert2GL()
  double StretchX() const { return stretch_x; }
  //! returns the y-scaling factor of the OpenGL map. \sa Convert2GL()
  double StretchY() const { return stretch_y; }

  /*! selects a filtering (simple subsampling) of the OpenGL map.
    method = TEXF_OFF, turns the filtering off, TEXF_2 does a factor
    2 subsampling, TEXF_4 applies TEXF_2 two times. */
  void Filter(unsigned int method);

  //! sets the pixel (x,y) of the OpenGL map to an RGB value
  void SetPixel(int x, int y, unsigned char r_val, unsigned char g_val, unsigned char b_val);
  
  //! returns the luminance value of pixel(x,y), band
  unsigned char GetPixel(int x, int y, int band);

  //! returns the x-size of the map
  unsigned int Width() const { return dx; }
  //! returns the y-size of the map
  unsigned int Height() const { return dy; }

  const char *GetClassName() const { return "o_OpenGLTexMap"; }
  virtual o_Property *MakeInstance() { return new o_OpenGLTexMap; }

  //! returns true if texture has be updated
  bool GetUpdatedFlag() { return updated; }
  
  //! sets or clears the update flag
  void SetUpdatedFlag(bool up) { updated=up; }
  
 private:
  /// pointer to the picture memory
  unsigned char *tm_p;
  unsigned int dx, dy;
  double stretch_x, stretch_y;
  bool updated;
};


/*! 
  \relates o_OpenGLTexMap

  CreateOpenGLMaps() is a function that generates OpenGL compliant
  texture maps from o_TextureMap properties. filter selects a 
  (pre) filtering as described by o_OpenGLTexMap::Filter() 
  */

bool CreateOpenGLMaps( o_Body *bp, unsigned int filter=TEXF_OFF );

#endif
