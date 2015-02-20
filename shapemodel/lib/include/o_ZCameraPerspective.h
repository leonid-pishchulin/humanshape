#include "o_ColorCamera.h" 
#include        <math.h>

#ifndef o_ZCameraPerspective_h
#define o_ZCameraPerspective_h

/*
Contains method P_render as a subsitution for other cameras P_render methods
Performs exact, perspective calculation of Z-buffer
*/

/*!
  \class o_ZCameraPerspective.h o_ZCameraPerspective.h
  \brief ColorCamera with perspective rendering

  ZCamera is a normal ColorCamera except for the fact that the method
  P_render is substituted to create exact !(no affine transform) depth 
  maps while rendering.
  */

class o_ZCameraPerspective : public o_ColorCamera {
  void	P_render( const o_World&, o_Triangle & );   


public:
  o_ZCameraPerspective();
  ~o_ZCameraPerspective();


};
#endif

