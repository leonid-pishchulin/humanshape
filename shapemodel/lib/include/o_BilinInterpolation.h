#ifndef	_o_BilinInterpolation_incl
#  define	_o_BilinInterpolation_incl
//
//	o_BilinInterpolation.h	- Header for oGeM
//
//	Oliver Grau, Jan./Mar.1993
//


#include	"o_BilinInterpolation.h"
#include	"o_Picture.h"

#define GRID_POINT_RANGE 1e-12

/*!
  \relates o_Picture o_BilinInterpolation.h

  Interpolate in image with mask image control. Mask image splits the image 
  into regions. Only pixel in the same image region are used for
  interpolation. Pixel with a predefined exclude_value are excluded from
  interpolation.
  
  Parameter:
  <UL>
  <LI> double ix, double iy : Image point coordinates in PCS</LI>
  <LI> const o_Picture<float> * org_img : o_Picture used for interpolation </LI>
  <LI> const Picture<short> * mask_img (optional) : o_Picture used as mask
  for conditional interpolation. Only Pixel with same value in mask_img
  are used for interpolation </LI>
  <LI> const double *exclude_val (optional) : All Pixel in org_img that
  have the value "*exclude_val" are not considered during interpolation. </LI>
  </UL>
  */

double BilinInterpolation(double ix, double iy, const o_Picture<float> *org_img, 
const o_Picture<short> * mask_img = NULL, const double *exclude_val = NULL);

#endif


