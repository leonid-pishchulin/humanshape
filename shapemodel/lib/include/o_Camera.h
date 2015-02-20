#ifndef	_o_realcamera
# define	_o_realcamera

//
//	o_Camera.h	- Camera class definitions
//
//	Oliver Grau, Jan./Feb.1993
//	Reinhard Koch, 28.4.94


#include "o_Picture.h"
//#include "CPicture.h"
#include "o_World.h"
#include "o_MetaCamera.h"

/*!
  \class o_Camera o_Camera.h
  \brief class for real world camera
  */

class	o_Camera : public o_MetaCamera {
 public:
  o_Camera( );
  ~o_Camera();
  virtual const   char    *GetClassName() const;

  /*! returns the width of the o_CPicture , if the camera 
    picture is set, otherwise returns the dimensions of the 
    camera plane. */
  virtual int     Width() const;
  //! returns the height of the o_CPicture.
  virtual int     Height() const;

  //! returns camera picture
  o_CPicture<unsigned char>	&GetPicture();

  //! return const camera picture
  const o_CPicture<unsigned char>     &GetPicture() const;

  //! returns monochrome picture computed from camera picture
  o_Picture<unsigned char>	&GetLuminancePicture();

  //! returns gradient in x from monochrome picture of Camera
  o_Picture<float>	&GetXgradientPicture();

  //! returns gradient in y from monochrome picture of Camera
  o_Picture<float>	&GetYgradientPicture();


  //! set camera picture
  void		SetPicture(o_CPicture<unsigned char>	&pic);

  /*! set camera picture. In oposite to SetPicture only a 
    reference to pic is stored. The picture object is not 
    deleted by the destructor! */
  void		SetPictureReference(o_CPicture<unsigned char>	&pic);

 protected:
  o_CPicture<unsigned char>	*pic;
  int	extern_allocated;		// if true *pic is not deleted by the
					// destructor
  o_Picture<unsigned char> *lumpic;
  o_Picture<float> *Xgradpic;
  o_Picture<float> *Ygradpic;


 private:
  void ConvertToLum();
  void Xgradient();
  void Ygradient();

};


#endif


