#ifndef	_incl_map
#define	_incl_map
//
//      o_TextureMap.h  - Header for oGeM maps
//
//      Oliver Grau, 15-JUL-1993
//      Ralf Toenjes 11.10.94 Min/Mag-TexFilter

#include	"oGeM.h"
#include	"o_Picture.h"
#include	"o_Property.h"
#include        "o_AffineMap.h"

class	o_Texel_f;

enum TEXTUREMODEL { MODULATE, DECAL, BLEND };

/*!
  \class o_TextureMap o_TextureMap.h
  \brief class for affine texture mapping
  */


class	o_TextureMap : public o_Property    {
	public:
#define DEFAULT_MIN_FILTER_MAP           NO_SPECIFICATION
#define DEFAULT_MAG_FILTER_MAP           NO_SPECIFICATION
		o_TextureMap();

		//! creates texture map as a link to o_CPicture
		o_TextureMap(o_CPicture<unsigned char> &pic);	

		/*! creates a new map with the given dimensions. In case
		  of monocrom images bands should be 1, and 3 for
		  RGB-colorspace */
		o_TextureMap( int dx, int dy, int bands );

		~o_TextureMap();

		const   char    *GetClassName() const ;

		/*! returns an infostring describing the object-values.
		  Possible modes are 0: long infostring, 1: short infostring */
		const   char    *GetInfostring(int mode);

		void	Resize(int dx,int dy, int bands);

		/*! returns a pixel from the map. Filtering method is defined
		  by SetTexFilter in o_AffineMap and dependent on affine
		  parameters. The coordinates x and y are normalized to the
		  interval [0,1]. */
		o_Texel_f	GetTextPixel( double x, double y, const o_AffineMap &affinemap);

		//! returns the o_CPicture
		o_CPicture<unsigned char> &GetPicture() const  {
		  return *ptex_cpic;
		}
		
		//! returns the x-dimension of texture map
		int GetNx() const {
		  return ptex_cpic->GetNx();
		}
		//! returns the x-dimension of texture map
 		int GetNy() const  {
		  return ptex_cpic->GetNy();
		}

		//! returns number of bands
		int GetNoOfBands() const  {
		 return ptex_cpic->GetNoOfBands();
		}

		//! returns value of pixel
		unsigned char GetPixel(int x, int y, int band_no = 0) const  {
		  return ptex_cpic->GetPixel( x, y, band_no);
		}

		//! returns bilinear interpolated value of pixel
		double GetPixel_bilin(double x, double y,
						 int band_no = 0) const  {
		  return ptex_cpic->GetPixel_bilin( x, y, band_no);
		}

		//! sets value of pixel
		void SetPixel(unsigned char value,int x,int y,int band_no = 0) {
		  ptex_cpic->SetPixel( value, x, y, band_no);
		}

		//! sets texture map as a link to o_CPicture
		void SetTextureMap(o_CPicture<unsigned char> &cpic)  {
		  ptex_cpic = &cpic;
		}

		//! fill band with value
		void FillBand(int band_no, unsigned char value)  {
		  ptex_cpic->FillBand(band_no, value);
		}

		//! copy texture map
		o_TextureMap &Copy(const o_TextureMap &map);

		virtual void    Copy(o_Property &org);
		virtual o_Property *MakeInstance();

		/*! sets minification filter method for initializing 
		  process  */
                void SetMinFilter( int filter_method );

		/*! returns minification filter method of (MIN_POINT,
		  MIN_BILINEAR, MIN_AVERAGE)  */   
                int  GetMinFilter(); 

		/*! sets magnification filter method for initializing 
		  process */
                void SetMagFilter( int filter_method );

		/*! returns magnification filter method of (MAG_POINT,
		  MAG_BILINEAR )  */
                int  GetMagFilter(); 
		
		/*! set texture model to one of the following:
		(MODULATE, DECAL, BLEND).  */
        	void  SetTextureModel(TEXTUREMODEL m);

		/*! returns one of: (MODULATE, DECAL, BLEND).  */
        	int  GetTextureModel();

  	protected:
		o_CPicture<unsigned char> *ptex_cpic;
	private:
                o_Texel_f       point(double x, double y);
                o_Texel_f       bilinear(double x, double y); 
                o_Texel_f       average(double x, double y,const o_AffineMap & affinemap);
                o_Texel_f       average_x_bilin_y(double x, double y,const o_AffineMap & affinemap);
                o_Texel_f       average_y_bilin_x(double x, double y, const o_AffineMap & affinemap);
                o_Texel_f       average_xy(double x, double y, const o_AffineMap & affinemap);  
		int texmap_is_linked;  // flag, texture map is linked
					// to a 'o_CPicture'
                int minfilter;
                int magfilter;
		TEXTUREMODEL mod;
};

#endif
