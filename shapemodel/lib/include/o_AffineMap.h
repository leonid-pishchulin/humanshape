#ifndef	_o_affinemap_incl
#  define	_o_affinemap_incl

//
//	o_AffineMap.h	- Header for oGeM affine mapping
//
//	Hellward Broszio, Apr. 1994
//      Ralf Toenjes 11.10.94       Min-/MagTextureFilter
//

#include	"oGeM.h"
#include	<vector>
using std::vector;

/*! 
  \class o_AffineMap o_AffineMap.h
  \brief class for affine mapping
  */

class	o_AffineMap : public o_MetaClass {
	public:
                /*! creates affine mapping, the parameters of affine
		  mapping are set to zero values. */
		o_AffineMap();

		~o_AffineMap();
                virtual const   char    *GetClassName() const;

		//! calculates for a target vector(h,v) the source vector(x,y).
		void Mapping(o_Vector2D &vec_xy, const o_Vector2D &vec_hv) const {

		  vec_xy.SetX( a[0]*vec_hv.GetX() +a[1]*vec_hv.GetY() +a[2] );
		  vec_xy.SetY( a[3]*vec_hv.GetX() +a[4]*vec_hv.GetY() +a[5] );

		}

		/*! calculates for the target coordinates h and v the source
		coordinates x and y. */	
		void Mapping(double *x, double *y, const double h, const double v) const {
		  *x =  a[0]*h + a[1]*v + a[2]; 
		  *y =  a[3]*h + a[4]*v + a[5]; 
		} 

		/*! sets affine parameters in dependence of the three
		  2D-vectors o_Vector2D in both lists. The mapping is
		  from the coordinates system (x,y) to the 
		  coordinate system (h,v). */
    void SetParameter(std::vector<o_Vector2D*> &vector_list_xy, 
                      std::vector<o_Vector2D*> &vector_list_hv);

		//! returns affine parameter Index (0 <= index <= 5).
		double GetParameter(int index) const;

		/*! Sets texture filter to minification or magnification out
		  of ( MIN_POINT, MAG_POINT, MIN_BILINEAR, MAG_BILINEAR,
		  MIN_AVERAGE ), \sa \link texture.html texture \endlink */
                void SetTexFilter(int method){ texture_filter_method = method;}

		/*!  Gets applied Min/MagFilter, \sa \link texture.html texture
		  \endlink */
                int GetTexFilter() const {return texture_filter_method;}

		/*! Calculates max. absolute affine parameters by using
		  viewplane-height and viewplane-width from given camera */
                void SetScaleParameter(int hei=0, int wid=0);

		//! Gets max. absolute value of affine parameter a0, a1
                int GetScale_x() const { return scale_x; }
		//! Gets max. absolute value of affine parameter a3, a4
                int GetScale_y() const { return scale_y; }

	protected:
		void SetParameter(const double *x, const double *y, 
				  const double *h, const double *v);
	private:
                int scale_x, scale_y;
                int texture_filter_method;  
		double a[6];

};


#endif
