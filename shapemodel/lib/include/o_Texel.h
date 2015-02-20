#ifndef	_incl_o_texel
#define	_incl_o_texel
//
//      o_Texel.h  - Header for oGeM texel
//
//      Ralf Toenjes, 18-MAR-1996
//

#include	"o_Color.h"


/*!
  \class o_Texel o_Texel.h
  \brief class for texel
  */

class	o_Texel : public o_Color {
	public:
		o_Texel() :o_Color() { };
		o_Texel(unsigned char ri, unsigned char gi, 
                        unsigned char bi, unsigned alphai = 255) 
		  :o_Color(ri, gi, bi){
		    alpha = alphai;
		};

		//! set r,g,b,alpha values
		void SetValue(unsigned char ri, unsigned char gi,
                              unsigned char bi, unsigned char alphai = 255) {
                        o_Color :: SetValue(ri, gi, bi);
                        alpha = alphai;
                };

		//! get value of the alpha band
                unsigned char GetAlpha() const { return alpha; }

	private:
		unsigned char alpha;
};


/*!
  \class o_Texel_f o_Texel.h
  \brief class for color
  */

class	o_Texel_f : public o_Color_f {
	public:
		//! creates texel object with values 0
		o_Texel_f() :o_Color_f() { };

		//! creates object and initialises values from tex
	        o_Texel_f(const o_Texel &tex) 
		  :o_Color_f(tex.GetR(), tex.GetG(), tex.GetB()) { 
		    alpha = tex.GetAlpha();
		};

		//! creates texel object with specified values
                o_Texel_f(float ri, float gi, float bi, float alphai = 1.0){ 
		    SetValue(ri,gi,bi,alphai);
		};

		//! set r,g,b,alpha values
		void SetValue(float ri, float gi, float bi, float alphai= 1.0){
		  o_Color_f :: SetValue(ri, gi, bi);
		  alpha = alphai;
		  ck_val();
                };

		//! set texel to val
		void SetValue( const o_Texel_f &val ) {
		  o_Color_f :: SetValue(val.GetR(), val.GetG(), val.GetB());
		  alpha = val.alpha;
		  ck_val();
                };

		//! assignment of a texel object
		o_Texel_f& operator=(const o_Texel_f& val) {
			SetValue( val );
			return *this;
		};

		//! get value of the alpha band
        	float GetAlpha() const { return alpha; } 

	private:
		void ck_val() {
                        if(alpha<0.0) alpha=0;
                        if(alpha>1.0) alpha=1; 
		};
		float alpha;
};

#endif

