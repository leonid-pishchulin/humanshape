#ifndef	_incl_texmap
#define	_incl_texmap
//
//      o_TextureMapper.h  - Header for oGeM texmap
//
//      Oliver Grau, 28-OCT-1993
//      Ralf Toenjes 11.10.94         pixelweise 3D Projektion statt affin

#include	"o_TriangleClip3D.h"
#include	"o_TextureMap.h"
#include	"o_TextureBinding.h"
#include	"o_Texel.h"
#include	"o_AffineMap.h"


/*!
  \class o_TextureMapper o_TextureMapper.h
  \brief class for texture mapping (only polygons)
  */

class	o_TextureMapper {
	public:
		o_TextureMap 	*map;
		o_TextureBinding	*bind;
		o_TriangleClip3D 	*walk;

		/*! returns texture value for the position x,y in the
		  image plane */
		o_Texel_f	GetPixel( int x, int y );

		/*! returns texture value for the position x,y in the
		  image plane by using perspective projection */
		o_Texel_f       GetPixel_3D( double h, double v,double zn, o_ColorCamera &cam);

		//! initialises internals
                void	Init(o_ColorCamera &cam, o_Triangle &pp, int Projection_3D_flag = 0);
		o_AffineMap affinemap;
        private:
                double factor_hx, factor_vx, factor_hy, factor_vy, konst_x, konst_y, konst_ax, konst_ay;
                double offset_x, offset_y;           
};

#endif




