#ifndef	_o_metacamera
# define	_o_metacamera

//
//	o_MetaCamera.h	- Camera class definitions
//
//	Oliver Grau, Jan./Feb.1993
//	Reinhard Koch Mar. 94
//	Oliver Grau, Feb.1998, revised docu and SetFocalLength


//#include "outputdev.h"
#include "o_World.h"

#ifndef TRUE
# define TRUE 1
#endif

#ifndef FALSE
# define FALSE 0
#endif

//
//	constants
//
const	double	z_infinite=(+1e38);


/*!
  \class o_MetaCamera o_MetaCamera.h
  \brief abstract camera class - do not instantiate
  */

class	o_MetaCamera : public o_LocalCoordinateSystem {
	public:
#define	CAM_DEFAULT_WID	     256
#define	CAM_DEFAULT_HEI	     256


    // rib:
    //
    // rd_dim is a hack for keeping compatibility to the previous
    // implementation of radial distortion stuff. The previous
    // implementation was based on k3 and k5 given in pel^-2 and
    // pel^-4 respectively (wrong in my opinion).
    //
    // Different dimensions for k3 and k5 (mm <--> pel) are not just a
    // simple scaling problem. This would be the case only if the
    // pixel size is quadratic. Thus we need the information about
    // the dimensions of k3, k5.
    // The base dimension for the radial distortion params is stored
    // in the private member rd_dim. v(k3) is given in rd_dim^-2 and
    // v5(k5) is given in rd_dim^-4 then.
    //
    // The ***new*** default for the base dimension is mm, because rd
    // params based on mm are more common (and meaningful!).
    //

    enum RD_DIM {
	e_rd_mm,	//k3 in mm^-2, k5 in mm^-4
	e_rd_pel	//k3 in pel^-2, k5 in pel^-4
    };


               	virtual const   char    *GetClassName() const;

		//! returns the C - vector in [mm], R3
		virtual const   o_Vector3D        &Origin() const;
		//! returns the H0,V0,A vectors for n=0,1,2
		virtual const   o_Vector3D        &Axis(int n) const;

		//! returns the camera parameter c-vector in [mm], R3
		const	o_Vector3D	&C() const { return c_vec;}
		//! returns the camera parameter a-vector (unit vector)
		const	o_Vector3D	&A() const { return a_vec;}
		//! returns the camera parameter h-vector in [pel], R3
		const	o_Vector3D	&H() const { return h_vec;}
		//! returns the camera parameter v-vector in [pel], R3
		const	o_Vector3D	&V() const { return v_vec;}

		//! returns the width (pixel) of the camera plane 
		virtual int	Width() const { return dx; }
		//! returns the height (pixel) of the camera plane 
		virtual int	Height() const { return dy; }

		/*! delivers 1 (TRUE), when coordinates (x,y) in <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> are inside camera plane and 0 (FALSE), when outside */
		int     InPicture(double x, double y) const
		{
		  // test if coordinates x,y are inside of picture
		  // test performed in pcs

		  if (x < 0 || x > dx-1 || y < 0 || y > dy-1)
		    return FALSE;
		  else
		    return TRUE;
		}

		/*! converts coordinates from <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> to <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#ics">ICS</A> coordinates */ 
		void    pcs2ics (double px, double py, 
				 double &ix, double &iy) const
		{
		  // convert pcs to ics coordinates
		    ix = px - 0.5*(Width() -1);
		    iy = py - 0.5*(Height() -1);
 		}
		
		/*! converts coordinates from <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#ics">ICS</A> to <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> coordinates.
		<A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> is the preferred picture coordinate system with origin in the upper left corner, <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#ics">ICS</A> is centered in the middle of the picture */
		void    ics2pcs (double ix, double iy, 
				 double &px, double &py) const
		{
		  // convert ics to pcs coordinates
		    px = ix + 0.5*(Width() -1);
		    py = iy + 0.5*(Height() -1);
 		}
	
		//! returns size of a pixel on the image plane in [mm/pel]
		double	PixelSizeX() const { return sx; }
		//! returns size of a pixel on the image plane in [mm/pel]
		double	PixelSizeY() const { return sy; }

		double	RadialDistortion() const { return v; }

		//! returns the radial distortion, term K3 in [RdDim()-2]
		double  RadialDistortion_K3() const { return v; }
		//! returns the radial distortion, term K3 in [RdDim()-4]
		double  RadialDistortion_K5() const { return v5; }
                //! Returns the base dimension of radial distortion params.
                // This could be one o_MetaCamera::e_rd_mm or 
                // o_MetaCamera::e_rd_pel
    o_MetaCamera::RD_DIM RdDim() const { return rd_dim; }

		//! returns the focal length in [mm]
                double  FocalLength() const ;

		//! returns the image plane vector in x-direction (unit vector)
		const   o_Vector3D &H0() const ;
		//! returns the image plane vector in y-direction (unit vector)
		const   o_Vector3D &V0() const ;

		//! returns the centerpoint of the image plane in [pel]
		double  CenterPointX() const ;
		//! returns the centerpoint of the image plane in [pel]
		double  CenterPointY() const ;

		/*! returns the line of sight as a unit vector. 'x' and 'y'
                are given in <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> */
		const   o_Vector3D &LineOfSight(double x, double y) const;

		//! returns homogeneous matrix of camera orientation
		const   o_HMatrix &Transformation() const;

		void	Project( const o_Vector3D &, double &, double & ) const ;
		void	Project( const o_Vector3D &, double &, double &, double & ) const ;
		/*! calculates the position of the projection of 'point' 
		in picture coordinates <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> of the camera . The results are stored to 'x' and 'y'. */
		virtual void	Project_PCS( const o_Vector3D &, double &, double & ) const ;

		/*! calculates the position of the projection of 'point' 
		in picture coordinates <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A> of the camera . The 
		results are stored to 'x' and 'y'.*/
		virtual void	Project_PCS( const o_Vector3D &, double &, double &, double & ) const ; 

                /*! calculates the position of the projection of 'point' in
		picture coordinates <A HREF="http://www.tnt.uni-hannover.de/project/3dmod/camera/central/cahv.html#pcs">PCS</A>
		of the camera . The results are stored to 'x' and 'y'. 'z'
		contains the depth of the 3D point. Note that when 'z' is
		less than the focal length (or negative), the 3D point is
		located behind the image plane. success will be set to
		false only if there is no numerical solution (e.g. success
		will be true for negative values of 'z'). When
		check_image_bounds is set to false, it is the callers
		responsibility to check whether or not the projected 2D
		point lies within image boundaries. When
		check_image_bounds is set to true, Project_PCS will fail
		when the projected point is outside the image boundaries*/
                void Project_PCS (const o_Vector3D &vec, double &x, double &y, double &z, bool &success, bool check_image_bounds = true) const ;


		//! apply transformation matrix 'tmat' to orientation
                void    ApplyTransform( const o_HMatrix &tmat );

		//! set parameter like in cam
		virtual void	Set_Camera(const o_MetaCamera &cam );
		
		/*! set camera with parameter set. ci in [mm], R3, ai (unit vector),
		hi in [pel], R3, vi in [pel], R3, width in [pel], height in [pel],
		rdist_k3 in [pel-2], rdist_k5 in [pel-4], 
		isx in [mm/pel], isy in [mm/pel] */
		virtual void	Set_Camera( const o_Vector3D &ci,
			const o_Vector3D &ai, const o_Vector3D &hi,
			const o_Vector3D &vi, int	width=CAM_DEFAULT_WID,
			int  height=CAM_DEFAULT_HEI, double rdist_k3=0.0,
			double isx=1, double isy=1, double rdist_K5=0.0,
			o_MetaCamera::RD_DIM dim = o_MetaCamera::e_rd_mm);

		virtual	void	Clear();

		//! set camera parameter with tsai parameter set
		void Set_Camera( o_Vector3D t, double tilt, double roll,
			double pan, 
			double xshift, double yshift,
			int	width, int height,
			double focal_len,
			double xscale =1.0,
			double isx=1, double isy=1,
			double	k3=0.0, double  k5=0.0,
			o_MetaCamera::RD_DIM dim = o_MetaCamera::e_rd_mm
			);

    void SetRadialDistortion (double k3, double k5, 
			      o_MetaCamera::RD_DIM dim = o_MetaCamera::e_rd_mm);

		//! Set the focal length in [mm]
                void  SetFocalLength(double nf) ;

		//! Set the x-center point (principal point) in [pel]
                void  SetCenterPointX(double cpx) ;

		//! Set the y-center point (principal point) in [pel]
                void  SetCenterPointY(double cpy) ;

		/*! writes out the camera parameters to wfp. If bin_flag is 0 the ascii-fileformat is taken, otherwise the binary format. If version is 10 the CAHV Verion 1.0 is written, otherwise CAHV Version 1.1 \sa <A HREF="../../cam.html">VirtualCamera-Fileformat</A> */
		void  Describe( write_file &wfp, int bin_flag=0, int version = 11) const ;
	
		/*! read camera paramters from file rfp. If no header is given
		at the first line, the method interpretes the file as 
		tsai-parameter file. If the file is in tsai file format
		the image dimensions must be known before. I.e. the image
		should be read in before (in case of real cameras).
		 \sa <A HREF="../../cam.html">VirtualCamera-Fileformat</A> */
                virtual int Readparameter( read_file &rfp );
                int   ReadTSAIparameter( read_file &rfp );
		
		/*! alternative method to read camera paramters form 
		  file name. If no header is given
		  at the first line, the method interpretes the file as 
		  tsai-parameter file. If the file is in tsai file format
		  the image dimensions must be known before. I.e. the image
		  should be read in before (in case of real cameras).
		  \sa <A HREF="../../cam.html">VirtualCamera-Fileformat</A> */
		virtual int Readparameter(const char *name);
		int   ReadTSAIparameter(const char *name);

	protected:
                 o_MetaCamera( );
		~o_MetaCamera() ;
	private:
		int	dx;	// Window (Picture) width
		int	dy;	// Window (Picture) height
		double	v;	// radial distortion
                o_Vector3D	c_vec,a_vec,h_vec,v_vec;
                double	sx,sy;	// Pixel size in the image plane
		double	v5;	// radial distortion
                RD_DIM rd_dim;
};


#endif






