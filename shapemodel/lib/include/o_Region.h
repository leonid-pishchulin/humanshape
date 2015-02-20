#ifndef	__inc_o_Region
#define	__inc_o_Region
//
//      Oliver Grau, 12-APR-1996
//
//	$Id: o_Region.h 433 2007-09-27 06:10:41Z thormae $


class	o_LabelImage;

/*!
  \class o_Region o_Region.h
  \brief pixel region 

  The objects of this class are build up from the class o_LabelImage .
  Each region is characterized by a unique label number (or value) in
  the corresponding label image. For stepping over all pixel positions
  an iterator can be used \sa o_RegionIterator .
  */

class o_Region
{
	public:
    //! is regularly instantiated by o_LabelImage::InitRegions()
		o_Region();

		virtual ~o_Region();
		virtual const   char *GetClassName() const;

		/*! sets the coordinates of the bounding box in the
		  corresponding	label image to x1,y1 ; x2,2 */
		int	BoundingBox( int &x1, int &y1, int &x2, int &y2 );

		//! returns the number of pixels of the region
		int	GetNumberOfPixels();

		/*! return the label value. That is the pixel value of the
		  region in the label image */
		long	GetLabelValue();

		//! returns the reference to the corresponding label image
		o_LabelImage	&GetLabel();

		/*! returns 1 if r is a neighboring region. If *x and *y 
		are given they contain the first location of a transmission
		of the regions in the label image are returned. If
		neighborhood8 is true a 8-neighborhood is considered */
		int	IsNeighbor(o_Region &r, int neighborhood8=0, int *x=0, int *y=0 );

	protected:
		friend class o_LabelImage;

		void	Init(o_LabelImage &lab, long labelval, int firstx, int firsty);
		void	CheckBoundingBox(int x, int y) {
			++pixels;
			if(x<xmin) xmin=x;
			if(y<ymin) ymin=y;
			if(x>xmax) xmax=x;
			if(y>ymax) ymax=y;
		};

		long	val;
		o_LabelImage	*label;
		int	xs,ys;
		int	xmin,ymin;
		int	xmax,ymax;
		int	pixels;
};

#endif

