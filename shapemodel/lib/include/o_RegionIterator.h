#ifndef __inc_o_RegionIterator
#define __inc_o_RegionIterator
//
//      Oliver Grau, 12-APR-1996
//
//      $Id: o_RegionIterator.h 433 2007-09-27 06:10:41Z thormae $

class	o_Region;
class	o_LabelImage;

/*!
  \class o_RegionIterator o_RegionIterator.h 
  \brief Iterator for pixel wise access to regions

  iterator for pixel wise access to image regions o_Region
  */

class o_RegionIterator {
	public:
    //! initalizes iterator for region reg
		o_RegionIterator(o_Region &reg);

		~o_RegionIterator();

		//! init iterator
		void	Init();

		/*! returns 0 if no more points are in the region or returns
		  1 and sets x,y to the	coordinates of the next point
		  belonging to the region. */
		int	Next(int &x, int &y);

	protected:
		o_Region *reg_ptr;
		int	x,y;
		int	x1,x2,y1,y2;
		o_LabelImage	*li;
};

#endif
