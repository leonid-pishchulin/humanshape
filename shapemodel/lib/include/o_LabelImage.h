#ifndef	__inc_o_LabelImage
#define	__inc_o_LabelImage
//
//	Oliver Grau, 12-APR-1996
//      Jochen Wingbermhle 1997
//      Sebastian Weik 28-Nov-1997
//
//	$Id: o_LabelImage.h 474 2007-10-11 07:01:34Z rhys $



#include "o_Region.h"
#include "o_Picture.h"
#include "o_CPicture.h"
#include "o_HashTable.h"
#include "o_Vector.h"
#include <list>
using std::list;

/*!
  \class o_LabelImage o_LabelImage.h
  \brief Is a (parent class of) o_CPicture<long>. The pixel values are
	interpreted as label numbers (or values)

  The class handles different means of creating label images from binary 
  images. After The label image has been created, several methods exist to
  examine the labelel regions.

  1. The method InitRegions() builds up an internal list with objects of from 
  the class o_Region. Currently there is no incrementel update possible! If 
  a new label is added to the label image or some regions are changed a call 
  to InitRegions() must be done. On destruction of a o_LabelImage object all
  associated region objects are deleted !
  
  2. The method FillNeighborList() fills an array of std::list<long*> with a list of neighbors labels
  for each region. Usage Recommended !

*/

class o_LabelImage : public o_CPicture<long>
{
 public:
  /*!
    creates a label image with the specified dimension (dx,dy)
  */
  o_LabelImage(int dx=0, int dy=0);
	
  /*!
    Destructor. Deletes all associated regions
  */
  ~o_LabelImage();

  /*!
    returns name of class
  */
  virtual const   char *GetClassName() const;

  /*!
    Scans all pixels and creates an associated object of the class
    o_Region for each label value. Before building up a new list all 
    region objects are deleted!   
  */
  void InitRegions();

  /*!
    returns the number of regions in the label image.
    InitRegions() must be called before
  */
  int	RegionListLen();

  /*!
    returns pointer to the n'th region object.
    InitRegions() must be called before
  */
  o_Region *GetRegion(int n);

  /*!
    returns pointer to the region object with lab as label value.
    InitRegions() must be called before
  */
  o_Region *GetRegionByValue(long lab);

  o_LabelImage & operator = (const o_LabelImage & img);
  o_LabelImage & operator = (const o_CPicture<long> & img);

  /*!
    Fills the Labelimage with labels according to connected regions (4-Point Neighborhood) 
    in the single band input image (pic). All points with a value below thres1 
    or above thres2 are considered as belonging to the background (label 0).
  */
  long CreateLabelImage4(o_CPicture<unsigned char> & pic,
			 unsigned char thresh1,unsigned char thresh2=255);
  /*!
    Fills the Labelimage with labels according to connected regions (8-Point Neighborhood) 
    in the single band input image (pic). All points with a value below thres1 
    or above thres2 are considered as belonging to the background (label 0).
  */
  long CreateLabelImage8(o_CPicture<unsigned char> & pic,
			 unsigned char thresh1,unsigned char thresh2=255);

  /*!
    Fills the Labelimage with labels according to connected regions (4-Point Neighborhood)
    in the single band input image (pic). All points with a value below thres1 
    or above thres2 are considered as belonging to the background (label 0).
    Works line recursively, thus also large images can be handled
  */
  long CreateLabelImage(o_Picture<unsigned char> & pic,
			unsigned char thresh1,unsigned char thresh2=255);

  /*!
    H� ? 0 wenn pic gleich gross wie labelimage , sonst 1
  */
  int  Add_Label(long number,int value,o_CPicture<unsigned char> &pic);

  /*!
    Clears region list
  */
  void ClearRegionList();

  /*!
    Returns the number of labels that have been previously created by CreateLabelImage4() or 
    CreateLabelImage8(). Does not return valid value otherwise.
  */
  long GetMaxLabel(){return number_of_labels;}

  /*!
    Struct describing the statistics of a single label
    sumx - sum of x coordinates of all pixels of label
    sumy - sum of y coordinates of all pixels of label
    npix - number of pixels of label
    m - 2D Vector describing the mean postion of label
  */
  typedef struct LabelStatistics {
    long sumx;
    long sumy;
    long npix;
    o_Vector2D m;
  } LabelStatistics;
		
  /*! Fills the member array that is an array of std::list<long*>. 
     Each element in the array describes
    the neighbors of that particular region. E.g. arr[4]->Size() gives the number of neighbored
    regions for the region with label number four.
    Returns the neighbor array
   */    
  std::list<long*> ** FillNeighborList();

  /*!
    Fills not only the array conatining the the neighbors but also an array with a structure
    LabelStatistics as described above. E.g. stc_arr[4]->m.GetX() yields the averaged x-coordinate
    of the region with label number 4
    Returns the neighbor array
  */
  std::list<long*> ** FillNeighborAndStatisticsList();

  /*!
    Returns the statistics array, see above. Has to be filled prior to quering
  */
  LabelStatistics ** GetStatisticsList(){return stc_arr;}
  
  /*!
    Searches the neighbor array for a region with exactly n neighbors. The label of the region is
    returned. array_length is the total number of regions
  */
  long GetRegionWithNNeighbors(long arr_length, long n);

		
 protected:

  friend class o_Region;

  o_HashTableSmpl<long,o_Region>	htabl;

  int	index;


  //  for region labeling J.W. 9/97;
		   
  /*!
    recursive labeling function for 4 neighborhood
  */
  void label4 (void);

  /*!
    recursive labeling function for 8 neighborhood
  */
  void label8 (void);

  /*!
    recursive labeling function for 4 neighborhood filling a horizontal line and then going recursive
  */
  void label4 (short x);

  int xpos,ypos;
  long max_region_size;
  long current_label;
  unsigned char current_value;
  long current_region_size;

  o_CPicture<unsigned char> *input_single_band_image;
  unsigned char lower_thres;
  unsigned char upper_thres;
  long biggest_region_label;

  long number_of_labels;
  //contains the number of labels od last labeling call

  /*!
    Inserts a pair of neighboring labels into neighborlistarray arr
    If the pair of neighbors already exists in list operation is skipped
  */
  void InsertNeighbor(long & lbl, long & nghb);

  std::list<long*> ** arr;
  //pointer to an array of lists containing the neighbors for a given region

  LabelStatistics ** stc_arr;
  //pointer to an array of Statistic structs containing some information on the regions

};



#endif



