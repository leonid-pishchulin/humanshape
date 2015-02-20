//
//
//	o_PointTable.h	header for class o_PointTable
//
//	Manfred Ribitzki	16-JUN-1995
//

#ifndef o_PointTable_incl
#define o_PointTable_incl

#include "o_Point.h"
#include "o_Vector.h"
#include "o_Triangle.h"
#include "o_Surface.h"
#include "o_World.h"
#include "o_HashTable.h"

/*!
  \class o_PointTable o_PointTable.h
  \brief A point container.

  You can retrieve points by its index ore vice versa.
  Note that the references are stored. All points in the table will get a link
  to the table. On deletion of the table a request for deleting is sent to the
  points.
  */

class o_PointTable : public o_MetaClass
{
private:
  //do not copy
  //! Constructs a point table from the points of Surf and its subsurfaces.
  o_PointTable (const o_PointTable &);		//there are no bodys

  o_PointTable &operator= (const o_PointTable &);

protected:
  o_HashTableSmpl<o_Point *, int> m_PointHash;
  o_Point **m_pPointVec;
  o_Vector3D *m_pOldPoints;	//the coordinates of the points at construction
                                //time or since the last call to SavePoints

  void GetTrianglePoints (o_Triangle &);
  void GetSurfacePoints (o_Surface &);
  void GetBodyPoints (o_Body &);
  void GetWorldPoints (o_World &);
  void Init();

public:
  //! Constructs a point table from a surfaces points.
  o_PointTable (o_Surface &);
  //! Constructs a point table from Bodys points.
  o_PointTable (o_Body &);
  //! Constructs a point table from the worlds points.
  o_PointTable (o_World &);

  ~o_PointTable ();

  const char *GetClassName() const;

  //! Returns reference to point, or NULL if out of bounds.
  o_Point *GetPoint (unsigned int nIndex);

  /*! Returns the index of the given point, or -1 if point does not 
    exist in the table. */
  int GetIndex (o_Point &);

  //! Gets the center of all points.
  void Center(o_Vector3D &Center);

  //! Transforms all points to new postion.
  void ApplyTransform (o_HMatrix &Mat);

  //! Returns the number of points in the table.
  unsigned long PointCount ();

  //! Saves all points in the table, so that they can be restored later.
  void SavePoints();

  //! Restores points from latest save.
  void RestorePoints ();
};


#endif
