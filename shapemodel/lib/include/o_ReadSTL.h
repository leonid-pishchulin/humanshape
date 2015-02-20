//
// o_ReadSTL.h
//
// Copyright (C) 2006 University of Adelaide
//
// Permission is granted to any individual or institution to use, copy, modify,
// and distribute this software, provided that this complete copyright and
// permission notice is maintained, intact, in all copies and supporting
// documentation.
//
// University of Hannover provides this software "as is" without
// express or implied warranty.
//
//
// History:
//
// Created: T. Thormaehlen  April 2006      Initial design and implementation.
// Updated:


#ifndef _o_ReadSTL_h
#define _o_ReadSTL_h

#include	"fileio.h"
#include	"o_World.h"

/*!
  \class o_ReadSTL o_ReadSTL.h
  \brief reads STL-file (Stereolithography) 

    The .STL File Format
    The .STL (stereolithography) file is the de-facto standard CAD representation for RP. It was established by 3D Systems in the late 80s. The .STL format of a CAD model is a faceted surface representation,i.e. a list of the triangular surfaces with no adjacency information. This is the standard input for most RP systems. There are two format for .STL files: binary and ASCII which differs as follow:<br><br>

    Binary .STL file <br>
    The binary .STL files format consists of an 80 bytes header used to describe the solid contained within the file, then 4 bytes represent the total number of facets in the  file. A facet is described as follow: the first 12 bytes (3 x 4 bytes) represent its normal, the next 36 bytes (3 x 3 x 4 bytes) represent its (three) vertices, then two unused bytes are padded to achieve a block size of 50 bytes.

    ASCII .STL files <br>
    ASCII files use keywords and are self explanatory. The ASCII .stl file must start with the lower case keyword solid and end with endsolid. Within these keywords are listings of individual triangles that define the faces of the solid model. Each individual triangle description defines a single normal vector directed away from the solid's surface followed by the xyz components for all three of the vertices. These values are all in Cartesian coordinates and are floating point values. The triangle values should all be positive and contained within the building volume. The normal vector is a unit vector of length one based at the origin. If the normals are not included then most software will generate them using the right hand rule. If the normal information is not included then the three values should be set to 0.0. There is a variety of errors in ASCII files that do not appear in binary files. For instance, it happens that keywords are either skipped of extraneous, hindering the extraction of data. Here's an example of an .STL ASCII file:<br>
    solid Solidname<br>
    facet normal 9.838605e-01 3.226734e-02 1.760037e-01<br>
      outer loop<br>
	vertex   -1.070000e+02 0.000000e+00 1.816000e+02<br>
	vertex   -1.060000e+02 0.000000e+00 1.760100e+02<br>
	vertex   -1.070000e+02 1.200000e+00 1.813800e+02<br>
      endloop<br>
    endfacet<br>
    facet normal 9.824255e-01 9.205564e-02 1.623759e-01<br>
      outer loop<br>
      vertex   -1.070000e+02 1.200000e+00 1.813800e+02<br>
      vertex   -1.060000e+02 0.000000e+00 1.760100e+02<br>
      [...]<br>
      endloop<br>
    endfacet<br>
    [...]<br>
    endsolid<br>
  */

class o_ReadSTL
{
  public:
    o_ReadSTL();
    ~o_ReadSTL();

  public:
    bool ReadFile(o_World & world, const char* filename);
};

#endif	

