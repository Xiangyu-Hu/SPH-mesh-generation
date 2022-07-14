#ifndef GEOMETRY_CALCULATOR_H
#define GEOMETRY_CALCULATOR_H

#include <array>
#include <vector>
#include "glbcls.h"
#include "polygonise.h"
#include "AABB_calculator.h"

using namespace std;
using namespace boost::mpi;
using namespace boost;

class GeometryCalculator {


public:
  GeometryCalculator(){};

  void Get_corner           (Real phii0[1+2*DIM_X][1+2*DIM_Y][1+2*DIM_Z], Real corner[1+2*DIM_X][1+2*DIM_Y][1+2*DIM_Z]);
  Real Get_cutcell_volume   (Real var[1+2*DIM_X][1+2*DIM_Y][1+2*DIM_Z], Real dx, Real dy, Real dz, Real dl);
  Real get_area_MQ          (Real phii[3][3][3], Real dx, Real dy, Real dz);
  Real get_area_tri         (Polygonise_MQ::XYZ p1, Polygonise_MQ::XYZ p2, Polygonise_MQ::XYZ p3);
  Real Get_surface_area     (Real var[1+2*DIM_X][1+2*DIM_Y][1+2*DIM_Z], Real dx, Real dy, Real dz, Real dl);
  int  mark_face            (Real cor[3][3]);    // mark if this is a full face
  int  check_face           (Real cor[2][2]);   // mark if this is a full face with a 4-point stencil
  int  mark_volume          (Real cor[3][3][3]); // mark if this is a full cell with a 27-point stencil
  int  check_volume         (Real cor[2][2][2]);// mark if this is a full cell with a 8-point stencil
  Real get_face             (int level_phi, Real multi, Real cor[2][2]);     //Get area of a sub-face with adaptive refinement
  Real get_subvolume        (int level_phi, Real multi, Real corner[2][2][2]);//Get volume of a sub-cell with adaptive refinement
  Real get_area             (Real phii[3][3]);     //Get volume fraction of a 2D cell from level set on cell corners
  Real get_volume           (Real phii[3][3][3]);  //Get volume fraction of a 3D cell from level set on cell corners
  Real Get_segment_length   (my_real min, my_real max, my_real v0, my_real v1);
};


#endif
