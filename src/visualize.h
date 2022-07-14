#ifndef _VISUAL_
#define _VISUAL_

#include "glbcls.h"
#include "Mypool.h"

using namespace std;
using namespace tbb;

class Visualize{

public:
  // variables


  // functions
  Visualize(){};
  void Set_scene(char *filename, int flag);
  void Draw_polygon(char *filename,vector<int> &f_vert,vector<double> &v,int j, Real val, Real min, Real max, int flag);
  void Draw_particle(char *filename, my_real pos, int id, Real val, Real min, Real max, Real tr);
  void Draw_boundary(char *filename, Real x_min, Real x_max, Real y_min, Real y_max, Real z_min, Real z_max);
};

#endif
