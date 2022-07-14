#include "AABB_calculator.h"
#include <iostream>

namespace AABB_calculator{

  bool Slab_segment(int d, const double min[3], const double max[3], const double v0[3], const double v1[3], double& f0, double& f1)
  {
    double f_dim_0, f_dim_1;

    f_dim_0 = (min[d] - v0[d])/(v1[d] - v0[d]);
    f_dim_1 = (max[d] - v0[d])/(v1[d] - v0[d]);

    if (f_dim_1 < f_dim_0){
      double tmp = f_dim_0;
      f_dim_0 = f_dim_1;
      f_dim_1 = tmp;
    }

    if (f_dim_1 < f0)
      return false;

    if (f_dim_0 > f1)
      return false;

    f0 = f_dim_0 >= f0 ? f_dim_0 : f0;
    f1 = f_dim_1 <= f1 ? f_dim_1 : f1;

    if (f0 > f1)
      return false;

    return true;
  }

  bool Segment_AABB_Intersection(const double min[3], const double max[3], const double v0[3], const double v1[3], double Pstart[3], double Pend[3])
  {
    double f0 = 0;
    double f1 = 1;

    if (!Slab_segment(0, min, max, v0, v1, f0, f1))
      return false;

    if (!Slab_segment(1, min, max, v0, v1, f0, f1))
      return false;

    if (!Slab_segment(2, min, max, v0, v1, f0, f1))
      return false;

    for (int i = 0; i < 3; ++i)
    {
      // std::cout<<"f0: "<<f0<<" | f1: "<<f1<<std::endl;
      Pstart[i] = (v1[i] - v0[i])*f0 + v0[i];
      Pend  [i] = (v1[i] - v0[i])*f1 + v0[i];
    }

    return true;
  }

}