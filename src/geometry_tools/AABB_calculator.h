#ifndef _AABB_CALCULATOR_H_
#define _AABB_CALCULATOR_H_

namespace AABB_calculator{

bool Slab_segment(int d, const double min[3], const double max[3], const double v0[3], const double v1[3], double& f0, double& f1);

bool Segment_AABB_Intersection(const double min[3], const double max[3], const double v0[3], const double v1[3], double Pstart[3], double Pend[3]);
}

#endif