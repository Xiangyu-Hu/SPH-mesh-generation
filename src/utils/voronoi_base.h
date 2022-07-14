#ifndef _VORONOI_BASE_H_
#define _VORONOI_BASE_H_
#include "glbcls.h"
#include "Mypool.h"
#include <boost/polygon/voronoi.hpp>
#include "voro++.hh"

using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;
using boost::polygon::x;
using boost::polygon::y;
using boost::polygon::low;
using boost::polygon::high;
using namespace std;
using namespace tbb;
using namespace boost;
using namespace boost::mpi;
using namespace voro;

struct Point {
  int a;
  int b;
  int c;
  Point(int x, int y, int z) : a(x), b(y), c(z) {}
};
struct Point_real {
  Real a;
  Real b;
  Real c;
  Point_real(Real x, Real y, Real z) : a(x), b(y), c(z) {}
};
struct Segment {
  Point p0;
  Point p1;
  Segment(int x1, int y1, int z1, int x2, int y2, int z2) : p0(x1, y1, z1), p1(x2, y2, z2) {}
};

namespace boost {
namespace polygon {

template <> struct geometry_concept<Point> {
  typedef point_concept type;
};

template <> struct point_traits<Point> {
  typedef int coordinate_type;

  static inline coordinate_type get(const Point& point, orientation_2d orient) {
    return (orient == HORIZONTAL) ? point.a : point.b;
  }
};

template <> struct geometry_concept<Segment> {
  typedef segment_concept type;
};

template <> struct segment_traits<Segment> {
  typedef int coordinate_type;
  typedef Point point_type;

  static inline point_type get(const Segment& segment, direction_1d dir) {
    return dir.to_int() ? segment.p1 : segment.p0;
  }
};
}  // polygon
}  // boost

class Vo_particle {

  private:
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & tag;
    if (tag == 1){ // for broadcasting VP coord
      ar & coord;
      ar & id;
    }
    else if (tag == 2){ // for gathering VP information
      ar & id;
      ar & mass;
      ar & mass_center;
      ar & energy;
    }
    else if (tag == 3){ // for gathering VP information
      ar & id;
      ar & mass;
      ar & coord;
    }
    else if (tag == 4){ // for gathering VP information
      ar & id;
      ar & h_min;
      ar & coord;
    }
  }
    
  public:
  // variables
  int      tag;
  int      id;
  int      color;
  Real     mass;
  Real     P;
  Real     rho;
  Real     vol;
  Real     h;
  Real     h_min;
  Real     timestep;
  Real     timestep_CVT;
  Real     energy;
  Real     error;
  my_real  coord;
  my_real  a;
  my_real  f;
  my_real  v;
  my_real  CVT_shift;
  my_real  mass_center;

  // functions
  Vo_particle (){};
  void     Initialize(int index);
  void     Reset();

};

class Voronoi_base {

  public:
  // variables
  voronoi_diagram<Real>             vd;

  // functions
  Voronoi_base (){};
};

#endif