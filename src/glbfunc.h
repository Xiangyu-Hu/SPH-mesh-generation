#ifndef GLBFUNC_H
#define GLBFUNC_H

#include "glbparam.h"
#include "glbcls.h"
#include<cstdlib>
#include <cmath>
#include <boost/serialization/string.hpp>

using namespace std;
using namespace tbb;
using namespace boost;

#define I_BITS 10
#define J_BITS 10
#define K_BITS 10
#define L_BITS (32 - (I_BITS + J_BITS + K_BITS))

#define I_MASK ((1 << I_BITS) - 1)
#define J_MASK ((1 << J_BITS) - 1)
#define K_MASK ((1 << K_BITS) - 1)
#define L_MASK ((1 << L_BITS) - 1)

#define I_SHIFT 0
#define J_SHIFT (I_SHIFT + I_BITS)
#define K_SHIFT (J_SHIFT + J_BITS)
#define L_SHIFT (K_SHIFT + K_BITS)

template<typename T>
typename std::enable_if<std::is_unsigned<T>::value, int>::type
inline constexpr signum(T x) {
    return T(0) < x;
}

template<typename T>
typename std::enable_if<std::is_signed<T>::value, int>::type
inline constexpr signum(T x) {
    return (T(0) < x) - (x < T(0));
}

//  the nth power
template<class T> T powern(T a, int n)
{
  T res = 1;
  if (n >= 0)
    for (int i=0; i<n; i++)
      res *= a;
  else if (n < 0){
    for (int i=0; i<(-1*n); i++)
      res *= a;
    res = 1/res;
  }
  return res;
};
//  Get the maximum
template<class T> T AMAX1(T a, T b)
{
  return (a >= b ? a : b);
};
//  Get the minimum
template<class T> T AMIN1(T a, T b)
{
  return (a <= b ? a : b);
};
//  multiply constant value to my_int or my_real
template<class T1, class T2> T1 my_multiply_const(T1 a, T2 b)
{
  T1 c;
  c.i = DIM_X==1 ? a.i*b : a.i;
  c.j = DIM_Y==1 ? a.j*b : a.j;
  c.k = DIM_Z==1 ? a.k*b : a.k;
  return c;
};
//  multiply self-defined data type to my_int or my_real
template<class T1, class T2> T1 my_multiply_data(T1 a, T2 b)
{
  T1 c;
  c.i = DIM_X==1 ? a.i*b.i : a.i;
  c.j = DIM_Y==1 ? a.j*b.j : a.j;
  c.k = DIM_Z==1 ? a.k*b.k : a.k;
  return c;
};
//  devide self-defined data type
template<class T1, class T2> T1 my_devide_data(T1 a, T2 b)
{
  T1 c;
  c.i = DIM_X==1 ? a.i/b.i : a.i;
  c.j = DIM_Y==1 ? a.j/b.j : a.j;
  c.k = DIM_Z==1 ? a.k/b.k : a.k;
  return c;
};
//  add constant value to my_int or my_real
template<class T1, class T2> T1 my_add_const(T1 a, T2 b)
{
  T1 c;
  c.i = DIM_X==1 ? a.i+b : a.i;
  c.j = DIM_Y==1 ? a.j+b : a.j;
  c.k = DIM_Z==1 ? a.k+b : a.k;
  return c;
};
//  add self-defined data type to my_int or my_real
template<class T1, class T2> T1 my_add_data(T1 a, T2 b)
{
  T1 c;
  c.i = DIM_X==1 ? a.i+b.i : a.i;
  c.j = DIM_Y==1 ? a.j+b.j : a.j;
  c.k = DIM_Z==1 ? a.k+b.k : a.k;
  return c;
};
//  minus self-defined data type
template<class T1, class T2> T1 my_minus_data(T1 a, T2 b)
{
  T1 c;
  c.i = DIM_X==1 ? a.i-b.i : a.i;
  c.j = DIM_Y==1 ? a.j-b.j : a.j;
  c.k = DIM_Z==1 ? a.k-b.k : a.k;
  return c;
};
//  set the value of my_int or my_real according to the same data type
template<class T> void my_set_data(T &a, T b)
{
  a.i = DIM_X==1 ? b.i : a.i;
  a.j = DIM_Y==1 ? b.j : a.j;
  a.k = DIM_Z==1 ? b.k : a.k;
};
//  set the value of my_int or my_real according to another constant
template<class T1, class T2> void my_set_const(T1 &a, T2 b)
{
  a.i = DIM_X==1 ? b : 1;
  a.j = DIM_Y==1 ? b : 1;
  a.k = DIM_Z==1 ? b : 1;
};
//  multiply all the members
template<class T1, class T2> void my_self_multiply(T1 a, T2 &b)
{
  b = 1;
  b *= (DIM_X==1 ? a.i : 1);
  b *= (DIM_Y==1 ? a.j : 1);
  b *= (DIM_Z==1 ? a.k : 1);
};
//  vector mold
template<class T1, class T2> void my_self_mold(T1 a, T2 &b)
{
  b = 0;
  b += (DIM_X==1 ? a.i*a.i : 0);
  b += (DIM_Y==1 ? a.j*a.j : 0);
  b += (DIM_Z==1 ? a.k*a.k : 0);
  b = sqrt(b);
};
//  add all the members
template<class T1, class T2> void my_self_add(T1 a, T2 &b)
{
  b = 0;
  b += (DIM_X==1 ? a.i : 0);
  b += (DIM_Y==1 ? a.j : 0);
  b += (DIM_Z==1 ? a.k : 0);
};
//  output
template<class T> void my_cout(T a, const char *title)
{
  cout<<"   "<<title<<"  "<<a.i<<"  "<<a.j<<"  "<<a.k<<"  "<<"\n";
};
//allocate & delete 3d martrix
template<class T> void allocate_3d_matrix(T*** &matrix, my_int res)
{
  matrix = new T**[res.i];
  for(int i=0; i<res.i; i++){
    matrix[i] = new T*[res.j];
    for(int j=0; j<res.j; j++){
      matrix[i][j] = new T[res.k];
    }
  }
};
// template <typename T> //switched order for deduction
// boost::multi_array<T, 3> allocate_3d_multi_array (int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, T value)
// {
//   boost::multi_array<T, 3> arr(boost::extents[range(xmin, xmax)][range(ymin, ymax)][range(zmin, zmax)]);
//   std::fill(arr.data(), arr.data() + arr.num_elements(), value);
//   return arr;
// };
//d. Sign of the first value is determined by the secomax_species_number value's sign
inline Real SIGN(Real a, Real b){
  if(b > 0.0) return fabs(a);
    return -fabs(a);
}
inline Real sgn(Real a)
{
  if(a > 0.0) return 1.0;
  return -1.0;
}
inline int sgn1(int a)
{
  if(a >= 0) return 1;
  return -1;
}
inline my_int get_cell_id(my_real coord, my_real dcell, my_int cell_start, my_int cell_end){
  my_int pos;

  pos.i = DIM_X==1 ? int(floor(coord.i/dcell.i)) : 0;
  pos.j = DIM_Y==1 ? int(floor(coord.j/dcell.j)) : 0;
  pos.k = DIM_Z==1 ? int(floor(coord.k/dcell.k)) : 0;

  pos.i = AMAX1(cell_start.i,AMIN1(pos.i,cell_end.i-1));
  pos.j = AMAX1(cell_start.j,AMIN1(pos.j,cell_end.j-1));
  pos.k = AMAX1(cell_start.k,AMIN1(pos.k,cell_end.k-1));

  return pos;
};

inline my_int get_subcell_id(my_real coord, my_real dcell, my_int nsubcell){
  my_int pos;

  pos.i = DIM_X==1 ? int(floor(coord.i/dcell.i)) : 0;
  pos.j = DIM_Y==1 ? int(floor(coord.j/dcell.j)) : 0;
  pos.k = DIM_Z==1 ? int(floor(coord.k/dcell.k)) : 0;

  pos.i = AMAX1(0,AMIN1(pos.i,nsubcell.i-1));
  pos.j = AMAX1(0,AMIN1(pos.j,nsubcell.j-1));
  pos.k = AMAX1(0,AMIN1(pos.k,nsubcell.k-1));

  return pos;
};

inline Real get_distance(my_real dr){
  Real dist = 0.0;
  dist += DIM_X==1 ? dr.i*dr.i : 0.0;
  dist += DIM_Y==1 ? dr.j*dr.j : 0.0;
  dist += DIM_Z==1 ? dr.k*dr.k : 0.0;

  dist = sqrt(dist);
  return dist;
};

inline void get_rgb_value(Real &r, Real &g, Real &b, Real value, Real min, Real max)
{
  Real f = (value - min)/(max - min);
  if (f < 0.) f = 0.;
  else if (f > 1.) f = 1.;

  Real a = (1. - f)/0.25;
  int  X = int(a);
  Real Y = a - X;
  r = g = b = 0.;

  if      (X == 0) {r=1; g=Y; b=0;}
  else if (X == 1) {r=1-Y; g=1; b=0;}
  else if (X == 2) {r=0; g=1; b=Y;}
  else if (X == 3) {r=0; g=1-Y; b=1;}
  else if (X == 4) {r=0; g=0; b=1;}

}

inline Real get_distance_2p(my_real a, my_real b){
  Real dist = 0.0;
  dist += DIM_X==1 ? pow(a.i-b.i,2) : 0.0;
  dist += DIM_Y==1 ? pow(a.j-b.j,2) : 0.0;
  dist += DIM_Z==1 ? pow(a.k-b.k,2) : 0.0;

  dist = sqrt(dist);
  return dist;
};

inline Real get_dot_product(my_real a, my_real b){
  Real val = 0.0;

  val += DIM_X==1 ? a.i*b.i : 0.0;
  val += DIM_Y==1 ? a.j*b.j : 0.0;
  val += DIM_Z==1 ? a.k*b.k : 0.0;

  return val;
};

inline my_real get_cross_product(my_real u, my_real v){
  my_real val; 

  val.i =  u.j*v.k - u.k*v.j;
  val.j = -u.i*v.k + u.k*v.i;
  val.k =  u.i*v.j - u.j*v.i;

  return val;
};

inline void my_shift_coordinate(my_real &p_coord, my_real coord, my_real domain, my_real box_l, my_real box_r)
{
  my_set_data (p_coord, coord);
#ifndef _MPI_
  #if PERI_DIM_X
  if ( coord.i > box_r.i + 1.e-20) p_coord.i -= domain.i;
  else if ( coord.i < box_l.i - 1.e-20) p_coord.i += domain.i;
  #endif
  #if PERI_DIM_Y
  if ( coord.j > box_r.j + 1.e-20) p_coord.j -= domain.j;
  else if ( coord.j < box_l.j - 1.e-20) p_coord.j += domain.j;
  #endif
  #if PERI_DIM_Z
  if ( coord.k > box_r.k + 1.e-20) p_coord.k -= domain.k;
  else if ( coord.k < box_l.k - 1.e-20) p_coord.k += domain.k;
  #endif
#endif
#ifdef _MPI_
  #if P_DIM_X == 1 && PERI_DIM_X == 1 && DIM_X == 1
  if (coord.i > box_r.i+1.e-20) p_coord.i -= domain.i;
  else if (coord.i < box_l.i-1.e-20) p_coord.i += domain.i;
  #endif
  #if P_DIM_Y == 1 && PERI_DIM_Y == 1 && DIM_Y == 1
  if (coord.j > box_r.j+1.e-20) p_coord.j -= domain.j;
  else if (coord.j < box_l.j-1.e-20) p_coord.j += domain.j;
  #endif
  #if P_DIM_Z == 1 && PERI_DIM_Z == 1 && DIM_Z == 1
  if (coord.k > box_r.k+1.e-20) p_coord.k -= domain.k;
  else if (coord.k < box_l.k-1.e-20) p_coord.k += domain.k;
  #endif
#endif
}
inline my_real my_check_position(my_real coord, my_real domain, my_real box_l, my_real box_r)
{
  my_real p_coord;
  my_set_data (p_coord, coord);
#ifndef _MPI_
  #if PERI_DIM_X
  if ( coord.i > box_r.i + 1.e-20) p_coord.i -= domain.i;
  else if ( coord.i < box_l.i - 1.e-20) p_coord.i += domain.i;
  #endif
  #if PERI_DIM_Y
  if ( coord.j > box_r.j + 1.e-20) p_coord.j -= domain.j;
  else if ( coord.j < box_l.j - 1.e-20) p_coord.j += domain.j;
  #endif
  #if PERI_DIM_Z
  if ( coord.k > box_r.k + 1.e-20) p_coord.k -= domain.k;
  else if ( coord.k < box_l.k - 1.e-20) p_coord.k += domain.k;
  #endif
#endif
#ifdef _MPI_
  #if P_DIM_X == 1 && PERI_DIM_X == 1 && DIM_X == 1
  if (coord.i > box_r.i+1.e-20) p_coord.i -= domain.i;
  else if (coord.i < box_l.i-1.e-20) p_coord.i += domain.i;
  #endif
  #if P_DIM_Y == 1 && PERI_DIM_Y == 1 && DIM_Y == 1
  if (coord.j > box_r.j+1.e-20) p_coord.j -= domain.j;
  else if (coord.j < box_l.j-1.e-20) p_coord.j += domain.j;
  #endif
  #if P_DIM_Z == 1 && PERI_DIM_Z == 1 && DIM_Z == 1
  if (coord.k > box_r.k+1.e-20) p_coord.k -= domain.k;
  else if (coord.k < box_l.k-1.e-20) p_coord.k += domain.k;
  #endif
#endif
  return p_coord;
}

inline Real my_get_ramping_factor (Real fac0, Real fac1, Real t, Real t0, Real t1)
{
  Real fac = 0.;
  if (t0>=t1) cout<<"<<<< ERROR in ramping time provided!!!!"<<endl;

  if (t <= t0){
    fac = fac0;
  }else if (t <= t1){
    fac = fac0 + (t-t0)/(t1-t0)*(fac1-fac0);
  }else{
    fac = fac1;
  }
  return fac;
}

// compress infor
inline unsigned int compress(int i, int j, int k, int l)
{
  if (i < 0 || i > (1 << I_BITS) - 1) {
    cout << "compress error: i is out of range" << endl;
    return 0xffffffff;
  }
  if (j < 0 || j > (1 << J_BITS) - 1) {
    cout << "compress error: j is out of range" << endl;
    return 0xffffffff;
  }
  if (k < 0 || k > (1 << K_BITS) - 1) {
    cout << "compress error: k is out of range" << endl;
    return 0xffffffff;
  }
  if (l < 0 || l > (1 << L_BITS) - 1) {
    cout << "compress error: l is out of range" << endl;
    return 0xffffffff;
  }

  return (i << I_SHIFT) | (j << J_SHIFT) | (k << K_SHIFT) | (l << L_SHIFT);
};
// uncompress infor
inline void uncompress(unsigned int data, int &i, int &j, int &k, int &l)
{
  i = (data >> I_SHIFT) & I_MASK;
  j = (data >> J_SHIFT) & J_MASK;
  k = (data >> K_SHIFT) & K_MASK;
  l = (data >> L_SHIFT) & L_MASK;
};

// This function returns true if the first pair is "less"
// than the second one according to some metric
// In this case, we say the first pair is "less" if the first element of the first pair
// is less than the first element of the second pair
inline bool pairCompare(const std::pair<Real, my_int>& firstElem, const std::pair<Real, my_int>& secondElem) {
  return firstElem.first < secondElem.first;
};

struct Greater
{
    template<class T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};
#endif
