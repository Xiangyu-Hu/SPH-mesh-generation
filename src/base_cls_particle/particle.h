#ifndef PARTICLE_H
#define PARTICLE_H
#include <cmath>
#include "glbcls.h"
#include "Mypool.h"

using namespace std;

class Particle_base{
private:

public:

// variables
  int      tag;
  int      id;
  int      color;
  int      level;        //tree level
  Real     mass;
  Real     p_mass;       // for load balance
  Real     timestep;
  Real     rho;          //density
  Real     P;            //pressure
  Real     vol;          //volumn
  Real     h;            //smooth length
  Real     scale;        //particle scale
  Real     c;
  my_real  coord;        //coordinate
  my_real  v;            //material velocity
  my_real  F1;           //force
  my_real  F2;           //force
  my_real  a;            //acceleration

  Particle_base(){};

// functions

  Real    Kernel_function(Real dist, Real h);
  Real    Kernel_function_1D(Real dist, Real h);
  Real    Kernel_function_2D(Real dist, Real h);
  Real    Kernel_function_3D(Real dist, Real h);
  Real    Derivative_kernel_function(Real dist, my_real dr, Real h);
  Real    Derivative_kernel_function_1D(Real dist, my_real dr, Real h);
  Real    Derivative_kernel_function_2D(Real dist, my_real dr, Real h);
  Real    Derivative_kernel_function_3D(Real dist, my_real dr, Real h);
  Real    Derivative_h_kernel_function(Real dist, Real h);
  void    Get_velocity_for_TGV();
  void    Set_particle_scale(SPH *sph, int &flag);
  void    Reinitialize_particle_info(int flag);
  void    Get_energy(Real &iEk, Real &iEp, Real &iEf);
#ifdef _MPI_
  void    Reset_color(int current_color);
  void    Shift_coordinate(my_real domain, my_real box_l, my_real box_r);
  void    Shift_coordinate_constrained(my_real domain, my_real box_l, my_real box_r, my_int dim);
  void    Get_patitioning_mass(SPH *sph, int flag, int n_neighbor);
  void    Normalize_particle_mass(SPH *sph);
#ifdef _WEIGHTED_PARTITION_
  void    Get_weighted_patitioning_mass(SPH *sph, int flag, int n_neighbor);
#endif
#endif
};

#endif
