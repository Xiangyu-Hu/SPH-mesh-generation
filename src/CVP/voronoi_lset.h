#ifndef _VORONOI_LSET_H_
#define _VORONOI_LSET_H_
#include "glbcls.h"
#include "voronoi_base.h" 
#include "Mypool.h"

class Voronoi_lset : public Voronoi_base{

  public:
  // variables
  int      total_num_positive_pkg;
  int      num_partition;
  int      iterate_num;
  int      convergenced;
  int      Max_iterate_number;
  Real     error_tolerance;
  Real     glbl_timestep;
  Real     glbl_timestep_CVT;
  Real     mass_target;
  Real     total_mass;
  Real     scale_coeff;
  Real     glbl_total_mass;
  Real     relax_ratio;
  Real     total_interface_area;
  Real     total_energy;
  Real     error_max;
  Real     error_record;

  my_real  domain;
  my_real  box_l;       //coordinate of the left corner of the domain
  my_real  box_r;
  my_real  particle_box;
  my_real  particle_box_l;
  my_real  particle_box_r;
  my_real  dpkg;
  Real     x_min, x_max;
  Real     y_min, y_max;
  Real     z_min, z_max;
  bool     P_X, P_Y, P_Z;

  p_Vo_particle * _vo_particle;
  concurrent_vector <p_Vo_particle> vo_particle;
  concurrent_vector <p_Vo_particle> vo_particle_ghost;// for the boundary condition 
  concurrent_vector <p_Vo_particle> vo_particle_record;
  
  std::vector<p_Levelset_package>   lset_positive_pkg;
  concurrent_vector<int>            contain_cutcell;

  Mypool<Vo_particle>               Vo_particlepool;

  // functions
  Voronoi_lset (){};
  void      Initialize(Levelset *level_set, int flag, communicator &world);
  my_real   Get_coord_of_VP_at_iRank(int iRank);
  void      Random_particle_distribution(Levelset *level_set, communicator &world);
  int       Get_CPU_id(my_real p_coord);
  void      Partitioning_for_initial_VP_distribution(Levelset *level_set, communicator &world);
  void      Set_target_mass(Levelset *level_set, communicator &world);
  void      Get_vp_position_and_hmin (std::vector<my_real> & vp_coords, std::vector<Real> &vp_scale, communicator &world);
  void      Ghost_VP_generation(Levelset *level_set, communicator &world);
  void      Get_VP_info_and_VD_generation(communicator &world);
  void      Update_position(communicator &world);
  void      Update_position_CVT(communicator &world);
  void      Check_partition_error(int &iter, communicator &world);
  void      Output_vtk_file(int i_iter, int n, communicator &world);
  void      Output_ghost_vtk_file(int i_iter, int n, communicator &world);
};

#endif