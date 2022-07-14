#ifndef _LEVEL_SET_PACKAGE_
#define _LEVEL_SET_PACKAGE_
#include "glbcls.h"
#include "lset_cell.h"
#include "boundary_condition.h"
#include "geometry_calculator.h"
#include <vector>
#include <array>

using namespace std;
using namespace boost::mpi;
using namespace boost;

typedef boost::multi_array<p_Levelset_cell, 3> array_lset_cell;
class   Levelset_package{

private:

public:
  // variable
  int level;
  my_int index;
  Levelset_cell lset_cell[ICPX][ICPY][ICPZ];
  array_lset_cell p_cell;

  Real total_mass; // calculate mass for SPH
  Real total_mass_surface;
  Real total_mass_segment;
  Real total_volume;
  Real total_area;
  Real total_length;
  Real maximum_phi;
  Real maximum_dl;
  Real middel_dl; 
  Real minimum_dl; 
  Real maximum_curv;
  Real minimum_curv;
  Real maximum_psi;
  int tag_convergence;
  int n_interface_cell;

  #ifdef _MPI_
  int color;
  #endif
  
  Levelset_package(){};
  // functions
  void Initialize (p_Levelset_levelinfo lset_level_info, int pkg_i, int pkg_j, int pkg_k);
  void Define_levelset_full_field(int pkg_i, int pkg_j, int pkg_k, my_real dpkg, my_real dcell, Levelset *level_set);
  my_real get_cell_position (int pkg_i, int pkg_j, int pkg_k, int cell_i, int cell_j, int cell_k, my_real dpkg, my_real dcell, my_real box_l);
  my_real get_cell_position_corner (int pkg_i, int pkg_j, int pkg_k, int cell_i, int cell_j, int cell_k, my_real dpkg, my_real dcell, my_real box_l);
  my_real get_pkg_position_corner (int pkg_i, int pkg_j, int pkg_k, my_real dpkg, my_real box_l);
  my_real get_pkg_position (int pkg_i, int pkg_j, int pkg_k, my_real dpkg, my_real box_l);
  void Copy_phi_value();
  void Set_cell_topology(p_Levelset_levelinfo lset_level_info);
  void Boundary(p_Levelset_levelinfo lset_level_info);
  void Boundary_curv(p_Levelset_levelinfo lset_level_info);
  void Get_normal(p_Levelset_levelinfo lset_level_info);
  void Get_BC_normal(p_Levelset_levelinfo lset_level_info);
  void Get_curvature(p_Levelset_levelinfo lset_level_info);
  void Clean_curvature();
  void Get_interface_tag(p_Levelset_levelinfo lset_level_info);
  void Get_extended_interface_tag(p_Levelset_levelinfo lset_level_info);
  void Get_singularity_tag(Levelset *level_set, my_real *Pos_singularity);
  void Get_segment_tag(Levelset *level_set);
  bool is_Segment_intersec_AABB (my_real p1, my_real p2, my_real min, my_real max);
  int  Return_package_interface_tag(p_Levelset_levelinfo lset_level_info);
  bool Is_cut_cell (int i, int j, int k);
  bool is_in_current_cell (my_real pos, my_real cell_corner, my_real dcell);
  bool is_in_current_pkg (my_real pos, my_real pkg_corner, my_real dpkg);
  void Get_phi_at_corners(int i, int j, int k, std::array<Real,8> &phi_corner_values);
  void Get_global_scale(p_Levelset_levelinfo lset_level_info);
  void Extend_psi (Levelset *level_set);
  void Update_extend_psi (Levelset *level_set);
  void Update_extend_psi_reinitialize (Levelset *level_set);
  void Reinitialize_psi_backup(Levelset *level_set);
  void Reinitialize_psi_reset(Levelset *level_set);
  void Reinitialize_psi_increment(Levelset *level_set);
  bool Check_covergence_psi(Levelset *level_set);
  void Get_psi_max(Levelset *level_set);
  void Redistribute_curv(Levelset *level_set);
  void Extend_curv_backup(p_Levelset_levelinfo lset_level_info);
  void Extend_curv(p_Levelset_levelinfo lset_level_info);
  void Update_extend_curv(Levelset *level_set);
  void Extend_curvature_for_smooth (Levelset *level_set);
  void Update_curvature_for_smooth (Levelset *level_set);
  void Update_extend_curv_for_smoothing(p_Levelset_levelinfo lset_level_info);
  bool Check_covergence_curv(Levelset *level_set);
  void Copy_curvature(Levelset *level_set);
  void Calculate_volume_mass(Levelset *level_set);
  void Output_level_set(Levelset *level_set, communicator &world, char *filename, int n);
};

#endif
