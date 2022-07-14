#ifndef _VORONOI_H_
#define _VORONOI_H_
#include "glbcls.h"
#include "voronoi_base.h" 
#include "Mypool.h"
#include "graph.h"
#include "visualize.h"

class Voronoi : public Voronoi_base{

  public:
  // variables
  int      num_partition;
  int      iterate_num;
  int      convergenced;
  int      Max_iterate_number;
  int      num_sph_particle;
  Real     error_tolerance;
  Real     glbl_timestep;
  Real     glbl_timestep_CVT;
  Real     mass_target;
  Real     total_mass;
  Real     scale_coeff;
  Real     glbl_total_mass;
  Real     relax_ratio;
  Real     enlarge_ratio;
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
  Real     x_min, x_max;
  Real     y_min, y_max;
  Real     z_min, z_max;
  bool     P_X, P_Y, P_Z;

  p_Particle                        *sph_particle;

  concurrent_vector <p_Vo_particle> vo_particle;
  concurrent_vector <p_Vo_particle> vo_particle_ghost;// for the boundary condition 
  concurrent_vector <p_Vo_particle> vo_particle_record;
  
  Mypool<Vo_particle>               Vo_particlepool;
  Mypool<p_Particle>                p_particlepool;
  Graphcls                          vo_graph;         // graph for partitioning
//  communicator                      world;
  Visualize                         visual;           // visualization tools
#if PERI_DIM != 0
  my_int   cross_boundary;
#endif

  // functions
  Voronoi (){};
  void      Initialize(SPH *sph, int flag, communicator &world);
  my_real   Get_coord_of_VP_at_iRank(int iRank);
  Real      Get_h_min_of_VP_at_iRank(int iRank);
  void      Set_VP_coords_and_hmin(std::vector<my_real> &vp_coords, std::vector<Real> &vp_scale, communicator &world);
  Real      Hex_close_packing_inbox(my_real box_l_, my_real box_, int dimX, int dimY, int dimZ, communicator &world);
  void      Random_particle_distribution(communicator &world);
  void      Random_particle_distribution_inbox(my_real box_l_, my_real box_, communicator &world);
  int       Get_CPU_id(my_real p_coord);
  void      Partitioning(SPH *sph, communicator &world);
  void      Set_target_mass(communicator &world);
  void      Get_VP_info_and_VD_generation(communicator &world);
  void      Update_position(communicator &world);
  void      Update_position_CVT(communicator &world);
  void      Update_vp_position_mean_velo(SPH *sph, communicator &world);
  void      Update_vp_position_mass_center(SPH *sph, communicator &world);
  void      Sync_VP_positions(SPH *sph, communicator &world);
  void      Shift_coordinate_constrained(my_int dim, SPH *sph, communicator &world);
  void      Check_partition_error(int &iter, communicator &world);
  void      Construct_graph(communicator &world);
  void      Add_edges(communicator &world);
  void      Output_vtk_file(int i_iter, int n, communicator &world);
  void      Output_plt_file(int i_iter, int n, communicator &world);
  void      Output_ghost_vtk_file(int i_iter, int n, communicator &world);
  void      Output_pov_ray_file(int i_iter, int n, communicator &world);
#if PERI_DIM != 0
  void      Shift_SPH_particle_position(SPH *sph, communicator &world);
  void      Shift_back_SPH_particle_position(SPH *sph, communicator &world);
#endif
};

#endif