#ifndef SPH_H
#define SPH_H
#include "glbcls.h"
#include "Mypool.h"
#include "parallelfunc.h"
#ifdef _INCPRS_
#include "particle_incprs.h"
#endif
#ifdef _CPRS_
#include "particle_cprs.h"
#endif
#ifdef _GSPH_
#include "particle_gsph.h"
#endif
#ifdef _ALE_
#include "particle_ale.h"
#endif
#ifdef _MESH_GENERATION_
#include "particle_mesh_generation.h"
#endif
#include "cell_list.h"
#include "level_infor.h"
#include "system.h"
#include "visualize.h"
#ifdef _MPI_
#include "graph.h"
#include "voronoi.h"
#endif

using namespace std;
using namespace boost;
using namespace boost::serialization;
using namespace boost::mpi;
using namespace tbb;

//-----------------------------------------
class SPH{

public:
// globle variables
  int      glbl_total_num_particle;
  int      glbl_total_num_color;
  int      num_level;
  int      Lmax;
  int      Lmin;
  int      need_for_refresh;
  my_int   glbl_num_particle;
  Real     glbl_timestep;      //current timestep size
  Real     glbl_max_v;
  vector<int> nparticle_in_each_cpu;

// local variables
  int      iterate_num;        //current iterate number
  int      i_iter;
  int      total_iterate_num;
  Real     local_timestep;     //current timestep size
  Real     run_time;           //current iterate time
  Real     timestep_shift;     //timestep size for redistributing particles

// globle variables
  Real     ini_vol;     //initial volumn
  Real     ini_h;       //initial smooth length
  Real     ini_scale;   //initial scale for base level
  Real     ini_rho;     //initial density
  Real     ini_mass;
  Real     P0;
  Real     max_v;
  Real     max_scale;
  Real     min_scale;
  Real     c;           //sound speed
  Real     glbl_Ek;     //kinetic energy
  Real     glbl_Etotal; //total energy
  Real     glbl_force_error;
  my_int   ini_num_cell;//number of cells at beginning
  my_real  domain;      //domain size in each direction
  my_real  box_l;       //coordinate of the left corner of the domain
  my_real  box_r;
  my_real  particle_box;
  my_real  particle_box_l;
  my_real  particle_box_r;
  my_real  p_size;      //initial_particle_size
  my_real  glbl_total_force;
  my_real  v_avg;
  
// local variables
  int      total_num_particle;
  int      num_buffer_particle;
  Real     Ek;          //kinetic energy
  Real     Ep;          //gravitational potential energy
  Real     Ef;          //pressure potential energy
  Real     Etotal;      //total energy

  p_Particle          *particle;
  Mypool<Particle>    particlepool;
  Mypool<Cell_list>   cell_listpool;
  p_Level_info        *level_info;

  Visualize          visual;           // visualization tools

  Real     *time_for_graph;
  Real     *time_for_partition;
  Real     *time_for_simulation;
  Real     *time_for_force_calculation;
  Real     *time_for_density;
  Real     *time_for_rest;
  Real     *time_for_update;
  Real     *time_for_reset_neighbor;
  Real     *time_for_mapping;
  Real     *time_for_bulid_local_map;
  Real     *time_for_map_particle_to_tree;
  Real     *time_for_refresh_neighbor;
  Real     *time_for_clear_cell_list;
  Real     *time_for_release_particle;
  Real     *time_for_update_every_level_info;
  Real     *time_for_communication0;
  Real     *time_for_communication1;
  Real     *time_for_communication2;
  Real     *time_for_communication_partition;
  Real     *time_for_communication_flag1;
  Real     *time_for_communication_flag2;
  Real     *time_for_exchangebuffer1;
  Real     *time_for_exchangebuffer2;
  Real     *time_for_total;
  timer    time_graph;
  timer    time_partition;
  timer    time_simulation;
  timer    time_force_calculation;
  timer    time_density;
  timer    time_rest;
  timer    time_update;
  timer    time_reset_neighbor;
  timer    time_mapping;
  timer    time_bulid_local_map;
  timer    time_map_particle_to_tree;
  timer    time_refresh_neighbor;
  timer    time_clear_cell_list;
  timer    time_release_particle;
  timer    time_update_every_level_info;
  timer    time_communication;
  timer    time_communication_partition;
  timer    time_exchangebuffer;
  timer    time_communication_flag;
  timer    time_total;
#if PERI_DIM != 0
  int      cross_boundary;
  concurrent_vector  <p_Particle> particle_peri; // periodical bc particles
#endif
#if SYM_DIM != 0
  concurrent_vector  <p_Particle> particle_sym; // symmetric bc particles
#endif
#ifdef _MPI_
  int       icount;
  int       target_processor;
  int       need_for_partition;
  int       need_for_partition_count;
  int       need_for_rebuild_graph;
  int       need_for_rebuild_local_tree;
  Real      exchange_local;
  Real      exchange_local_old;
  Real      error_local;
  Real      error_global;
  Real      error_tolerance;
  Real      exchange_total;
  Real      exchange_max;
  Real      exchange_min;
  Real      exchange_avg;
  Real      glbl_total_mass;
  Real      total_mass;
  Real      total_mass_old;
  Real      imbalance_local;
  Real      imbalance_global;
  Real      imbalance_tolerance;
  Real      imbalance_min;
  Real      imbalance_max;
  Real      imbalance_avg;
  Real      imbalance_total;
  Real      migrate_total;
  Real      migrate_local;
  Real      migrate_max;
  Real      migrate_min;
  Real      migrate_avg;
#ifdef _WEIGHTED_PARTITION_
  Real      weight;
  Real      t_force;
  Real      t_neighbor;
  Real      t_simulation;
  Real      t_force_old;
  Real      t_neighbor_old;
  Real      t_simulation_old;
  Real      glbl_total_dist_number;
  Real      glbl_total_neighbor_size;
#endif
  my_real   mass_center;
  my_real   glbl_mass_center;
  my_real   local_box_l; //coordinate of the left corner of the local tree
  my_real   local_box_r; //coordinate of the right corner of the local tree 
  concurrent_vector  <p_Particle> particle_buffer;
  concurrent_vector  <p_Particle> particle_ghost; // particles that are found to be buffers for neighboring domain
  concurrent_vector  <int>particle_ghost_start; // the starting ghost particle in each time that data exhcange happens
  Graphcls           sph_graph;   // graph for buffer particle communication
  Voronoi            sph_voronoi; // used for voronoi partitioning
  Mypool<Color_list> color_listpool;
#endif
#ifdef _MEM_CHECK_
  Real     local_mem;
  Real     mem_particle;
  Real     mem_cell;
  Real     mem_list;
  Real     glbl_mem_max;
  Real     glbl_mem_total;
#endif
  Real     safe_guard;

  SPH(){};

// functions
  void Output_pov_ray_file(int n, communicator &world);
  void Output_every_level_dat(int n, communicator &world);
  void Output_plt_file(int n, int flag, communicator &world);
  void Output_restart_file();
  void Initialize_system(communicator &world);
  void Define_level_infor (communicator &world);
  void Define_level_infor_tight (communicator &world);
  void Set_current_time(Real Time);
  void Set_timestep (communicator &world);
  Real Get_timestep ();
  void Reset_timestep(Real timestep_new);
  void Reset_cell_list_info(communicator &world);
  void Map_the_particle_to_tree(communicator &world);
  void Update_every_level_info(int flag, communicator &world);
  void Reset_particle_neighbor_info(communicator &world);
  void Refresh_neighbor_info(int flag, communicator &world);
  void Refresh_particle_scale(communicator &world);
  void Refresh_particle_infor(communicator &world, int flag1, int flag2);
  void Add_coarse_level(communicator &world);
  void Add_finer_level(communicator &world);
  void Read_infile();
  void Get_energy(communicator &world);
  void Check_total_force(communicator &world);
  void Output_runtime_info(int flag, communicator &world);
#ifdef _SCLL_
  void Build_every_level_subcell_list(communicator &world);
#endif
#if SYM_DIM != 0
  void Release_memory_of_sbc_paticle(communicator &world);
  void Construct_symmetric_BC_particles(communicator &world);
  void Refresh_sbc_particle_info(communicator &world, int flag);
  void Map_the_sym_particle_to_tree(communicator &world);
#endif
#if PERI_DIM != 0
  void Release_memory_of_pbc_paticle(communicator &world);
  void Construct_periodical_BC_particles(communicator &world, int flag);
  void Refresh_pbc_particle_info(communicator &world, int flag);
  void Map_the_peri_particle_to_tree(communicator &world);
#endif
#ifdef _MPI_
  #if defined(_DIST_NUMBER_) || defined(_MASS_) || defined(_WEIGHTED_PARTITION_)
  Real Get_total_p_mass_local(communicator &world);
  #endif
  #ifdef _WEIGHTED_PARTITION_
  Real Get_global_neighbor_size(communicator &world, int flag);
  Real Get_global_dist_number(communicator &world);
  void Get_weighted_patitioning_mass(int flag, communicator &world);
  #endif
  void Calculate_total_partitioning_mass(communicator &world);
  void Get_patitioning_mass(int flag, communicator &world);
  void Normalize_particle_mass(communicator &world);
  void Release_memory_of_buffer_paticle(communicator &world);
  void Release_memory_of_ghost_paticle(communicator &world);
  void Update_exchange_infor(communicator &world);
  void Construct_graph(communicator &world);
  void Construct_graph_modified(communicator &world);
  void Partitioning(communicator &world);
  void Check_need_for_partition(communicator &world);
  void Update_vp_position_mean_velo(communicator &world);
  void Update_vp_position_mass_center(communicator &world);
  void Migrate_SPH_particle(communicator &world);
  void Data_communication(int flag, communicator &world);
  void Communication_between_processor_pair(int flag, int iter, communicator &world);
  void Migration_between_processor_pair(int iter, concurrent_vector <p_Particle> &incoming_particles, communicator &world);
  void Map_the_buffer_particle_to_tree(communicator &world);
  void Find_bound_in_coarsest_level (communicator &world);
#if PERI_DIM != 0
  void Check_PBC_particle_positions (communicator &world);
#endif
  void Traverse_and_construct_local_tree (communicator &world, int flag);
  void Build_local_tree(communicator &world, int flag);
  void Output_global_info(int flag, communicator &world);
  void Output_color_list_dat(int n, communicator &world);
#endif
#ifdef _MEM_CHECK_
  void Check_memory_consumption(communicator &world);
#endif
};
#endif
