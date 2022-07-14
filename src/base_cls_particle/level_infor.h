#ifndef LEVEL_INFO_H
#define LEVEL_INFO_H
#include <cmath>
#include "glbcls.h"
#include "glbfunc.h"
#include "Mypool.h"
#include "boost/multi_array.hpp"
#ifdef _MPI_
#include "graph.h"
#endif

using namespace std;
using namespace boost::mpi;
using namespace boost;

typedef boost::multi_array<DTAG, 3> array_tag;
typedef boost::multi_array<p_Cell_list, 3> array_cell_list;
#ifdef _MPI_
typedef boost::multi_array<p_Color_list, 3> array_color_list;
#endif
#ifdef _NARROW_BAND_GRAPH_
typedef boost::multi_array<tbb::atomic<int>, 3> array_tag_atomic;
#endif

class Level_info{
  int                total_num_particle;
  int                total_nun_color;

public:

// variables
  int                level;
  int                Lmin;
  int                Lmax;
  int                total_num_cell;
  int                num_leaf;
  int                num_leaf_particle;
  my_int             num_cell;            //number of cells in each directions
  my_int             glbl_num_cell;
  my_int             glbl_cell_start;
  my_int             glbl_cell_end;
  my_int             cell_start;
  my_int             cell_end;
#if PERI_DIM != 0 || SYM_DIM != 0
  my_int             num_cell_margin;     //number of cells in periodical boundary buffer
  Mypool<Particle>   *p_particlepool;
#endif
  my_real            dcell;               //cell size in each directions
  my_real            domain;              //domain size in each direction
  my_real            box_l;               //coordinate of the left corner of the domain
  my_real            box_r;
  Real               scale;
  array_tag          exist_cell_list;     //record the cell list
  array_tag          exist_leaf_particle;
#ifdef _MPI_
  my_real            local_box_l;         //coordinate of the local left corner of the domain
  my_real            local_box_r;
  // for input and output matrix with mpi
  array_tag          exchange_exist_cell_list;
  array_tag          exist_exchange_particle;
#ifdef _NARROW_BAND_GRAPH_
  array_tag_atomic   list_is_narrow_band;
#endif
  Mypool<Color_list> *p_color_listpool; 
  array_color_list   color_list;
#endif
#ifdef _SCLL_
  my_real            dsubcell;
  my_int             nsubcell;
  int                total_num_subcell;
  my_int             subcell_start;
  my_int             subcell_end;
#endif
  concurrent_vector  <p_Particle> leaf_particle;
  Mypool<Cell_list>  *p_cell_listpool; 
  // every level Cell_listpool point to the main memory pool
  array_cell_list    table_cell_list;
  //point to different cell list;
  Level_info(){};  

// functions
  void Initialize (int i, SPH *sph, communicator &world);
  void Reinitialize_cell_list(communicator &world, SPH* sph);
  void Allocate_memory();
  void Initial_cell_list_and_reset_tags();
  int  Get_level();
  void Reset_cell_list_info(communicator &world);
  void Reset_tags(communicator &world);
  void Update_exist_status();
  void Update_leaf_particle_status();
  void Update_cell_list(Particle *current_particle, communicator &world);
  void Update_every_level_info(int flag, p_Level_info *level_infos, int child_level, communicator &world);
#ifndef _SCLL_
  void Refresh_neighbor_info(SPH *sph, Particle *current, my_int cell_id, int flag);
#else
  void Refresh_neighbor_info(SPH *sph, Particle *current, my_int scell_id, int flag);
#endif
  void Interaction(Particle *current, Particle *neighbor, my_real dr, Real dist, int flag, int status);
  void Output_level(communicator &world, char *filename, int n);
#ifdef _SCLL_
  void Set_sub_cell_range();
  void Build_subcell_list (communicator &world);
#endif
#if SYM_DIM != 0
  void Construct_symmetric_BC_particles(concurrent_vector  <p_Particle> &particle_sym);
#endif
#if PERI_DIM != 0
  void Construct_periodical_BC_particles_1(concurrent_vector  <p_Particle> &particle_peri);
  #ifdef _MPI_
  void Construct_periodical_BC_particles_2(concurrent_vector  <p_Particle> &particle_peri, communicator &world);
  #endif
#endif

#ifdef _MPI_
  #if defined(_DIST_NUMBER_) || defined(_WEIGHTED_PARTITION_)
    #ifndef _SCLL_
  void Get_patitioning_mass(SPH *sph, Particle_base *current, my_int cell_id, int flag);
    #else
  void Get_patitioning_mass(SPH *sph, Particle_base *current, my_int scell_id, int flag);
    #endif
  #endif
  #ifdef _NARROW_BAND_GRAPH_
  void Reset_narrow_band_list(communicator &world);
  void Find_narrow_band_cells();
  #endif
  void Allocate_color_list(communicator &world);
  void Reset_color_list();
  void Get_local_cell_start_end(SPH *sph);
  void Initialize_color_infor();
  void Calculate_local_tree_start_end(my_int &start, my_int &end);
  void Shift_cell_index(int &t, int &s, int &m, int i, int j, int k);
  void Init_exchange_status();
  void Update_exchange_status(p_Level_info *level_infos, int parent_level);
  void Merge_exchange_status();
  void Update_every_level_color_info(p_Level_info *level_infos, int child_level);
  #if PERI_DIM != 0
  void Update_PBC_color_info();
  void Update_cell_list_and_shift_coord(Particle *current_particle, communicator &world);
  #endif
  void Traverse_tree_and_construct_graph(Graphcls *sph_graph);
  void Traverse_tree_and_construct_graph_modified(Graphcls *sph_graph);
  void Exchange_graph_topology(communicator &world, serialization_vector <unsigned int> &exchange_vector);
  void Exchange_buffer_between_pairs(communicator &world,int target_processor, concurrent_vector  <p_Particle> &particle_ghost, int flag);
  void Output_color_list(communicator &world, char *filename, int n);
#endif
#ifdef _MEM_CHECK_
  void Check_memory_consumption(communicator &world, long long &mem_cell_temp, long long &mem_list_temp); 
#endif
};
#endif
