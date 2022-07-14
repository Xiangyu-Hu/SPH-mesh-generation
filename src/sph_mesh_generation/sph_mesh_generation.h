#ifndef _SPH_MG_H
#define _SPH_MG_H
#include "sph.h"
#include "mesh_binary_writer.h"
#include "particle_binary_writer.h"
#include "particle_mesh_generation.h"
#include "graph_mesh.h"
#include "level_set.h"

class SPH_mesh_generation : public SPH {

public:

  // variables
  int  firststep;
  Real gamma;
  Real WVT_coeff;
  Real APD_coeff;

  Real damp_coeff_start;
  Real damp_coeff_end;
  Real damp_ramping_start;
  Real damp_ramping_end;

  Real nu_coeff_start;
  Real nu_coeff_end;
  Real nu_ramping_start;
  Real nu_ramping_end;
  my_real body_force;

  Real angle_thres;
  Real t_for_post;
  int  n_post;

  #ifdef _MPI_
  Real sigma_coeff_start;
  Real sigma_coeff_end;
  Real sigma_ramping_start;
  Real sigma_ramping_end;
  #endif

  Real tet_delete_thres;

  #ifdef _READ_SDF_
  char   sdf_file[256];
  my_int cells_to_cut;
  #endif
  
  int v_reini_start;
  int v_reini_end;
  int v_reini_change;
  
//   Real max_h;
//   Real min_h;
  
  int max_num_iteration;
  
  int local_id_record;
  
  int total_num_particle_mesh;
  int glbl_total_num_particle_mesh;
  int current_glbl_total_num_particle_mesh;
  
  my_int resolution;
  
  int  count_full;
  int  max_count_full;
  int  max_count_relax;
  Real fill_coeff;
  bool fill;
  int  current_stage;

  int total_num_particle_surface;
  int glbl_total_num_particle_surface;
  int current_glbl_total_num_particle_surface_old;
  int current_glbl_total_num_particle_surface;
  tbb::atomic<bool> find_surface_particle;
  
  int total_num_particle_singularity;
  int glbl_total_num_particle_singularity;
  int current_glbl_total_num_particle_singularity;
  std::vector<int> num_particles_per_singularity;
  tbb::atomic<bool> find_singularity_particle;

  int total_num_particle_segment;
  int glbl_total_num_particle_segment;
  int current_glbl_total_num_particle_segment;
  tbb::atomic<bool> find_segment_particle;
  
  concurrent_vector <p_Particle> particle_singularity;
  concurrent_vector <p_Particle> particle_segment;
  concurrent_vector <p_Particle> particle_surface;
  concurrent_vector <p_Particle> ghost_particle;
  concurrent_vector <p_Particle> particle_bound;
  concurrent_vector <p_Particle> dummy_particle;
  concurrent_vector <p_Particle> real_particle;
  
  Levelset  level_set;

  Graphmeshcls mesh;

  Mesh_binary_writer     mesh_writer;

  Particle_binary_writer particle_writer;
// functions
  SPH_mesh_generation                       (Initialization &Ini, communicator &world);

  void Set_timestep                         (communicator &world, int stage);
  void Load_case                            (communicator &world);
  void Load_rst_case                        (communicator &world);
  void Load_data_for_post                   (communicator &world);
  void Initialize_case                      (communicator &world);
  void Get_num_particles                    (communicator &world);
  void Pre_run                              (communicator &world);
  void Recalculate_local_particle_number    (communicator &world);
  void Allocate_particles                   (communicator &world);
  void Asign_local_particle_index           (communicator &world);
  void Set_local_ghost_index                (communicator &world);
  void Random_particle_distribution         (communicator &world);
  void SPH_run_mesh_generation              (communicator &world);
  void Get_particle_info                    (communicator &world);
  void Set_particle_scale                   (communicator &world);
  void Find_and_map_particle                (communicator &world);
  void Check_and_pop_surface_particle       (communicator &world);
  void Check_total_number_of_particle       (communicator &world);
  void Remap_particle                       (communicator &world);
  void Reset_particle_array                 (communicator &world);
  void Remap_this_particle_to_positive_phase(p_Particle particle);
  void Prepare_for_the_current_step         (communicator &world);
  void Reset_particle_force                 (communicator &world, int stage);
  void Get_kernel_summation                 (communicator &world);
  void Accumulate_particle_force            (communicator &world, int stage);
  void Map_surface_particle_force           (communicator &world);
  void Update_velocity_half_step            (communicator &world);
  void Update_coord_full_step               (communicator &world);
  void Free_dummy_particles                 (std::vector<p_Particle> &particle_total, communicator &world);
  void Reconstruce_dummy_particles          (std::vector<p_Particle> &particle_total, communicator &world);
  void Post_processing                      (int _nout, communicator &world);
  void Output_plt_file                      (int n, int flag, communicator &world);
  void Output_paraview_file                 (int n, int flag, communicator &world);
  #if SYM_DIM != 0
  void Output_SBC_plt_file                  (int n, communicator &world);
  #endif
  void Output_simulation_infor              (int n, communicator &world);
  void Save_restart_file                    (communicator &world);
  void Save_particle_post_file              (int n, communicator &world);
  #ifdef _MPI_
  void Random_particle_distribution_within_R(Real R, communicator &world);
  void Get_particle_info_for_ghosts         (communicator &world);
  #if defined (_CVP_LSET_) || defined (_READ_VP_)
  void Random_particle_distribution_within_VPi (communicator &world);
  #endif
  #endif
};

#endif
