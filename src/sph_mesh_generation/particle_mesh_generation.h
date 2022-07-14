#ifndef PARTICLE_MG_H
#define PARTICLE_MG_H
#include <cmath>
#include "Mypool.h"
#include "particle.h"
#include "glbcls.h"

using namespace std;

class Particle_mesh_generation : public Particle_base{
private:
#ifdef _MPI_
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & tag;
    if (tag == 0){
      ar & id;
      ar & type;
      ar & h;
      ar & scale;
      ar & coord;
      ar & color;
      ar & v;
      ar & nu;
      ar & mass;
      ar & P;
      ar & rho;
      ar & h;
      ar & vol;
    }else if (tag == 1){
      ar & id;
      ar & sum;
    }else if (tag == 2){
      ar & id;
      ar & v;
      ar & nu;
    }else if (tag == 3){
      ar & coord;
      ar & color;
      ar & id;
      ar & type;
      ar & p_mass;
      ar & scale;
      ar & h;
      ar & level;     
      ar & rho;
      ar & vol;
      ar & mass;
      ar & P;
      ar & safe_guard;
    }else if (tag == 4){
    }
  }
#endif
public:

  int      type;
  int      local_id;
  Real     sum;
  Real     phi;
  Real     curv;
  Real     safe_guard;
  Real     h_nominal;
  Real     nu;
  my_real  norm;
  my_real  coord_corr;
// variables
#if PERI_DIM != 0 || SYM_DIM != 0
  concurrent_vector <p_Particle> copy;
#endif
  concurrent_vector <p_Particle> neighbor; // pointer to all the neighbors

  Particle_mesh_generation();

// functions
  void    Initialize                             (SOLVER *sph, int id_, int type_, int rank_);
  void    Clearup                                ();
  void    Reset_neighbor_info                    ();
  void    Add_neighbor                           (Particle *current_particle);
  void    Refresh_neighbor_info                  (SPH *sph, int flag);
  void    Set_timestep                           ();
  void    Set_timestep                           (int stage, SOLVER *sph);
  void    Set_local_index                        (int index_);
  void    Reset_particle_force                   (SOLVER *sph, int stage);
  void    Calculate_particle_infor               (SOLVER *sph);
  void    Calculate_particle_infor               (Levelset *level_set);
  void    Set_particle_info                      (SOLVER *sph);
  void    Get_dispersed_position                 (Real    dx);
  void    Find_and_map_particle                  (SOLVER *sph);
  void    Find_and_map_surface_particle          (SOLVER *sph);
  void    Find_and_map_segment_particle          (SOLVER *sph);
  void    Get_extended_char_cell_tag_and_idx     (SOLVER *sph, int &tag_interface, int &tag_characteristic, int &idx_characteristic);
  void    Get_level_set_char_cell_tag_and_idx    (SOLVER *sph, int &tag_interface, int &tag_characteristic, int &idx_characteristic);
  void    Get_current_particle_level_set_cell_tag(SOLVER *sph, int &tag_interface, int &tag_characteristic);
  void    Find_characteristic_particle           (SOLVER *sph);
  void    Map_characteristic_particle            (SOLVER *sph);
  void    Pop_surface_particle                   (SOLVER *sph);
  void    Get_kernel_summation                   (SOLVER *sph);
  void    Accumulate_particle_force_real         (SOLVER *sph, int stage);
  void    Accumulate_particle_force_surf         (SOLVER *sph, int stage);
  void    Accumulate_particle_force_seg          (SOLVER *sph, int stage);
  void    Map_surface_particle_force             (SOLVER *sph);
  void    Update_velocity_half_step              (SOLVER *sph);
  void    Update_coord_full_step                 (SOLVER *sph);
  void    Local_voronoi_diagram                  (SOLVER *sph);
  bool    Reconstruce_dummy_particles            (SOLVER *sph, p_Particle dummy_particle);
#if PERI_DIM != 0
  void    Set_bc_particle_info                   (p_Particle copy_particle, my_real shift);
  void    Refresh_pbc_particle_info              (p_Particle copy_particle, int flag);
#endif
#if SYM_DIM != 0
  void    Set_bc_particle_info                   (p_Particle copy_particle, my_real shift, my_real pos);
  void    Refresh_sbc_particle_info              (p_Particle copy_particle, int flag, SPH *sph);
#endif
#ifdef _MPI_
  void    Set_particle_info                      (p_Particle income_particle);
  void    Set_buffer_particle_info               (p_Particle temp, int flag);
#endif
};

#endif