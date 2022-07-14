#ifndef _LEVEL_SET_
#define _LEVEL_SET_
#include "glbcls.h"
#include "Mypool.h"
#include "lset_level_infor.h"
#include "lset_cell.h"
#include "lset_package.h"
#include "lset_binary_writer.h"
#ifdef _MPI_
#include "voronoi_lset.h"
#endif

using namespace std;
using namespace boost::mpi;
using namespace boost;

class   Levelset{

private:

public:
  // variables
  int num_level;
  int glbl_total_num_particle;
  int glbl_total_num_pkg;
  int glbl_num_interface_cell;
  int Lmin;
  int Lmax;
  int n_out;
  
  my_int ini_num_cell;
  
  my_real domain;              //domain size in each direction
  my_real box_l;               //coordinate of the left corner of the domain
  my_real box_r;
  
  int num_singularity;
  my_real *singularity;

  int num_segment;
  std::pair <my_real, my_real> *segment;

  Real convergence_error;
  Real total_mass;         // for global particle  
  Real total_mass_surface; // for surface particle
  Real total_mass_segment;
  Real total_volume;       // for initial smooth length and density
  Real total_area;
  Real total_length;
  Real maximum_phi;        // for global smooth length
  Real maximum_psi;        // for global smooth length
  Real maximum_dl;         // user defined
  Real middel_dl;          // user defined
  Real minimum_dl;         // user defined
  Real maximum_curv;       // for determine the smooth length on the surface
  Real minimum_curv;       // for determine the smooth length on the surface

  #ifdef _READ_SDF_
  char   sdf_file[256];
  my_int cells_to_cut;

  int      num_artifact_region;
  Real    *radi_artifact_region;
  my_real *center_artifact_region;
  #endif

  #if defined (_MPI_)
  Voronoi_lset lset_voronoi;
  Real         error_tolerance;
    #ifdef _CVP_LSET_INIT_
    int    num_color_tmp;
    #endif
  #endif
  
  // user defined values;
  Real fac_maximum_dl;
  Real fac_middel_dl;
  Real fac_minimum_dl;
  Real maximum_curv_artificial;

  p_Levelset_levelinfo *lset_level_info;
  std::vector<p_Levelset_package> lset_pkg;
  std::vector<p_Levelset_package> lset_bc_pkg;
  std::vector<p_Levelset_package> interface_pkg;

  Mypool<Levelset_package> lset_pkg_pool;
  
  Level_set_binary_writer     mesh_writer;

  Levelset();

  // functions
  void Initialize(SOLVER *sph, communicator &world);
  Real F_phi (my_real pos);
  void Calculate_target_density_field(SOLVER *sph, communicator &world);
  void Copy_phi_value(communicator &world);
  void Calculate_level_set(communicator &world);
  Real Get_scale(Real curv, Real phi);
  void Find_interface_packages(communicator &world);
  void Get_extended_cell_tags(communicator &world);
  int  Get_positive_pkg (std::vector<p_Levelset_package> &lset_positive_pkg);
  Real Get_phi_at_position(my_real coord);
  void Calculate_psi(communicator &world);
  void Extend_psi(communicator &world);
  void Reinitialize_psi(communicator &world, int &tag_convergence);
  void Get_psi_max(communicator &world);
  void Redistribute_curv(communicator &world);
  void Recalculate_max_curv(communicator &world);
  void Smoothing_curvarure_field (communicator &world);
  void Extend_curvature(communicator &world);
  void Smooth_curvature(int &tag_convergence, communicator &world);
  void Calculate_global_effective_curv(communicator &world);
  void Calculate_total_volume_mass(communicator &world);
  void Output_level_set_dat(int n, communicator &world);
  void Output_level_set_vti(int n, communicator &world);
  void Save_lset_restart_file(char* filename, SOLVER *sph, communicator &world);
  void Load_lset_restart_file(ifstream &load, SOLVER *sph, communicator &world);
  #if defined (_MPI_) && (defined (_CVP_LSET_) || defined (_CVP_LSET_INIT_))
  void CVP_for_initial_partitioning (std::vector<my_real> & vp_coords, std::vector<Real> &vp_scale, communicator &world);
  void Save_VP_coord_hmin (std::vector<my_real> & vp_coords, std::vector<Real> &vp_scale, communicator &world);
  #endif
};

#endif