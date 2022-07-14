#include "parallelfunc_mesh.h"
#include "glbfunc.h"

/*----------------------------------------------------------------------------------------------------------------------------*/
class ApplyGetGlobalScale{
  std::vector<p_Levelset_package> my_pkg;
public:
  Real maximum_phi;
  Real maximum_curv;
  Real minimum_curv;
  void operator()( const blocked_range<int>& r ) {
    std::vector<p_Levelset_package> my_pkg_ = my_pkg;

    for( int num=r.begin(); num!=r.end(); ++num) 
    {
      maximum_phi = AMAX1 (maximum_phi, my_pkg_[num]->maximum_phi);
      maximum_curv = AMAX1 (maximum_curv, my_pkg_[num]->maximum_curv);
      minimum_curv = AMIN1 (minimum_curv, my_pkg_[num]->minimum_curv);
    }
  }

  ApplyGetGlobalScale(ApplyGetGlobalScale &x, split)
  {
    my_pkg = x.my_pkg;
    maximum_phi = 0.;
    maximum_curv = 0.;
    minimum_curv = 1.e20;
  }

  void join( const ApplyGetGlobalScale& y)
  {
    maximum_phi = AMAX1 (maximum_phi, y.maximum_phi);
    maximum_curv = AMAX1 (maximum_curv, y.maximum_curv);
    minimum_curv = AMIN1 (minimum_curv, y.minimum_curv);
  }

  ApplyGetGlobalScale( std::vector<p_Levelset_package> &lset_pkg)
  {
    my_pkg = lset_pkg;
    maximum_phi = 0.;
    maximum_curv = 0.;
    minimum_curv = 1.e20;
  }
};

void Parallel_get_global_scale(int n, std::vector<p_Levelset_package> &lset_pkg, Levelset *level_set)
{
  ApplyGetGlobalScale GetGlobalScale(lset_pkg);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<int>(0, n), GetGlobalScale, ap);
  level_set->maximum_phi = GetGlobalScale.maximum_phi;
  level_set->maximum_curv = GetGlobalScale.maximum_curv;
  level_set->minimum_curv = GetGlobalScale.minimum_curv;
}

/*----------------------------------------------------------------------------------------------------------------------------*/
class ApplyGetGlobalPSIMAX{
  std::vector<p_Levelset_package> my_pkg;
public:
  Real maximum_psi;
  void operator()( const blocked_range<int>& r ) {
    std::vector<p_Levelset_package> my_pkg_ = my_pkg;

    for( int num=r.begin(); num!=r.end(); ++num) 
    {
      maximum_psi = AMAX1 (maximum_psi, my_pkg_[num]->maximum_psi);
    }
  }

  ApplyGetGlobalPSIMAX(ApplyGetGlobalPSIMAX &x, split)
  {
    my_pkg = x.my_pkg;
    maximum_psi = 0.;
  }

  void join( const ApplyGetGlobalPSIMAX& y)
  {
    maximum_psi = AMAX1 (maximum_psi, y.maximum_psi);
  }

  ApplyGetGlobalPSIMAX( std::vector<p_Levelset_package> &interface_pkg)
  {
    my_pkg = interface_pkg;
    maximum_psi = 0.;
  }
};

void Parallel_get_psi_max(int n, std::vector<p_Levelset_package> &interface_pkg, Levelset *level_set)
{
  ApplyGetGlobalPSIMAX GetGlobalPsiMax(interface_pkg);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<int>(0, n), GetGlobalPsiMax, ap);
  level_set->maximum_psi = GetGlobalPsiMax.maximum_psi;
}

/*----------------------------------------------------------------------------------------------------------------------------*/
class ApplyGetTotalVolumeMass{
  std::vector<p_Levelset_package> my_pkg;
public:
  Real total_volume;
  Real total_area;
  Real total_length;
  Real total_mass_segment;
  Real total_mass_surface;
  Real total_mass;
  void operator()( const blocked_range<int>& r ) {
    std::vector<p_Levelset_package> my_pkg_ = my_pkg;

    for( int num=r.begin(); num!=r.end(); ++num) 
    {
      total_volume       += my_pkg_[num]->total_volume;
      total_area         += my_pkg_[num]->total_area;
      total_length       += my_pkg_[num]->total_length;
      total_mass_segment += my_pkg_[num]->total_mass_segment;
      total_mass_surface += my_pkg_[num]->total_mass_surface;
      total_mass         += my_pkg_[num]->total_mass;
    }
  }

  ApplyGetTotalVolumeMass(ApplyGetTotalVolumeMass &x, split)
  {
    my_pkg = x.my_pkg;
    total_volume       = 0.;
    total_area         = 0.;
    total_length       = 0.;
    total_mass_segment = 0.;
    total_mass_surface = 0.;
    total_mass         = 0.;
  }

  void join( const ApplyGetTotalVolumeMass& y)
  {
      total_volume       += y.total_volume;
      total_area         += y.total_area;
      total_length       += y.total_length;
      total_mass_segment += y.total_mass_segment;
      total_mass_surface += y.total_mass_surface;
      total_mass         += y.total_mass;
  }

  ApplyGetTotalVolumeMass( std::vector<p_Levelset_package> &lset_pkg)
  {
    my_pkg = lset_pkg;
    total_volume       = 0.;
    total_area         = 0.;
    total_length       = 0.;
    total_mass_segment = 0.;
    total_mass_surface = 0.;
    total_mass         = 0.;
  }
};

void Parallel_get_total_volume_mass(int n, std::vector<p_Levelset_package> &lset_pkg, Levelset *level_set)
{
  ApplyGetTotalVolumeMass GetTotalVolumeMass(lset_pkg);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<int>(0, n), GetTotalVolumeMass, ap);
  level_set->total_mass         = GetTotalVolumeMass.total_mass;
  level_set->total_mass_surface = GetTotalVolumeMass.total_mass_surface;
  level_set->total_mass_segment = GetTotalVolumeMass.total_mass_segment;
  level_set->total_volume       = GetTotalVolumeMass.total_volume;
  level_set->total_area         = GetTotalVolumeMass.total_area;
  level_set->total_length       = GetTotalVolumeMass.total_length;
}

/*----------------------------------------------------------------------------------------------------------------------------*/
class ApplyGetNumInterfaceCell{
  std::vector<p_Levelset_package> my_pkg;
public:
  int n_interface_cell;
  void operator()( const blocked_range<int>& r ) {
    std::vector<p_Levelset_package> my_pkg_ = my_pkg;

    for( int num=r.begin(); num!=r.end(); ++num) 
    {
      n_interface_cell += my_pkg_[num]->n_interface_cell;
    }
  }

  ApplyGetNumInterfaceCell(ApplyGetNumInterfaceCell &x, split)
  {
    my_pkg = x.my_pkg;
    n_interface_cell = 0.;
  }

  void join( const ApplyGetNumInterfaceCell& y)
  {
    n_interface_cell += y.n_interface_cell;
  }

  ApplyGetNumInterfaceCell( std::vector<p_Levelset_package> &interface_pkg)
  {
    my_pkg = interface_pkg;
    n_interface_cell = 0.;
  }
};

void Parallel_get_num_of_interface_cell(int n, std::vector<p_Levelset_package> &interface_pkg, Levelset *level_set, int &num_interface_cell)
{
  ApplyGetNumInterfaceCell GetNumInterfaceCell(interface_pkg);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<int>(0, n), GetNumInterfaceCell, ap);
  num_interface_cell = GetNumInterfaceCell.n_interface_cell;
}