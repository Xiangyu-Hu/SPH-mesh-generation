#ifndef _LEVEL_SET_LEVEL_INFO_
#define _LEVEL_SET_LEVEL_INFO_

#include "glbcls.h"
#include "glbfunc.h"
#include "lset_package.h"
#include "lset_cell.h"
#include "Mypool.h"
#include "parallelfunc_mesh.h"
#include "boost/multi_array.hpp"

using namespace std;
using namespace boost::mpi;
using namespace boost;
using namespace tbb;

typedef boost::multi_array<Real, 3> array_real;
typedef boost::multi_array<DTAG, 3> array_tag;
typedef boost::multi_array<p_Levelset_package, 3> array_lset_pkg;

class   Levelset_levelinfo{

private:

public:
  //variables
  int level;
  int Lmin;
  int Lmax;
  int glbl_total_num_pkg;

  Real scale;
  Real dl;

  my_int num_pkg;
  my_int glbl_num_pkg;
  my_int glbl_pkg_start;
  my_int glbl_pkg_end;
  my_int pkg_start;
  my_int pkg_end;
  
  my_int             num_pkg_margin;     //number of packages in periodical boundary buffer

  my_real dpkg;
  my_real dcell;
  my_real domain;
  my_real box_l;
  my_real box_r;

  array_lset_pkg table_lset_pkg_list;

  Mypool<Levelset_package> *p_lset_pkg_pool;

  Levelset_levelinfo(){};
  
  // functions
  void Initialize (int i, Levelset *level_set, communicator &world);
  void Set_cell_topology (Levelset *level_set, communicator &world);
  void Find_interface_packages (Levelset *level_set, communicator &world);
  void Get_extended_cell_tags (Levelset *level_set, communicator &world);
  void Allocate_memory();
  void Initial_pkg_list_and_reset_tags(Levelset *level_set);
  void Define_levelset_full_field(Levelset *level_set);
  void Get_id_pkg_cell (my_real pos, my_int &id_pkg, my_int &id_cell);
  void Search_for_char_cell_within_dx (my_real coord, Real rad, int &tag_interface, int &tag_characteristic, int &idx_characteristic);
  int  Get_unique_cell_id (my_int id_pkg, my_int id_cell);
  void Output_level_set(Levelset *level_set, communicator &world, char *filename, int n);
  void Write_level_set_rstfile(Levelset *level_set, communicator &world, char *filename);
  void Load_level_set_rstfile(ifstream &load, Levelset *level_set, communicator &world);
  #ifdef _READ_SDF_
  void Load_SDF_file(Levelset *level_set, communicator &world);
  #endif
};

#endif