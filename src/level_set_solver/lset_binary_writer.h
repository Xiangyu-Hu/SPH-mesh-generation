#ifndef _LSET_BINARY_WRITER_H_
#define _LSET_BINARY_WRITER_H_
#include <sstream>
#include "glbcls.h"
#include "glbfunc.h"
#include "writer.h"
#include "Mypool.h"

using namespace std;
using namespace boost;
using namespace boost::mpi;

class Level_set_binary_writer : public Writer{

  public:
  // variables
  std::map<std::string, bool> print_info;

  Level_set_binary_writer      ();
  // functions

  void    Output_pvti_file_lset (int n_out, p_Levelset_levelinfo lset_level_info, communicator &world);
  void    Output_pvd_file_lset  (int n_out, p_Levelset_levelinfo lset_level_info, communicator &world);
  void    Output_vti_file_lset  (int n_out, p_Levelset_levelinfo lset_level_info, communicator &world);
  void    Output_lset           (int n_out, p_Levelset_levelinfo lset_level_info, communicator &world);
};

#endif