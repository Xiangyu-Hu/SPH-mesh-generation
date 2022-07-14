#ifndef _PARTICLE_BINARY_WRITER_H_
#define _PARTICLE_BINARY_WRITER_H_
#include <sstream>
#include "glbcls.h"
#include "writer.h"
#include "Mypool.h"

using namespace std;
using namespace boost;
using namespace boost::mpi;

class Particle_binary_writer : public Writer{

  public:
  // variables
  std::map<std::string, bool> print_info;

  Particle_binary_writer      ();
  // functions

  void    Output_pvtu_file(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  void    Output_pvd_file (int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  void    Output_vtu_file (int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  void    Output_particle (int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
};

#endif