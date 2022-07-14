#ifndef _MESH_BINARY_WRITER_H_
#define _MESH_BINARY_WRITER_H_
#include <sstream>
#include "glbcls.h"
#include "writer.h"
#include "Mypool.h"

using namespace std;
using namespace boost;
using namespace boost::mpi;

class Mesh_binary_writer : public Writer{

  public:
  // variables
  std::map<std::string, bool> print_info;

  Mesh_binary_writer      ();
  // functions

  void    Output_pvtu_file_tris (int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  void    Output_pvd_file_tris  (int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  void    Output_pvtu_file_tets (int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  void    Output_pvd_file_tets   (int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  void    Output_vtu_file_tris  (int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  void    Output_vtu_file_tets  (int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  void    Output_mesh           (int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
};

#endif