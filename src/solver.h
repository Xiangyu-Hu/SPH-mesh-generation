#ifndef SOLVER_H
#define SOLVER_H
#include <cmath>
#include "glbcls.h"
#include "system.h"
#ifdef _INCPRS_
#include "sph_incprs.h"
#endif
#ifdef _CPRS_
#include "sph_cprs.h"
#endif
#ifdef _GSPH_
#include "sph_gsph.h"
#endif
#ifdef _ALE_
#include "sph_ale.h"
#endif
#ifdef _MESH_GENERATION_
#include "sph_mesh_generation.h"
#endif
#ifdef _EULARIAN_CAT_
#include "Eu_simulator.h"
#endif

using namespace std;
using namespace boost::mpi;

class Solver{
  int      iterate_num;        //current iterate number
  int      total_iterate_num;
  Real     output_dt;
  Real     timestep;           //current timestep size
  Real     run_time;           //current iterate time
  Real     integral_time;      //integral time in current inteval

  SOLVER current_solver;

public:  
  Solver(Initialization &Ini, communicator &world);
  void Load_solver(Initialization &Ini, communicator &world);
  void Set_current_time(Real Time);
  void MRSPH_Solver(Initialization &Ini, int n, communicator &world);
  void SPH_Mesh_Generation_Solver (Initialization &Ini, int n, communicator &world);
  void Eularian_Cat_Solver (Initialization &Ini, int n, communicator &world);
};

#endif
