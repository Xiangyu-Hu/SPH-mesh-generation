#ifndef SYSTEM_H
#define SYSTEM_H
#include "glbcls.h"

using namespace std;
using namespace boost::mpi;

// system parameters for global control
class Initialization{

public:
  Real     start_time, end_time;
  int      n_thread;
  Real     output_dt;
  int      num_output;           //num_output_file
  Initialization(communicator &world);
  Real Get_start_time();
  Real Get_dt();
  int  Get_num_timestep();
};

#endif
