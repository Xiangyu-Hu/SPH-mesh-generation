#include <cmath>
#include "glbcls.h"
#include "system.h"
#include "solver.h"

using namespace std;
using namespace tbb;
using namespace boost;
using namespace boost::mpi;
/***************************************************/
/*                                                 */
/*                Program Main Prototype           */
/*                  Writing by Ji, Zhe             */
/*                  Date: Jan 13, 2016             */
/*                                                 */
/***************************************************/

int main(int argc, char **argv)
{
  // mpi initialization
  environment env(argc,argv);
  communicator world;

  //system parameters
  Initialization Ini(world);

  //Initialize TBB
#ifdef _MPI_
  #ifdef _CVP_LSET_INIT_
    task_scheduler_init init();
  #else
    task_scheduler_init init(4);
  #endif
#else
  task_scheduler_init init();
#endif
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::trunc);
  out<<"<<<<<PID "<<world.rank()<<" started!\n";
  out.close();
#endif

  //set solver
  Solver solver(Ini, world);

  solver.Load_solver (Ini, world);
  
  return 0;
}
