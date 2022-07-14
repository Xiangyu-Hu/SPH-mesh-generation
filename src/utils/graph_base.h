#ifndef GRAPH_BASE_H
#define GRAPH_BASE_H
#include "glbcls.h"
#include "Mypool.h"

using namespace std;
using namespace boost;
using namespace boost::mpi;

class Graphcls_base {

  public:
  // variables
  Graph          graph;
  tbb::mutex     graph_Mutex;
  int            total_edges;

  Graphcls_base(){};
  // functions

  void     Remove_all_edges(communicator &world);
  void     Remove_all_verts(communicator &world);
  void     Add_edge(int color_i, int color_j);
};

#endif
