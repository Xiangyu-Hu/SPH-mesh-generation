#ifndef GRAPH_H
#define GRAPH_H
#include "glbcls.h"
#include "Mypool.h"
#include "graph_base.h"

using namespace std;
using namespace boost;
using namespace boost::mpi;

class Graphcls : public Graphcls_base{

  public:
  // variables
  int            num_level;
  int            target_processor;
  int            graph_iteration;
  int            maximal_degree;
  Real           *time_for_exchange_graph_topology;
  Real           *time_for_construct_graph;
  Real           *time_for_edge_coloring;
  timer          time_exchange_graph_topology;
  timer          time_construct_graph;
  timer          time_edge_coloring;
  p_Level_info   *level_info;
  // colored edge pool
  vector <pair<int,int>> edge_pool_colored;
  // color of each edge
  vector <int>   edge_color;
  tbb::atomic<int>   **graph_matrix;

  Graphcls(){};
  // functions
  void     Initialize(communicator &world);
  void     Update_graph_SPH(SPH *sph, communicator &world);
  void     Update_graph_SPH_modified(SPH *sph, communicator &world);
  void     Handle_graph_edge_coloring(communicator &world);

};

#endif
