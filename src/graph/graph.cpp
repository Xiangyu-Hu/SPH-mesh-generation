#include "boost/multi_array.hpp"
#include "level_infor.h"
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
#include "glbfunc.h"
#include"graph.h"

/***************************************************/
/*                                                 */
/*     Functions defined in class "Graphls"        */
/*                                                 */
/***************************************************/

//-------------------------------------------------------
// Initialization
//-------------------------------------------------------
void Graphcls::Initialize(communicator &world)
{
  time_for_exchange_graph_topology  = new Real[world.size()+1];
  time_for_construct_graph          = new Real[world.size()+1];
  time_for_edge_coloring            = new Real[world.size()+1];
  for (int i = 0; i < world.size()+1; i++){
    time_for_exchange_graph_topology[i]  = 0.;
    time_for_construct_graph[i]          = 0.;
    time_for_edge_coloring[i]            = 0.;
  }
  if(world.rank()==0){
    for(int i=0; i<world.size(); i++)
      add_vertex(i,graph);
  }

  // initialize graph matrix
  graph_matrix = new tbb::atomic<int>*[world.size()];
  for (int i = 0; i < world.size(); i++){
    graph_matrix[i] = new tbb::atomic<int>[world.size()];
    for (int j = 0; j < world.size(); j++){
        graph_matrix[i][j] = 0;
    }
  }

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<class Graphcls is initialized\n";
  out.close();
#endif
}
//-------------------------------------------------------
// Update graph modified version
//-------------------------------------------------------
void Graphcls::Update_graph_SPH_modified(SPH *sph, communicator &world)
{
  num_level    = sph->num_level;
  level_info   = sph->level_info;

  serialization_vector <unsigned int> exchange_vector;
  exchange_vector.Vector.clear();

    time_exchange_graph_topology.restart();
  if(world.rank()==0){

    Remove_all_edges(world);

    serialization_vector <unsigned int> *gather_exchange_vector;
    gather_exchange_vector = new serialization_vector <unsigned int>[world.size()];
#ifdef _NARROW_BAND_GRAPH_
    for( int i=0; i< num_level; i++)
      level_info[i]->Find_narrow_band_cells();
#endif
    for( int i=0; i< num_level; i++)
      level_info[i]->Initialize_color_infor();

    for( int ilevel=0; ilevel < num_level; ilevel++){
      static affinity_partitioner ap;
      parallel_for( blocked_range<int>(0, world.size()),
               [&](const blocked_range<int>& r){
        for(int i=r.begin(); i!=r.end(); ++i){
          gather_exchange_vector[i].Vector.clear();
          gather_exchange_vector[i].mem_size = 0;
          gather_exchange_vector[i].tag = 0;
        }
      }, ap);

      gather(world,exchange_vector,gather_exchange_vector,0);

      for( int i=1; i<world.size(); i++){
        my_int st;
        my_set_const(st, 0);
        st.i = int(gather_exchange_vector[i].Vector[0]);
        st.j = int(gather_exchange_vector[i].Vector[1]);
        st.k = int(gather_exchange_vector[i].Vector[2]);

        static affinity_partitioner ap;
        parallel_for( blocked_range<int>(3, int(gather_exchange_vector[i].Vector.size())),
                 [&](const blocked_range<int>& r){
          for(int j=r.begin(); j!=r.end(); ++j){
            int i_t, j_t, k_t, level_temp;
            int i_tt, j_tt, k_tt; i_tt = j_tt = k_tt=0;
            unsigned int data;
            
            data = gather_exchange_vector[i].Vector[j];
            if (data == 0xffffffff) {
              cout << "compress error" << endl;world.abort(-1);
            }
            // uncompress
            uncompress(data, i_t, j_t, k_t, level_temp);

            std:: pair <int,int> pair_t;
            pair_t.first = i;
            pair_t.second = ilevel;

            level_info[ilevel]->Shift_cell_index(i_tt, j_tt, k_tt, i_t+st.i, j_t+st.j, k_t+st.k);
            level_info[ilevel]->color_list[i_tt][j_tt][k_tt]->Add_color(pair_t);
          }
        }, ap);
      }
    }
#if PERI_DIM != 0
    for( int i = 0; i< num_level; i++)
      level_info[i]->Update_PBC_color_info();
#endif

    for( int i = num_level-2; i >= 0; i--)
      level_info[i]->Update_every_level_color_info(level_info, i+1);

    delete [] gather_exchange_vector;
  }else{
#ifdef _NARROW_BAND_GRAPH_
    for( int i=0; i< num_level; i++)
      level_info[i]->Find_narrow_band_cells();
#endif
    for( int i=0; i< num_level; i++){
      exchange_vector.Vector.clear();
      exchange_vector.mem_size = 0;
      exchange_vector.tag = 0;

      level_info[i]->Exchange_graph_topology(world,exchange_vector);

      exchange_vector.mem_size = int(exchange_vector.Vector.size());
      exchange_vector.tag = 1;

      gather(world,exchange_vector,0);
    }
  }
    time_for_exchange_graph_topology[world.rank()] += time_exchange_graph_topology.elapsed();
    time_construct_graph.restart();

  if (world.rank()==0){
    // initialize graph matrix
    for (int i = 0; i < world.size(); i++){
      for (int j = 0; j < world.size(); j++){
          graph_matrix[i][j] = 0;
      }
    }

    // traverse the tree and establish global topology graph
    for( int i=0; i < num_level; i++){
//      level_info[i]->Traverse_tree_and_construct_graph(this);
      level_info[i]->Traverse_tree_and_construct_graph_modified(this);
    }

    for (int i = 0; i < world.size(); i++){
      for (int j = 0; j < world.size(); j++){
        if (graph_matrix[i][j] == 1){
          Add_edge(i, j);
        }
      }
    }

    total_edges = int(num_edges(graph));
  }

    time_for_construct_graph[world.rank()] += time_construct_graph.elapsed();
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_graph_SPH finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// Update graph
//-------------------------------------------------------
void Graphcls::Update_graph_SPH(SPH *sph, communicator &world)
{
  num_level    = sph->num_level;
  level_info   = sph->level_info;

  serialization_vector <unsigned int> exchange_vector;
  exchange_vector.Vector.clear();

    time_exchange_graph_topology.restart();
  if(world.rank()==0){

    Remove_all_edges(world);

    serialization_vector <unsigned int> *gather_exchange_vector;
    gather_exchange_vector = new serialization_vector <unsigned int>[world.size()];

    for( int i=0; i< num_level; i++)
      level_info[i]->Initialize_color_infor();

    for( int ilevel=0; ilevel < num_level; ilevel++){
      static affinity_partitioner ap;
      parallel_for( blocked_range<int>(0, world.size()),
               [&](const blocked_range<int>& r){
        for(int i=r.begin(); i!=r.end(); ++i){
          gather_exchange_vector[i].Vector.clear();
          gather_exchange_vector[i].mem_size = 0;
          gather_exchange_vector[i].tag = 0;
        }
      }, ap);
      gather(world,exchange_vector,gather_exchange_vector,0);
      
      for( int i=1; i<world.size(); i++){
        my_int st;
        my_set_const(st, 0);
        st.i = int(gather_exchange_vector[i].Vector[0]);
        st.j = int(gather_exchange_vector[i].Vector[1]);
        st.k = int(gather_exchange_vector[i].Vector[2]);

        static affinity_partitioner ap;
        parallel_for( blocked_range<int>(3, int(gather_exchange_vector[i].Vector.size())),
                 [&](const blocked_range<int>& r){
          for(int j=r.begin(); j!=r.end(); ++j){
            int i_t, j_t, k_t, level_temp;
            int i_tt, j_tt, k_tt; i_tt = j_tt = k_tt=0;
            unsigned int data;
            
            data = gather_exchange_vector[i].Vector[j];
            if (data == 0xffffffff) {
              cout << "compress error" << endl;world.abort(-1);
            }
            // uncompress
            uncompress(data, i_t, j_t, k_t, level_temp);

            std:: pair <int,int> pair_t;
            pair_t.first = i;
            pair_t.second = ilevel;

            level_info[ilevel]->Shift_cell_index(i_tt, j_tt, k_tt, i_t+st.i, j_t+st.j, k_t+st.k);
            level_info[ilevel]->color_list[i_tt][j_tt][k_tt]->Add_color(pair_t);
          }
        }, ap);
      }
    }
#if PERI_DIM != 0
    for( int i = 0; i< num_level; i++)
      level_info[i]->Update_PBC_color_info();
#endif
    for( int i = num_level-2; i >= 0; i--)
      level_info[i]->Update_every_level_color_info(level_info, i+1);
    delete [] gather_exchange_vector;
  }else{
  
    for( int i=0; i< num_level; i++){
      exchange_vector.Vector.clear();
      exchange_vector.mem_size = 0;
      exchange_vector.tag = 0;

      level_info[i]->Exchange_graph_topology(world,exchange_vector);

      exchange_vector.mem_size = int(exchange_vector.Vector.size());
      exchange_vector.tag = 1;

      gather(world,exchange_vector,0);
    }
  }
    time_for_exchange_graph_topology[world.rank()] += time_exchange_graph_topology.elapsed();
    time_construct_graph.restart();

  if (world.rank()==0){
    // traverse the tree and establish global topology graph
    for( int i=0; i < num_level; i++)
      level_info[i]->Traverse_tree_and_construct_graph(this);

    total_edges = int(num_edges(graph));
  }

    time_for_construct_graph[world.rank()] += time_construct_graph.elapsed();
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_graph_SPH finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// Edge coloring
//-------------------------------------------------------
void Graphcls::Handle_graph_edge_coloring(communicator &world)
{
    time_edge_coloring.restart();

  if(world.rank()==0){
    maximal_degree = 0;

    edge_pool_colored.clear();
    vector <pair<int,int>>(edge_pool_colored).swap(edge_pool_colored);
    edge_color.clear();
    vector <int>(edge_color).swap(edge_color);

    std::pair<VertexIterator, VertexIterator> vi = vertices(graph);
    ParallelGetDegree(vi.first,vi.second,&graph,maximal_degree);

    // perform edge coloring
    size_t colors = edge_coloring(graph, get(edge_bundle, graph));

    // record the coloring number
    graph_iteration = int(colors);

    if(graph_iteration < maximal_degree && graph_iteration > (maximal_degree + 1)){
      cout<<"The edge coloring number is wrong !!!"<<endl;
      cout<<"Maximal degree :"<<" "<<maximal_degree<<" "<<"color number :"<<" "<<int(colors)<<" "<<"edge number :"<<" "<<num_edges(graph)<<endl;
    }
    // broadcast the coloring number
    broadcast(world,graph_iteration,0);

    std::pair<EdgeIterator, EdgeIterator> ei = edges(graph);

    //store the edge index and edge color
    for (EdgeIterator eit = ei.first; eit != ei.second; ++eit){
      EdgeDescriptor ed = *eit;

      // store the current edge infor in pool
      std::pair<int,int> pair_index;
      pair_index.first = int(graph[source(ed,graph)]);
      pair_index.second = int(graph[target(ed,graph)]);
      edge_pool_colored.push_back(pair_index);
      edge_color.push_back(graph[ed]);
    }
  }else{
    broadcast(world,graph_iteration,0);
  }

    time_for_edge_coloring[world.rank()] += time_edge_coloring.elapsed();
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Handle_graph_edge_coloring finished\n";
  out.close();
#endif
}
