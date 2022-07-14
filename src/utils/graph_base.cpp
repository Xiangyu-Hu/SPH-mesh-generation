#include "graph_base.h"

//-------------------------------------------------------
// remove all edges
//-------------------------------------------------------
void Graphcls_base::Remove_all_edges(communicator &world)
{
  // remove all the edges
  std::pair<EdgeIterator, EdgeIterator> ei = edges(graph);

  EdgeIterator next = ei.first;
  for (EdgeIterator eit = ei.first; eit != ei.second; eit = next){
    ++next;
    EdgeDescriptor ed = *eit;
    remove_edge(ed,graph);
  }

  if(int(num_edges(graph)) != 0){
    cout<<"The removing edge procedures are wrong !!!"<<endl;
    world.abort(-1);
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Remove_all_edges finished\n";
  out.close();
#endif
}

//-------------------------------------------------------
// remove all vertexes
//-------------------------------------------------------
void Graphcls_base::Remove_all_verts(communicator &world)
{
  VertexIterator vi, vi_end, next;
  boost::tie( vi, vi_end ) = vertices (graph);
  
  for (next = vi; vi != vi_end; vi = next) {
    next++;
    remove_vertex(*vi, graph);
  }

  if(int(num_vertices(graph)) != 0){
    cout<<"The removing vertices procedures are wrong !!!"<<endl;
    world.abort(-1);
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Remove_all_verts finished\n";
  out.close();
#endif
}

//-------------------------------------------------------
// add edge
//-------------------------------------------------------
void Graphcls_base::Add_edge(int color_i, int color_j)
{
  //scoped lock in order to release lock automatically
  tbb::mutex::scoped_lock lock(graph_Mutex);
  std::pair<EdgeDescriptor, bool> test_1 = edge(color_i,color_j,graph);
  std::pair<EdgeDescriptor, bool> test_2 = edge(color_j,color_i,graph);
  if (test_1.second == false && test_2.second == false){
    add_edge(color_i,color_j,graph);
//    cout<<"Edge "<<color_i<<" & "<<color_j<<"added\n";
  }
}