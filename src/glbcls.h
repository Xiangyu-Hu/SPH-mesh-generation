#ifndef GLBCLS_H
#define GLBCLS_H

//-----------------------------------------
//  Self-defined data types
//----------------------------------------
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <string.h>
#include <fstream>
#include <vector>
#include "glbparam.h"
#include "tbb/tick_count.h"
#include "tbb/atomic.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/blocked_range.h"
#include "tbb/blocked_range2d.h"
#include "tbb/blocked_range3d.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/scalable_allocator.h"
#include <boost/mpi.hpp>
#include <boost/config.hpp>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
//#ifdef _MPI_  // JZ09092018::commented for mesh generation
#include <boost/serialization/utility.hpp>
#include <boost/serialization/string.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/edge_coloring.hpp>
#include <boost/graph/properties.hpp>
//#endif

using namespace std;
using namespace tbb;
using namespace boost;
using namespace boost::serialization;

// JZ09092018::added for mesh generation
class   Graphcls_base;
class   Voronoi_base;
class   Vo_particle;
class   Writer;

class   Particle_base;
class   Cell_list;
class   Level_info;
class   SPH;
class   Solver;
class   Initialization;
class   Visualize;
typedef Cell_list* p_Cell_list;
typedef Level_info* p_Level_info;
typedef char DTAG;  //define DTAG as data type to tag cells;

#ifdef _INCPRS_
class   SPH_incprs;
class   Particle_incprs;
typedef Particle_incprs* p_Particle;
#endif
#ifdef _CPRS_
class   SPH_cprs;
class   Particle_cprs;
typedef Particle_cprs* p_Particle;
#endif
#ifdef _GSPH_
class   SPH_gsph;
class   Particle_gsph;
typedef Particle_gsph* p_Particle;
#endif
#ifdef _ALE_
class   SPH_ale;
class   Particle_ale;
typedef Particle_ale* p_Particle;
#endif
#ifdef _MESH_GENERATION_
class   GeometryCalculator;
class   Graphmeshcls;
class   Voronoi_construction;
class   Mesh_binary_writer;
class   Particle_binary_writer;
class   SPH_mesh_generation;
class   Particle_mesh_generation;
typedef Particle_mesh_generation* p_Particle;
class   Levelset;
class   Levelset_levelinfo;
class   Levelset_package;
class   Levelset_cell;
class   Level_set_binary_writer;
typedef Levelset_cell* p_Levelset_cell;
typedef Levelset_package* p_Levelset_package;
typedef Levelset_levelinfo* p_Levelset_levelinfo;
#endif
#ifdef _EULARIAN_CAT_
class   Eularian_cat;
class   Eularian_levelinfo;
class   Eularian_package;
class   Eularian_cell;
class   Fluidcomprs_solver;
class   Fluidcomprs_STsolver;
class   Stencil;
class   Boundary_Condition;
class   Boundary_Condition_cat;
class   Material;
typedef Eularian_cell* p_Eularian_cell;
typedef Eularian_package* p_Eularian_package;
typedef Eularian_levelinfo* p_Eularian_levelinfo;
#endif

typedef adjacency_list<listS, vecS, undirectedS, size_t, size_t> Graph;
typedef adjacency_list<listS, vecS, bidirectionalS, size_t, size_t> Graph_Vo;
typedef graph_traits<Graph>::vertex_descriptor  VertexDescriptor;
typedef graph_traits<Graph>::edge_descriptor   EdgeDescriptor;
typedef graph_traits<Graph>::vertex_iterator  VertexIterator;
typedef graph_traits<Graph>::edge_iterator  EdgeIterator;
typedef graph_traits<Graph>::in_edge_iterator  InEdgeIterator;
typedef graph_traits<Graph>::out_edge_iterator  OutEdgeIterator;
typedef graph_traits<Graph>::adjacency_iterator  AdjacencyIterator;

#ifdef _MPI_ // JZ09092018::modified for mesh generation
class   Graphcls;
class   Voronoi;
class   Color_list;
typedef Color_list*  p_Color_list;
typedef Vo_particle* p_Vo_particle;
  #ifdef _MESH_GENERATION_
  class Voronoi_lset;
  #endif
#endif

#ifdef _DOUBLE_
typedef double Real;
#endif

#ifdef _SINGLE_
typedef float  Real;
#endif

class my_int{
  private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & i;
    ar & j;
    ar & k;
  }
  public:
  int i; int j; int k;
  my_int(){i = j = k = 0;};
};

class my_real{
  private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & i;
    ar & j;
    ar & k;
  }
  public:
  Real i; Real j; Real k;
  my_real(){i = j = k = 0.;};
};

class my_tensor{
  private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & t00;
    ar & t01;
    ar & t02;
    ar & t10;
    ar & t11;
    ar & t12;
    ar & t20;
    ar & t21;
    ar & t22;
  }
  public:
    Real t00; Real t01; Real t02; Real t10; Real t11; Real t12; Real t20; Real t21; Real t22;
  my_tensor(){t00 = t01 = t02 = t10 = t11 = t12 = t20 = t21 = t22 = 0.;};
};

class my_sym_tensor{
  private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & t00;
    ar & t01;
    ar & t02;
    ar & t11;
    ar & t12;
    ar & t22;
  }
  public:
    Real t00; Real t01; Real t02; Real t11; Real t12; Real t22;
  my_sym_tensor(){t00 = t01 = t02 = t11 = t12 = t22 = 0.;};
};

#ifdef _MPI_
typedef struct{
  int index;
  int color;
}VertexProperty;

typedef struct{
  int index;
  int color;
}EdgeProperty;

//serialization 3d martrix
template<class T> class serialization_3d_martrix{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    for( int i=0; i<Res.i; i++)
      for( int j=0; j<Res.j; j++)
        for( int k=0; k<Res.k; k++)
        {
          ar & Martix[i][j][k];
        }
    }
public:

  T*** Martix;
  my_int Res;
  serialization_3d_martrix(T*** &matrix, my_int res);
};
//serialization 3d martrix
template<class T> class serialization_3d_martrix_v1{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    for( int i=0; i<Res.i; i++)
      for( int j=0; j<Res.j; j++)
        for( int k=0; k<Res.k; k++)
        {
          ar & Martix[i][j][k];
        }
    }
public:

  T*** Martix;
  my_int Res;
  
  serialization_3d_martrix_v1() { };
  
  void memory_allocate(my_int res);
  void memory_release();
};
//serialization vector
template<class T> class serialization_vector{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    ar & mem_size;
    for(int i=0; i<mem_size; i++)
    {
      if(tag == 1)
        ar & Vector[i];
      else
      {
        T temp;
        ar & temp;
        Vector.push_back(temp);
      }
    }
    }
public:

  int tag;
  int mem_size;
  concurrent_vector <T> Vector;

  serialization_vector() {tag = mem_size = 0;};
};

template<class T> class serialization_vector_pointer{
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    ar & Vector;
    }
public:
  std:: vector <T> Vector;
};
#endif
#endif
