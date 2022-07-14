#ifndef VORONOI_CONSTRUCTION_H
#define VORONOI_CONSTRUCTION_H
#include "glbcls.h"
#include "glbfunc.h"
#include "Mypool.h"
#include "voronoi_base.h"
#ifdef _MESH_GENERATION_
#include "particle_mesh_generation.h"
#include "graph_mesh.h"
#endif

using namespace std;
using namespace boost;

// graph based triangulation
class Voronoi_construction : public Voronoi_base{

  public:
  // variables

  Voronoi_construction    (){};
  // functions
  void Initialize         (){};
  #ifdef _MESH_GENERATION_
  void Local_VD_generation(Particle *sp, Graphmeshcls *mesh);
  #endif
};

#endif