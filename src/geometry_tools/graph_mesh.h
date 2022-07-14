#ifndef GRAPH_MESH_H
#define GRAPH_MESH_H
#include "glbcls.h"
#include "tetgen.h"
#include "Mypool.h"
#include "graph_base.h"

using namespace std;
using namespace boost;
using namespace boost::mpi;

inline bool pairCompare_intint(const std::pair<int, int>& firstElem, const std::pair<int, int>& secondElem) {
  return firstElem.first <= secondElem.first;
}

inline int Get_color (int index_, int *off_set, int nproc, int total)
{
  if (index_ < 0 || index_ >= total) return -1;
  
  for (int i = 0; i < nproc-1; i++){
    if (index_ >= off_set[i] && index_ < off_set[i+1])
      return i;
    else if (i+1 == nproc-1 && index_ < total)
      return nproc-1;
  }
  
  return -2;
}

class int4_graph{

public:

  int i;  
  int j;
  int k;
  int l;
  mutable Real vol;
  mutable Real aspect;
  mutable Real radius_ratio;
  // mutable Real dihedangle[4];
  mutable Real mindihedangle;

  int4_graph() { };
  int4_graph (const int4_graph& std4)
  {  
    i = std4.i;
    j = std4.j;
    k = std4.k;
    l = std4.l;
    vol = std4.vol;
    aspect = std4.aspect;
    radius_ratio = std4.radius_ratio;
    mindihedangle = std4.mindihedangle;
    // dihedangle[0] = std4.dihedangle[0];
    // dihedangle[1] = std4.dihedangle[1];
    // dihedangle[2] = std4.dihedangle[2];
    // dihedangle[3] = std4.dihedangle[3];
  };

  bool operator < (const int4_graph& std4) const 
  {  
    if(i < std4.i)
      return  i < std4.i;
    else if(i == std4.i && j < std4.j)
      return  j < std4.j;
    else if(i == std4.i && j == std4.j && k < std4.k)
      return  k < std4.k;
    else if(i == std4.i && j == std4.j && k == std4.k && l < std4.l)
      return  l < std4.l;
    else
      return false;
  };

  bool operator == (const int4_graph& std4) const 
  {  
    return  (i == std4.i) && (j == std4.j) && (k == std4.k) && (l == std4.l);  
  };

  bool operator > (const int4_graph& std4) const 
  {  

    if(i > std4.i)
      return  i > std4.i;
    else if(i == std4.i && j > std4.j)
      return  j > std4.j;
    else if(i == std4.i && j == std4.j && k > std4.k)
      return  k > std4.k;
    else if(i == std4.i && j == std4.j && k == std4.k && l > std4.l)
      return  l > std4.l;
    else
      return false;
  };
};

class int3_graph{

public:

	int i; 	
	int j;
	int k;
  mutable Real vol;
  mutable Real aspect;

	int3_graph() { };
	int3_graph (const int3_graph& std3)
	{  
		i = std3.i;
		j = std3.j;
		k = std3.k;
    vol = std3.vol;
    aspect = std3.aspect;
	};

	bool operator < (const int3_graph& std3) const 
	{  
		if(i < std3.i)
			return  i < std3.i;
		else if(i == std3.i && j < std3.j)
			return  j < std3.j;
		else if(i == std3.i && j == std3.j && k < std3.k)
			return  k < std3.k;
		else
			return false;
	};

	bool operator == (const int3_graph& std3) const 
	{  
		return  (i == std3.i) && (j == std3.j) && (k == std3.k);  
	};

	bool operator > (const int3_graph& std3) const 
	{  

		if(i > std3.i)
			return  i > std3.i;
		else if(i == std3.i && j > std3.j)
			return  j > std3.j;
		else if(i == std3.i && j == std3.j && k > std3.k)
			return  k > std3.k;
		else
			return false;
	};
};

class int2_graph{

public:

	int i; 	
	int j;

	int2_graph() { };
	int2_graph (const int2_graph& std2)
	{  
		i = std2.i;
		j = std2.j;
	};

	bool operator < (const int2_graph& std2) const 
	{  
		if(i < std2.i)
			return  i < std2.i;
		else if(i == std2.i && j < std2.j)
			return  j < std2.j;
		else
			return false;
	};

	bool operator == (const int2_graph& std2) const 
	{  
		return  (i == std2.i) && (j == std2.j);  
	};

	bool operator > (const int2_graph& std2) const 
	{  
		if(i > std2.i)
			return  i > std2.i;
		else if(i == std2.i && j > std2.j)
			return  j > std2.j;
		else
			return false;
	};
};

class int1_graph{

public:

  int i;

  int1_graph() { };
  int1_graph (const int1_graph& std1)
  {  
    i = std1.i;
  };

  bool operator < (const int1_graph& std1) const 
  {  
    if(i < std1.i)
      return  i < std1.i;
    else
      return false;
  };

  bool operator == (const int1_graph& std1) const 
  {  
    return  i == std1.i;  
  };

  bool operator > (const int1_graph& std1) const 
  {  
    if(i > std1.i)
      return  i > std1.i;
    else
      return false;
  };
};

class Graphmeshcls : public Graphcls_base{

  public:
  // variables
  // quality section
  // Tris
  Real maximum_angle;
  Real minimum_angle;
  Real maximum_tri_aspect_ratio;
  Real maximum_edge;
  Real minimum_edge;
  Real maximum_area;
  Real minimum_area;
  Real G_min;
  Real G_avg;
  Real minimum_angle_avg;
  int  ntri_sm_10;
  int  ntri_sm_20;
  int  ntri_sm_30;
  int  ntri_sm_40;
  int  glbl_num_tris;

  // Tets
  Real maximum_dihedangle;
  Real minimum_dihedangle;
  Real minimum_dihedangle_avg;
  Real minimum_tet_aspect_ratio;
  Real maximum_tet_aspect_ratio;
  Real average_tet_aspect_ratio;
  Real maximum_tetradius;
  Real minimum_tetradius;
  Real average_tetradius;
  Real maximum_volume;
  Real minimum_volume;
  Real average_volume;
  int  ntet_sm_10;
  int  ntet_sm_20;
  int  ntet_sm_30;
  int  ntet_sm_40;
  int  glbl_num_tets;

  int  glbl_communication_tets;
  int  glbl_communication_tris;
  Real glbl_communication_volume;
  
  std::set<int3_graph> tris;

  std::set<int4_graph> tets;
  std::set<int4_graph> tets_error;

  Graphmeshcls(){};
  // functions
  void    Initialize           (int nvert, communicator &world);
  #ifdef _MESH_GENERATION_
  void    Reconstruct_tris                        (std::vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  void    Tri_quality_statistics                  (std::vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  void    Reconstruct_tets                        (std::vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  void    Tet_quality_statistics                  (std::vector<p_Particle> particle_total, SOLVER *sph, communicator &world);
  bool    Check_whether_this_tet_should_be_delete (int stage, p_Particle particle_1, p_Particle particle_2, p_Particle particle_3, p_Particle particle_4, int iRank, SOLVER *sph);
  void    Write_tri_quality_plt                   (int n, Real run_time, communicator &world);
  void    Write_tet_quality_plt                   (int n, Real run_time, communicator &world);
  #endif
};

#endif