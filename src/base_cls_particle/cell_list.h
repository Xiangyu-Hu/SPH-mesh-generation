#ifndef CELL_LIST_H
#define CELL_LIST_H
#include <cmath>
#include "glbcls.h"
#include "Mypool.h"

using namespace std;

#ifdef _MPI_
class Color_list{

public:
  // variables
  concurrent_vector <std::pair <int,int>> color_list; // record the color infor
  Color_list(){};

  // functions
  void Clear_data ();
  void Add_color (std::pair <int,int> current_color);
};
#endif

class Cell_list{
  int      level;

public:

// variables
#ifndef _SCLL_
  concurrent_vector <p_Particle> particle_list;
#else
  concurrent_vector <std::pair <int,p_Particle>> particle_list;
  std::vector <int> start;
  std::vector <int> end;
#endif
  Cell_list(){};

// functions
  void Initialize (Level_info *level_info);
  void Reset_tags ();
#ifndef _SCLL_
  void Add_particle (Particle *current_particle);
#else
  void Add_particle_and_key (std::pair <int,p_Particle> current_particle);
  void Find_cell_start_and_end (int num_cell);
  void Build_subcell_list (int total_num_subcell);
#endif
};

#endif
