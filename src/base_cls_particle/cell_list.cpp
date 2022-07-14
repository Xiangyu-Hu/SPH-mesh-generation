#include "parallelfunc.h"
#ifdef _INCPRS_
#include "particle_incprs.h"
#endif
#ifdef _CPRS_
#include "particle_cprs.h"
#endif
#ifdef _GSPH_
#include "particle_gsph.h"
#endif
#ifdef _ALE_
#include "particle_ale.h"
#endif
#ifdef _MESH_GENERATION_
#include "particle_mesh_generation.h"
#endif
#include "level_infor.h"
#include "cell_list.h"

#ifdef _MPI_
/***************************************************/
/*                                                 */
/*    Functions defined in class "Color_list"      */
/*                                                 */
/***************************************************/

//-------------------------------------------------------
// reset the counter to 0
//-------------------------------------------------------
void Color_list::Clear_data()
{
  color_list.clear();
//  concurrent_vector <std::pair <int,int>>(color_list).swap(color_list);
//  color_list.shrink_to_fit();
//  if (int(color_list.capacity()) != 0){
//    cout<<"color list size is wrong!!!\n";
//  }
}

//-------------------------------------------------------
// add color pair to color list
//-------------------------------------------------------
void Color_list::Add_color (std::pair <int,int> current_color)
{
  color_list.push_back(current_color);
}
#endif

/***************************************************/
/*                                                 */
/*     Functions defined in class "Cell_list"      */
/*                                                 */
/***************************************************/

//-------------------------------------------------------
// initialze all the necessary parameters for cell_list
//-------------------------------------------------------
void Cell_list::Initialize
(Level_info *level_info)
{
  level = level_info->level;
  particle_list.clear();
#ifdef _SCLL_
  start.clear();
  end.clear();
#endif
}
//-------------------------------------------------------
// reset the counter to 0
//-------------------------------------------------------
void Cell_list::Reset_tags()
{
  particle_list.clear();
#ifdef _SCLL_
  start.clear();
  end.clear();
#endif
//  concurrent_vector <p_Particle>(particle_list).swap(particle_list);
//  particle_list.shrink_to_fit();
//  if (int(particle_list.capacity()) != 0){
//    cout<<"particle list size is wrong!!!\n";
//  }
}
#ifndef _SCLL_
//-------------------------------------------------------
// add particle to cell_list
//-------------------------------------------------------
void Cell_list::Add_particle (Particle *current_particle)
{
  particle_list.push_back(current_particle);
}
#else
//-------------------------------------------------------
// add particle and index of subcell to cell_list
//-------------------------------------------------------
void Cell_list::Add_particle_and_key (std::pair <int,p_Particle> current_particle)
{
  particle_list.push_back(current_particle);
}
//-------------------------------------------------------
// find sub cell start and end
//-------------------------------------------------------
void Cell_list::Find_cell_start_and_end (int num_cell)
{
  for(int i = 0; i < num_cell; i++){
    start.push_back(-1);
    end.push_back(-1);
  }
  int num_particle = particle_list.size();

  for(int i = 0; i < num_particle; i++){
    int hash = particle_list[i].first;
    if(i == 0 || hash != particle_list[i-1].first)
    {
      start[hash] = i;
      if ( i > 0)
        end[ particle_list[i-1].first] = i;
    }
    if (i == num_particle - 1)
    {
      end[hash] = i+1;
    }
  }
}
//-------------------------------------------------------
// build sub_cell list
//-------------------------------------------------------
void Cell_list::Build_subcell_list (int total_num_subcell)
{

  if (particle_list.size() > 0){

    sort_by_key(particle_list, total_num_subcell);

    Find_cell_start_and_end (total_num_subcell);
  }
}
#endif
