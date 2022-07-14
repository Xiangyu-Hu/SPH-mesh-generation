#include "level_infor.h"
#ifdef _INCPRS_
#include "sph_incprs.h"
#include "particle_incprs.h"
#endif
#ifdef _CPRS_
#include "sph_cprs.h"
#include "particle_cprs.h"
#endif
#ifdef _GSPH_
#include "sph_gsph.h"
#include "particle_gsph.h"
#endif
#ifdef _ALE_
#include "particle_ale.h"
#include "sph_ale.h"
#endif
#ifdef _MESH_GENERATION_
#include "particle_mesh_generation.h"
#include "sph_mesh_generation.h"
#endif
#include "cell_list.h"
#include "boost/multi_array.hpp"
typedef boost::multi_array_types::extent_range mrange;

/***************************************************/
/*                                                 */
/*     Functions defined in class "Level_info"     */
/*                                                 */
/***************************************************/

//-------------------------------------------------------
// initialze all the necessary parameters for level_info
//-------------------------------------------------------
void Level_info::Initialize(int i, SPH *sph, communicator &world)
{
  level              = i;
  Lmin               = sph->Lmin;
  Lmax               = sph->Lmax;
  total_num_particle = sph->glbl_total_num_particle;
  total_nun_color    = sph->glbl_total_num_color;
  p_cell_listpool    = &(sph->cell_listpool);
#ifdef _MPI_
  p_color_listpool   = &(sph->color_listpool);
#endif
#if PERI_DIM != 0 || SYM_DIM != 0
  p_particlepool     = &(sph->particlepool);
  my_set_const (num_cell_margin, 0);
#endif
  num_leaf_particle  = 0;
  leaf_particle.clear();
  Real multi = powern(Real(SCALE_RATIO),level);

  Real box_size = sph->domain.i;
#if DIM_Y == 1
  box_size = AMIN1(sph->domain.j, box_size);
#endif
#if DIM_Z == 1
  box_size = AMIN1(sph->domain.k, box_size);
#endif
  if ( sph->ini_scale/multi > box_size){
    cout<<"ERROR: Current level cannot be initialized!!\n";
    exit(0);
  }
  my_set_data  (domain, sph->domain);
  my_set_data  (box_l, sph->box_l);
  my_set_data  (box_r, sph->box_r);
  my_set_const (num_cell, 1);
  my_set_data  (num_cell, my_multiply_const(sph->ini_num_cell, multi));
  my_self_multiply (num_cell, total_num_cell);
  dcell = my_devide_data (sph->domain, num_cell);
#if DIM_X == 1
  scale = dcell.i;
#endif
#if DIM_Y == 1
  scale = AMIN1(dcell.j, scale);
#endif
#if DIM_Z == 1
  scale = AMIN1(dcell.k, scale);
#endif

  cell_start.i = cell_start.j = cell_start.k = 0;
  cell_end.i = cell_end.j = cell_end.k = 1;
  my_set_data (cell_end, num_cell);

  int dL = i - Lmin;
  int dN = powern(SCALE_RATIO,dL);

#if PERI_DIM_X == 1 || SYM_DIM_X == 1
  cell_start.i -= dN;
  cell_end.i += dN;
  num_cell_margin.i = dN;
#endif
#if PERI_DIM_Y == 1 || SYM_DIM_Y == 1
  cell_start.j -= dN;
  cell_end.j += dN;
  num_cell_margin.j = dN;
#endif
#if PERI_DIM_Z == 1 || SYM_DIM_Z == 1
  cell_start.k -= dN;
  cell_end.k += dN;
  num_cell_margin.k = dN;
#endif

  glbl_cell_start.i = cell_start.i;
  glbl_cell_start.j = cell_start.j;
  glbl_cell_start.k = cell_start.k;
  glbl_cell_end.i = cell_end.i;
  glbl_cell_end.j = cell_end.j;
  glbl_cell_end.k = cell_end.k;
  my_set_const (glbl_num_cell, 1);
  glbl_num_cell.i = glbl_cell_end.i - glbl_cell_start.i;
  glbl_num_cell.j = glbl_cell_end.j - glbl_cell_start.j;
  glbl_num_cell.k = glbl_cell_end.k - glbl_cell_start.k;

#ifdef _MPI_
  #if PERI_DIM_X == 1 && P_DIM_X == 1
  cell_start.i -= num_cell.i;
  cell_end.i += num_cell.i;
  #endif
  #if PERI_DIM_Y == 1 && P_DIM_Y == 1
  cell_start.j -= num_cell.j;
  cell_end.j += num_cell.j;
  #endif
  #if PERI_DIM_Z == 1 && P_DIM_Z == 1
  cell_start.k -= num_cell.k;
  cell_end.k += num_cell.k;
  #endif

  Get_local_cell_start_end(sph);

  Allocate_color_list(world);
#endif
  Allocate_memory();
  Initial_cell_list_and_reset_tags(); // reset tags to zero

#ifdef _SCLL_
  int factor = 3;
  dsubcell = my_multiply_const(dcell, 1./factor);
  my_set_const (nsubcell, factor);
  my_self_multiply(nsubcell, total_num_subcell);
  subcell_start = my_multiply_data (cell_start, nsubcell);
  subcell_end   = my_multiply_data (cell_end, nsubcell);
#endif

  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Level: "<<level<<" is initialized\n";
    cout<<"<<<<< Total number of cell    : "<<total_num_cell<<"\n";
#ifdef _SCLL_
    cout<<"<<<<< Total number of subcell : "<<total_num_subcell<<"\n";
#endif
    cout<<"<<<<< Scale size              : "<<scale<<"\n";
    cout<<"<<<<< cell start              : "<<glbl_cell_start.i<<" "<<glbl_cell_start.j<<" "<<glbl_cell_start.k<<"\n";
    cout<<"<<<<< cell end                : "<<glbl_cell_end.i<<" "<<glbl_cell_end.j<<" "<<glbl_cell_end.k<<"\n";
    cout<<"**********************************************************\n";
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Level "<<level<<" is initialized\n";
  out.close();
#endif
}
//--------------------------------------------------
// expand cell_list due to the newly 
// added coarse level
//--------------------------------------------------
void Level_info::Reinitialize_cell_list(communicator &world, SPH* sph)
{
  cell_start.i = cell_start.j = cell_start.k = 0;
  cell_end.i = cell_end.j = cell_end.k = 1;
  my_set_data (cell_end, num_cell);

  int dL = level - Lmin;
  int dN = powern(SCALE_RATIO,dL);

#if PERI_DIM != 0 || SYM_DIM != 0
  my_set_const (num_cell_margin, 0);
#endif

#if PERI_DIM_X == 1 || SYM_DIM_X == 1
  cell_start.i -= dN;
  cell_end.i += dN;
  num_cell_margin.i = dN;
#endif
#if PERI_DIM_Y == 1 || SYM_DIM_Y == 1
  cell_start.j -= dN;
  cell_end.j += dN;
  num_cell_margin.j = dN;
#endif
#if PERI_DIM_Z == 1 || SYM_DIM_Z == 1
  cell_start.k -= dN;
  cell_end.k += dN;
  num_cell_margin.k = dN;
#endif

  glbl_cell_start.i = cell_start.i;
  glbl_cell_start.j = cell_start.j;
  glbl_cell_start.k = cell_start.k;
  glbl_cell_end.i = cell_end.i;
  glbl_cell_end.j = cell_end.j;
  glbl_cell_end.k = cell_end.k;
  my_set_const (glbl_num_cell, 1);
  glbl_num_cell.i = glbl_cell_end.i - glbl_cell_start.i;
  glbl_num_cell.j = glbl_cell_end.j - glbl_cell_start.j;
  glbl_num_cell.k = glbl_cell_end.k - glbl_cell_start.k;

#ifdef _MPI_
  #if PERI_DIM_X == 1 && P_DIM_X == 1
  cell_start.i -= num_cell.i;
  cell_end.i += num_cell.i;
  #endif
  #if PERI_DIM_Y == 1 && P_DIM_Y == 1
  cell_start.j -= num_cell.j;
  cell_end.j += num_cell.j;
  #endif
  #if PERI_DIM_Z == 1 && P_DIM_Z == 1
  cell_start.k -= num_cell.k;
  cell_end.k += num_cell.k;
  #endif

  Get_local_cell_start_end(sph);

  Allocate_color_list(world);
#endif
  Allocate_memory();
  Initial_cell_list_and_reset_tags(); // reset tags to zero
#ifdef _SCLL_
  subcell_start = my_multiply_data (cell_start, nsubcell);
  subcell_end   = my_multiply_data (cell_end, nsubcell);
#endif

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Level "<<level<<" is reinitialized\n";
  out.close();
#endif
}
//-------------------------------------------------------
// Allocate_memory
//-------------------------------------------------------
void Level_info::Allocate_memory()
{
  exist_cell_list.resize
    (boost::extents[mrange(cell_start.i,cell_end.i)][mrange(cell_start.j,cell_end.j)][mrange(cell_start.k,cell_end.k)]);

  table_cell_list.resize
    (boost::extents[mrange(cell_start.i,cell_end.i)][mrange(cell_start.j,cell_end.j)][mrange(cell_start.k,cell_end.k)]);

  exist_leaf_particle.resize
    (boost::extents[mrange(cell_start.i,cell_end.i)][mrange(cell_start.j,cell_end.j)][mrange(cell_start.k,cell_end.k)]);
#ifdef _MPI_
  exchange_exist_cell_list.resize
    (boost::extents[mrange(cell_start.i,cell_end.i)][mrange(cell_start.j,cell_end.j)][mrange(cell_start.k,cell_end.k)]);

  exist_exchange_particle.resize
    (boost::extents[mrange(cell_start.i,cell_end.i)][mrange(cell_start.j,cell_end.j)][mrange(cell_start.k,cell_end.k)]);

#ifdef _NARROW_BAND_GRAPH_
  list_is_narrow_band.resize
    (boost::extents[mrange(cell_start.i,cell_end.i)][mrange(cell_start.j,cell_end.j)][mrange(cell_start.k,cell_end.k)]);
#endif
#endif
}
//-------------------------------------------------------
// reset tags to 0 for initialization
//-------------------------------------------------------
void Level_info::Initial_cell_list_and_reset_tags()
{
  static affinity_partitioner ap;

  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
                [&](const blocked_range3d<int>& r){
    for(int i=r.pages().begin(); i!=r.pages().end(); ++i){
      for(int j=r.rows().begin(); j!=r.rows().end(); ++j){
        for(int k=r.cols().begin(); k!=r.cols().end(); ++k){
          exist_cell_list[i][j][k]          = 0;
          exist_leaf_particle[i][j][k]      = 0;
#ifdef _MPI_
          exist_exchange_particle[i][j][k]  = 0;
          exchange_exist_cell_list[i][j][k] = 0;
  #ifdef _NARROW_BAND_GRAPH_
          list_is_narrow_band[i][j][k] = 0;
  #endif
#endif
          if (NULL == table_cell_list[i][j][k]){
            table_cell_list[i][j][k]   = p_cell_listpool->malloc();
          }
          table_cell_list[i][j][k]->Initialize(this);
        }
      }
    }
  }, ap);
}
//-------------------------------------------------------
// reset cell list info
//-------------------------------------------------------
void Level_info::Reset_cell_list_info(communicator &world)
{
  my_int start;
  my_int end;
  start.i = cell_start.i;
  start.j = cell_start.j;
  start.k = cell_start.k;
  end.i   = cell_end.i;
  end.j   = cell_end.j;
  end.k   = cell_end.k;

  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(start.i, end.i, start.j, end.j, start.k, end.k),
           [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++)
          table_cell_list[i][j][k]->Reset_tags();
  }, ap);
}
//-------------------------------------------------------
// Update the cell_list in current level
//-------------------------------------------------------
void Level_info::Update_cell_list(Particle *current_particle, communicator &world)
{
  my_real coord_shift = my_minus_data (current_particle->coord, box_l);
  my_int  pos = get_cell_id (coord_shift, dcell, cell_start, cell_end);
#ifdef _SCLL_
  int key = 0;
  my_real coord_shift_sub_cell = my_minus_data (coord_shift, my_multiply_data(dcell,pos));

  my_int  pos_sub = get_subcell_id (coord_shift_sub_cell, dsubcell, nsubcell);
  key = (pos_sub.k*DIM_Z*nsubcell.j + pos_sub.j*DIM_Y)*nsubcell.i + pos_sub.i*DIM_X;

  std::pair <int,p_Particle> current_particle_pair;
  current_particle_pair.first = key;
  current_particle_pair.second = current_particle;
  table_cell_list[pos.i][pos.j][pos.k]->Add_particle_and_key(current_particle_pair);
#else
  table_cell_list[pos.i][pos.j][pos.k]->Add_particle(current_particle);
#endif
  current_particle->level = level;
}
//-------------------------------------------------------
// reset tags to 0
//-------------------------------------------------------
void Level_info::Reset_tags(communicator &world)
{
  my_int start;
  my_int end;
  start.i = cell_start.i;
  start.j = cell_start.j;
  start.k = cell_start.k;
  end.i   = cell_end.i;
  end.j   = cell_end.j;
  end.k   = cell_end.k;

  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(start.i, end.i, start.j, end.j, start.k, end.k),
           [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
          exist_cell_list[i][j][k] = 0;
          exist_leaf_particle[i][j][k] = 0;
#ifdef _MPI_
          exist_exchange_particle[i][j][k] = 0;
#endif
        }
  }, ap);
}
//-------------------------------------------------------
// update the exist tag system for current level
//-------------------------------------------------------
void Level_info::Update_exist_status()
{
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
                [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++)
          if ((int)(table_cell_list[i][j][k]->particle_list.size()) > 0)
            exist_cell_list[i][j][k] = 1;
  }, ap);
}
//-------------------------------------------------------
// update the exist tag system for current level
//-------------------------------------------------------
void Level_info::Update_leaf_particle_status()
{
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
                [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++)
          if ((int)(table_cell_list[i][j][k]->particle_list.size()) > 0)
            exist_leaf_particle[i][j][k] = 1;
  }, ap);
}
//-------------------------------------------------------
// update every level info for cell list date structure 
// according to the chirldren infomation
//-------------------------------------------------------
void Level_info::Update_every_level_info(int flag, p_Level_info *level_infos, int child_level, communicator &world)
{
  array_cell_list chirldren_cell_list       = level_infos[child_level]->table_cell_list;
  array_tag       chirldren_exist_cell_list = level_infos[child_level]->exist_cell_list;

  int i_c = DIM_X==1 ? SCALE_RATIO : 1;
  int j_c = DIM_Y==1 ? SCALE_RATIO : 1;
  int k_c = DIM_Z==1 ? SCALE_RATIO : 1;

  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
                [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++){
      for(int j = r.rows().begin(); j < r.rows().end(); j++){
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
          int ii = i*SCALE_RATIO;
          int jj = j*SCALE_RATIO;
          int kk = k*SCALE_RATIO;
          for(int t=0; t<i_c; t++){
            for(int s=0; s<j_c; s++){
              for(int m=0; m<k_c; m++){
                if (chirldren_exist_cell_list[ii+t][jj+s][kk+m] == 1){
                  int num_particle_in_chirldren = int(chirldren_cell_list[ii+t][jj+s][kk+m]->particle_list.size());
                  for(int icyc=0; icyc<num_particle_in_chirldren; icyc++){
#ifdef _SCLL_
                    Particle *current = chirldren_cell_list[ii+t][jj+s][kk+m]->particle_list[icyc].second;
                    int key = 0;
                    my_int pos;
                    pos.i = i; pos.j = j; pos.k = k;
                    my_real coord_shift = my_minus_data (current->coord, box_l);
                    my_real coord_shift_sub_cell = my_minus_data (coord_shift, my_multiply_data(dcell, pos));
                    my_int  pos_sub = get_subcell_id (coord_shift_sub_cell, dsubcell, nsubcell);
                    key = (pos_sub.k*DIM_Z*nsubcell.j + pos_sub.j*DIM_Y)*nsubcell.i + pos_sub.i*DIM_X;

                    std::pair <int,p_Particle> current_particle_pair;
                    current_particle_pair.first = key;
                    current_particle_pair.second = current;
                    table_cell_list[pos.i][pos.j][pos.k]->Add_particle_and_key(current_particle_pair);
#else
                    table_cell_list[i][j][k]->Add_particle(chirldren_cell_list[ii+t][jj+s][kk+m]->particle_list[icyc]);
#endif
                  }
                }
              }
            }
          }
        }
      }
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_level "<<level<<" information finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// refresh the neighbor information for each particle
//-------------------------------------------------------
#ifndef _SCLL_
void Level_info::Refresh_neighbor_info(SPH *sph, Particle *current, my_int cell_id, int status)
{
  int     i_c = DIM_X==1 ? 1 : 0;
  int     j_c = DIM_Y==1 ? 1 : 0;
  int     k_c = DIM_Z==1 ? 1 : 0;

  if (exist_cell_list[cell_id.i][cell_id.j][cell_id.k] == 1 && exist_leaf_particle[cell_id.i][cell_id.j][cell_id.k] == 1){
    if (int(table_cell_list[cell_id.i][cell_id.j][cell_id.k]->particle_list.size()) == 0){
      cout<<"Warning: the tag system for cell list is not consistent !!!"<<endl;
      exit(0);
    }
  }
    for(int t=AMAX1((cell_id.i-i_c),cell_start.i); t<=AMIN1((cell_id.i+i_c),(cell_end.i-1)); t++){
      for(int s=AMAX1((cell_id.j-j_c),cell_start.j); s<=AMIN1((cell_id.j+j_c),(cell_end.j-1)); s++){
        for(int m=AMAX1((cell_id.k-k_c),cell_start.k); m<=AMIN1((cell_id.k+k_c),(cell_end.k-1)); m++){
          if (exist_cell_list[t][s][m] == 1){
            int p_list_size = table_cell_list[t][s][m]->particle_list.size();
            if (p_list_size == 0){
              cout<<"Warning: the tag system for cell list is not consistent !!!"<<endl;
              exit(0);
            }
            for(int jcyc=0;jcyc<p_list_size;jcyc++){
              Particle *neighbor = table_cell_list[t][s][m]->particle_list[jcyc];
              int flag = 0;
              my_real dr = my_minus_data(current->coord, neighbor->coord);
              Real  dist = get_distance(dr);
              if( (dist <= (current->h*CUT_OFF + 1.e-10)) || (dist <= (neighbor->h*CUT_OFF + 1.e-10)) ){
                current->Add_neighbor(neighbor);
                if(neighbor->level != level){
                  neighbor->Add_neighbor(current);
                }
              }
            }
          }
        }
      }
    }
}
#else
void Level_info::Refresh_neighbor_info(SPH *sph, Particle *current, my_int scell_id, int status)
{
    Real   interaction_range = current->h*CUT_OFF;
    my_int sub_cell_range;
    sub_cell_range.i = DIM_X == 1 ? int ((interaction_range+dsubcell.i-1.e-10)/dsubcell.i) : 0;
    sub_cell_range.j = DIM_Y == 1 ? int ((interaction_range+dsubcell.j-1.e-10)/dsubcell.j) : 0;
    sub_cell_range.k = DIM_Z == 1 ? int ((interaction_range+dsubcell.k-1.e-10)/dsubcell.k) : 0;

    my_int  sj; my_set_const (sj, 0);
    my_int  s_start, s_end;
    s_start.i = AMAX1((scell_id.i-sub_cell_range.i),subcell_start.i);
    s_start.j = AMAX1((scell_id.j-sub_cell_range.j),subcell_start.j);
    s_start.k = AMAX1((scell_id.k-sub_cell_range.k),subcell_start.k);
    s_end.i   = AMIN1((scell_id.i+sub_cell_range.i+1),subcell_end.i);
    s_end.j   = AMIN1((scell_id.j+sub_cell_range.j+1),subcell_end.j);
    s_end.k   = AMIN1((scell_id.k+sub_cell_range.k+1),subcell_end.k);

    for(sj.i=s_start.i; sj.i<s_end.i; sj.i++){
      for(sj.j=s_start.j; sj.j<s_end.j; sj.j++){
        for(sj.k=s_start.k; sj.k<s_end.k; sj.k++){

          my_int  sub_j;
          my_int  targ_cell_id;

          targ_cell_id.i = DIM_X == 1 ?  floor(Real(sj.i)/nsubcell.i) : 0;
          targ_cell_id.j = DIM_Y == 1 ?  floor(Real(sj.j)/nsubcell.j) : 0;
          targ_cell_id.k = DIM_Z == 1 ?  floor(Real(sj.k)/nsubcell.k) : 0;

          sub_j.i = DIM_X == 1 ? sj.i - nsubcell.i*targ_cell_id.i : 0;
          sub_j.j = DIM_Y == 1 ? sj.j - nsubcell.j*targ_cell_id.j : 0;
          sub_j.k = DIM_Z == 1 ? sj.k - nsubcell.k*targ_cell_id.k : 0;

            if (exist_cell_list[targ_cell_id.i][targ_cell_id.j][targ_cell_id.k] == 1){
              if (int(table_cell_list[targ_cell_id.i][targ_cell_id.j][targ_cell_id.k]->particle_list.size()) == 0){
                cout<<"Warning: the tag system for cell list is not consistent !!!"<<endl;
                exit(0);
              }
              int  key_j = (sub_j.k*DIM_Z*nsubcell.j + sub_j.j*DIM_Y)*nsubcell.i + sub_j.i*DIM_X;
              int  start_index = table_cell_list[targ_cell_id.i][targ_cell_id.j][targ_cell_id.k]->start[key_j];
              if (start_index != -1){
                int end_index = table_cell_list[targ_cell_id.i][targ_cell_id.j][targ_cell_id.k]->end[key_j];
                for (int i = start_index; i < end_index; i ++){
                  Particle *neighbor = table_cell_list[targ_cell_id.i][targ_cell_id.j][targ_cell_id.k]->particle_list[i].second;
                  my_real dr = my_minus_data(current->coord, neighbor->coord);
                  Real  dist = get_distance(dr);
                  if((dist <= (current->h*CUT_OFF + 1.e-10)) || (dist <= (neighbor->h*CUT_OFF + 1.e-10))){
                    if ( current->h > neighbor->h ){
                      current->Add_neighbor(neighbor);
                      neighbor->Add_neighbor(current);
                    }else if ( current->h == neighbor->h ){
                      current->Add_neighbor(neighbor);
                    }
                  }
                }
              }
          }
        }
      }
    }
}
#endif
//-------------------------------------------------------
// interaction between two neighbor particle 
// for the neighbor information
//-------------------------------------------------------
void Level_info::Interaction
(Particle *current, Particle *neighbor, my_real dr, Real dist, int flag, int status)
{
  if( (dist <= (current->h*CUT_OFF + 1.e-10)) || (dist <= (neighbor->h*CUT_OFF + 1.e-10)) )
  {
    if(flag == 1)
    {
      current->Add_neighbor(neighbor);
      neighbor->Add_neighbor(current);

    }else
    {
      current->Add_neighbor(neighbor);
    }
  }
}
//--------------------------------------------------
// output current level 
//--------------------------------------------------
void Level_info::Output_level(communicator &world, char *filename, int n)
{
  ofstream out(filename, ios::app);

  //defining header for tecplot(plot software)
  out<<"\n";
  out<<"title='View'"<<"\n";
  out<<"variables=";
#if DIM_X
  out<<"x, ";
#endif
#if DIM_Y
  out<<"y, ";
#endif
#if DIM_Z
  out<<"z, ";
#endif

  int LEN;
#ifndef _NARROW_BAND_GRAPH_
  LEN = 6;
  out<<"level, color, leaf, exist, exist_exc, np"<<"\n";
#else
  LEN = 7;
  out<<"level, color, leaf, exist, exist_exc, np, narrow_band"<<"\n";
#endif

  my_int start;
  my_int end;
 
  start.i = cell_start.i;
  start.j = cell_start.j;
  start.k = cell_start.k;
  end.i   = cell_end.i;
  end.j   = cell_end.j;
  end.k   = cell_end.k;

  out<<"zone t='sub_filed_"<<n<<"_"<<world.rank()<<"_"<<level-Lmin<<"'"
     <<"  i= "<<end.i-start.i+1
     <<"  j= "<<end.j-start.j+1
     <<"  k= "<<end.k-start.k+1
     <<"  DATAPACKING=BLOCK, VARLOCATION=([";
  int pos_s = DIM_X+DIM_Y+DIM_Z+1;
  out<<pos_s<<"-";
  out<<2*pos_s -1 + LEN-1<<"]=CELLCENTERED) SOLUTIONTIME="<<n<<"\n";

#if DIM_X
  for(int k=start.k; k<=end.k; k++){
    for(int j=start.j; j<=end.j; j++){
      for(int i=start.i; i<=end.i; i++){
        out<<dcell.i*i+box_l.i<<" ";
      }
      out<<"\n";
    }
  }
#endif
#if DIM_Y
  for(int k=start.k; k<=end.k; k++){
    for(int j=start.j; j<=end.j; j++){
      for(int i=start.i; i<=end.i; i++){
        out<<dcell.j*j+box_l.j<<" ";
      }
      out<<"\n";
    }
  }
#endif
#if DIM_Z
  for(int k=start.k; k<=end.k; k++){
    for(int j=start.j; j<=end.j; j++){
      for(int i=start.i; i<=end.i; i++){
        out<<dcell.k*k+box_l.k<<" ";
      }
      out<<"\n";
    }
  }
#endif

  // level
  for(int k=start.k; k<end.k; k++){
    for(int j=start.j; j<end.j; j++){
      for(int i=start.i; i<end.i; i++){
        out<<level<<" ";
      }
      out<<"\n";
    }
  }

  // color
  for(int k=start.k; k<end.k; k++){
    for(int j=start.j; j<end.j; j++){
      for(int i=start.i; i<end.i; i++){
        out<<world.rank()<<" ";
      }
      out<<"\n";
    }
  }

  // leaf
  for(int k=start.k; k<end.k; k++){
    for(int j=start.j; j<end.j; j++){
      for(int i=start.i; i<end.i; i++){
        if (exist_leaf_particle[i][j][k] == 1)
          out<<"1 ";
        else if (exist_leaf_particle[i][j][k] == 0)
          out<<"0 ";
      }
      out<<"\n";
    }
  }

  // exist
  for(int k=start.k; k<end.k; k++){
    for(int j=start.j; j<end.j; j++){
      for(int i=start.i; i<end.i; i++){
        if (exist_cell_list[i][j][k] == 1)
          out<<"1 ";
        else if (exist_cell_list[i][j][k] == 0)
          out<<"0 ";
      }
      out<<"\n";
    }
  }

  // exist_exc
  for(int k=start.k; k<end.k; k++){
    for(int j=start.j; j<end.j; j++){
      for(int i=start.i; i<end.i; i++){
#ifdef _MPI_
        if (exist_exchange_particle[i][j][k] == 1)
          out<<"1 ";
        else if (exist_exchange_particle[i][j][k] == 0)
          out<<"0 ";
#endif
#ifndef _MPI_
        out<<"0 ";
#endif
      }
      out<<"\n";
    }
  }

  // np
  for(int k=start.k; k<end.k; k++){
    for(int j=start.j; j<end.j; j++){
      for(int i=start.i; i<end.i; i++){
        out<<table_cell_list[i][j][k]->particle_list.size()<<" ";
      }
      out<<"\n";
    }
  }
#ifdef _NARROW_BAND_GRAPH_
  // np
  for(int k=start.k; k<end.k; k++){
    for(int j=start.j; j<end.j; j++){
      for(int i=start.i; i<end.i; i++){
        if (list_is_narrow_band[i][j][k] == 1)
          out<<"1 ";
        else if (list_is_narrow_band[i][j][k] == 0)
          out<<"0 ";
      }
      out<<"\n";
    }
  }
#endif
  out.close();

}
#ifdef _SCLL_
//--------------------------------------------------
// build and sort subcell list
//--------------------------------------------------
void Level_info::Build_subcell_list (communicator &world){
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
                [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
          table_cell_list[i][j][k]->Build_subcell_list(total_num_subcell);
        }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Build_subcell_list in "<<level<<" finished\n";
  out.close();
#endif
}
#endif
#if SYM_DIM != 0
//--------------------------------------------------
// construct symmetric particles in buffer area
//--------------------------------------------------
void Level_info::Construct_symmetric_BC_particles
(concurrent_vector  <p_Particle> &particle_sym)
{
  // find cells that are outside of physical domain
  // and copy all the particles belonging to corresponding
  // cell to current sym_boundary cells
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
                [&](const blocked_range3d<int>& r){
    int t, s, m; t = s = m = 0; int tag = 0;
    my_real pos; pos.i = pos.j = pos.k = 0.;
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
          if (i < -0.5 || i > num_cell.i-0.5)
            {t = num_cell.i + i/abs(i)*num_cell.i - i-1; tag = 1; pos.i = i/abs(i);}
          else { t = i; pos.i = 0.;}
          if (j < -0.5 || j > num_cell.j-0.5) 
            {s = num_cell.j + j/abs(j)*num_cell.j - j-1; tag = 1; pos.j = j/abs(j);}
          else {s = j; pos.j = 0.;}
          if (k < -0.5 || k > num_cell.k-0.5) 
            {m = num_cell.k + k/abs(k)*num_cell.k - k-1; tag = 1; pos.k = k/abs(k);}
          else {m = k; pos.k = 0.;}
          int p_list_size = table_cell_list[t][s][m]->particle_list.size();
          if (tag == 1 && p_list_size != 0){
            for(int icyc=0; icyc<p_list_size; icyc++){
#ifndef _SCLL_
              p_Particle cp = table_cell_list[t][s][m]->particle_list[icyc];
#else
              p_Particle cp = table_cell_list[t][s][m]->particle_list[icyc].second;
#endif
              p_Particle particle_new = new Particle;

	      my_real shift;
	      
	      if (pos.i < -0.5) shift.i = 2.*box_l.i - cp->coord.i;
	      else if (pos.i > 0.5) shift.i = 2.*box_r.i - cp->coord.i;
	      else shift.i = cp->coord.i;
	      if (pos.j < -0.5) shift.j = 2.*box_l.j - cp->coord.j;
	      else if (pos.j > 0.5) shift.j = 2.*box_r.j - cp->coord.j;
	      else shift.j = cp->coord.j;
	      if (pos.k < -0.5) shift.k = 2.*box_l.k - cp->coord.k;
	      else if (pos.k > 0.5) shift.k = 2.*box_r.k - cp->coord.k;
	      else shift.k = cp->coord.k;

              particle_new->Set_bc_particle_info(cp, shift, pos);
              particle_sym.push_back(particle_new);
#ifndef _SCLL_
              table_cell_list[i][j][k]->Add_particle(particle_new);
#else
              int key = 0;
              my_int  cell_id; cell_id.i = i; cell_id.j = j; cell_id.k = k;
              my_real coord_shift = my_minus_data (particle_new->coord, box_l);
              my_real coord_shift_sub_cell = my_minus_data (coord_shift, my_multiply_data(dcell,cell_id));

              my_int  pos_sub = get_subcell_id (coord_shift_sub_cell, dsubcell, nsubcell);
              key = (pos_sub.k*DIM_Z*nsubcell.j + pos_sub.j*DIM_Y)*nsubcell.i + pos_sub.i*DIM_X;

              std::pair <int,p_Particle> current_particle_pair;
              current_particle_pair.first = key;
              current_particle_pair.second = particle_new;
              table_cell_list[i][j][k]->Add_particle_and_key(current_particle_pair);
#endif
            }
          }
          tag = 0;
        }
  }, ap);
  
}
#endif
#if PERI_DIM != 0
//--------------------------------------------------
// construct periodical particles in buffer area
//--------------------------------------------------
void Level_info::Construct_periodical_BC_particles_1
(concurrent_vector  <p_Particle> &particle_peri)
{
  // find cells that are outside of physical domain
  // and copy all the particles belonging to corresponding
  // cell to current peri_boundary cells
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
                [&](const blocked_range3d<int>& r){
    int t, s, m; t = s = m = 0; int tag = 0;
    my_real pos; pos.i = pos.j = pos.k = 0.;
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
#ifndef _MPI_
          if (i < -0.5 || i > num_cell.i-0.5)
            {t = i - i/abs(i)*num_cell.i; tag = 1; pos.i = i/abs(i)*domain.i;}
          else { t = i; pos.i = 0.;}
          if (j < -0.5 || j > num_cell.j-0.5) 
            {s = j - j/abs(j)*num_cell.j; tag = 1; pos.j = j/abs(j)*domain.j;}
          else {s = j; pos.j = 0.;}
          if (k < -0.5 || k > num_cell.k-0.5) 
            {m = k - k/abs(k)*num_cell.k; tag = 1; pos.k = k/abs(k)*domain.k;}
          else {m = k; pos.k = 0.;}
#endif
#ifdef _MPI_
  #if PERI_DIM_X == 1
    #if P_DIM_X == 0
          if (i < -0.5 || i > num_cell.i-0.5)
            {t = i - i/abs(i)*num_cell.i; tag = 1; pos.i = i/abs(i)*domain.i;}
          else { t = i; pos.i = 0.;}
    #elif P_DIM_X == 1
          t = i; pos.i = 0.;
    #endif
  #endif
  #if PERI_DIM_Y == 1
    #if P_DIM_Y == 0
          if (j < -0.5 || j > num_cell.j-0.5)
            {s = j - j/abs(j)*num_cell.j; tag = 1; pos.j = j/abs(j)*domain.j;}
          else {s = j; pos.j = 0.;}
    #elif P_DIM_Y == 1
          s = j; pos.j = 0.;
    #endif
  #endif
  #if PERI_DIM_Z == 1
    #if P_DIM_Z == 0
          if (k < -0.5 || k > num_cell.k-0.5)
            {m = k - k/abs(k)*num_cell.k; tag = 1; pos.k = k/abs(k)*domain.k;}
          else {m = k; pos.k = 0.;}
    #elif P_DIM_Z == 1
          m = k; pos.k = 0.;
    #endif
  #endif
#endif
          int p_list_size = table_cell_list[t][s][m]->particle_list.size();
          if (tag == 1 && p_list_size != 0){
            for(int icyc=0; icyc<p_list_size; icyc++){
#ifndef _SCLL_
              p_Particle cp = table_cell_list[t][s][m]->particle_list[icyc];
#else
              p_Particle cp = table_cell_list[t][s][m]->particle_list[icyc].second;
#endif
              p_Particle particle_new = new Particle;
              particle_new->Set_bc_particle_info(cp, pos);
              particle_peri.push_back(particle_new);

#ifndef _SCLL_
              table_cell_list[i][j][k]->Add_particle(particle_new);
#else
              int key = 0;
              my_int  cell_id; cell_id.i = i; cell_id.j = j; cell_id.k = k;
              my_real coord_shift = my_minus_data (particle_new->coord, box_l);
              my_real coord_shift_sub_cell = my_minus_data (coord_shift, my_multiply_data(dcell,cell_id));

              my_int  pos_sub = get_subcell_id (coord_shift_sub_cell, dsubcell, nsubcell);
              key = (pos_sub.k*DIM_Z*nsubcell.j + pos_sub.j*DIM_Y)*nsubcell.i + pos_sub.i*DIM_X;

              std::pair <int,p_Particle> current_particle_pair;
              current_particle_pair.first = key;
              current_particle_pair.second = particle_new;
              table_cell_list[i][j][k]->Add_particle_and_key(current_particle_pair);
#endif
            }
          }
          tag = 0;
        }
  }, ap);
  
}
#ifdef _MPI_
//--------------------------------------------------
// construct periodical particles in buffer area
//--------------------------------------------------
void Level_info::Construct_periodical_BC_particles_2
(concurrent_vector  <p_Particle> &particle_peri, communicator &world)
{
  // find particles missing when partitioning results
  // features self interacting
  int i_c = (DIM_X==1 && PERI_DIM_X == 1 && P_DIM_X == 1) ? 1 : 0;
  int j_c = (DIM_Y==1 && PERI_DIM_Y == 1 && P_DIM_Y == 1) ? 1 : 0;
  int k_c = (DIM_Z==1 && PERI_DIM_Z == 1 && P_DIM_Z == 1) ? 1 : 0;

  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
                [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
          int tag = 0;
          my_real pos; pos.i = pos.j = pos.k = 0.;
          for( int r= AMAX1((i-i_c),cell_start.i); r<= AMIN1((i+i_c),(cell_end.i-1)); r++)
            for( int s= AMAX1((j-j_c),cell_start.j); s<= AMIN1((j+j_c),(cell_end.j-1)); s++)
              for( int t= AMAX1((k-k_c),cell_start.k); t<= AMIN1((k+k_c),(cell_end.k-1)); t++){
                int np_i = table_cell_list[r][s][t]->particle_list.size();
                if ( np_i > 0) tag = 1;
              }
              if (tag == 1){
                for (int rr = -1*i_c*num_cell.i; rr <= i_c*num_cell.i; rr+=num_cell.i)
                  for (int ss = -1*j_c*num_cell.j; ss <= j_c*num_cell.j; ss+=num_cell.j)
                    for (int tt = -1*k_c*num_cell.k; tt <= k_c*num_cell.k; tt+=num_cell.k){
                      if (rr+ss+tt != 0){
                        int id_i_t = i + rr;
                        int id_j_t = j + ss;
                        int id_k_t = k + tt;
                        if (id_i_t >= cell_start.i && id_i_t < cell_end.i &&
                            id_j_t >= cell_start.j && id_j_t < cell_end.j &&
                            id_k_t >= cell_start.k && id_k_t < cell_end.k){
                          pos.i = -1 * Real(rr/num_cell.i) * domain.i;
                          pos.j = -1 * Real(ss/num_cell.j) * domain.j;
                          pos.k = -1 * Real(tt/num_cell.k) * domain.k;

                          int np_j = table_cell_list[id_i_t][id_j_t][id_k_t]->particle_list.size();
                          if ( np_j > 0){
                            for(int icyc=0; icyc<np_j; icyc++){
#ifndef _SCLL_
                              p_Particle cp = table_cell_list[id_i_t][id_j_t][id_k_t]->particle_list[icyc];
#else
                              p_Particle cp = table_cell_list[id_i_t][id_j_t][id_k_t]->particle_list[icyc].second;
#endif
                              p_Particle particle_new = new Particle;
                              particle_new->Set_bc_particle_info(cp, pos);
                              particle_peri.push_back(particle_new);
                            }
                          }
                        }
                      }
                    }
              }
        }
  }, ap);
}
#endif
#endif
#ifdef _MPI_
#if defined(_DIST_NUMBER_) || defined(_WEIGHTED_PARTITION_)
//-------------------------------------------------------
// cout the number of dist calculation
// for current particle
//-------------------------------------------------------
#ifndef _SCLL_
void Level_info::Get_patitioning_mass(SPH *sph, Particle_base *current, my_int cell_id, int status)
{

  int     i_c = DIM_X==1 ? 1 : 0;
  int     j_c = DIM_Y==1 ? 1 : 0;
  int     k_c = DIM_Z==1 ? 1 : 0;

  if(current->level == level){
    for(int t=AMAX1((cell_id.i-i_c),cell_start.i); t<=AMIN1((cell_id.i+i_c),(cell_end.i-1)); t++){
      for(int s=AMAX1((cell_id.j-j_c),cell_start.j); s<=AMIN1((cell_id.j+j_c),(cell_end.j-1)); s++){
        for(int m=AMAX1((cell_id.k-k_c),cell_start.k); m<=AMIN1((cell_id.k+k_c),(cell_end.k-1)); m++){
          if (exist_cell_list[t][s][m] == 1){
            current->p_mass += table_cell_list[t][s][m]->particle_list.size();
          }
        }
      }
    }
  }
}
#else
void Level_info::Get_patitioning_mass(SPH *sph, Particle_base *current, my_int scell_id, int status)
{
  Real   interaction_range = current->h*CUT_OFF;
  my_int sub_cell_range;
  sub_cell_range.i = DIM_X == 1 ? int ((interaction_range+dsubcell.i-1.e-10)/dsubcell.i) : 0;
  sub_cell_range.j = DIM_Y == 1 ? int ((interaction_range+dsubcell.j-1.e-10)/dsubcell.j) : 0;
  sub_cell_range.k = DIM_Z == 1 ? int ((interaction_range+dsubcell.k-1.e-10)/dsubcell.k) : 0;

  my_int  sj; my_set_const (sj, 0);
  my_int  s_start, s_end;
  s_start.i = AMAX1((scell_id.i-sub_cell_range.i),subcell_start.i);
  s_start.j = AMAX1((scell_id.j-sub_cell_range.j),subcell_start.j);
  s_start.k = AMAX1((scell_id.k-sub_cell_range.k),subcell_start.k);
  s_end.i   = AMIN1((scell_id.i+sub_cell_range.i+1),subcell_end.i);
  s_end.j   = AMIN1((scell_id.j+sub_cell_range.j+1),subcell_end.j);
  s_end.k   = AMIN1((scell_id.k+sub_cell_range.k+1),subcell_end.k);

  for(sj.i=s_start.i; sj.i<s_end.i; sj.i++){
    for(sj.j=s_start.j; sj.j<s_end.j; sj.j++){
      for(sj.k=s_start.k; sj.k<s_end.k; sj.k++){

        my_int  sub_j;
        my_int  targ_cell_id;
        targ_cell_id.i = DIM_X == 1 ?  floor(Real(sj.i)/nsubcell.i) : 0;
        targ_cell_id.j = DIM_Y == 1 ?  floor(Real(sj.j)/nsubcell.j) : 0;
        targ_cell_id.k = DIM_Z == 1 ?  floor(Real(sj.k)/nsubcell.k) : 0;

        sub_j.i = DIM_X == 1 ? sj.i - nsubcell.i*targ_cell_id.i : 0;
        sub_j.j = DIM_Y == 1 ? sj.j - nsubcell.j*targ_cell_id.j : 0;
        sub_j.k = DIM_Z == 1 ? sj.k - nsubcell.k*targ_cell_id.k : 0;

        if (exist_cell_list[targ_cell_id.i][targ_cell_id.j][targ_cell_id.k] == 1){
          if (int(table_cell_list[targ_cell_id.i][targ_cell_id.j][targ_cell_id.k]->particle_list.size()) == 0){
            cout<<"Warning: the tag system for cell list is not consistent !!!"<<endl;
            exit(0);
          }
          int  key_j = (sub_j.k*DIM_Z*nsubcell.j + sub_j.j*DIM_Y)*nsubcell.i + sub_j.i*DIM_X;
          int  start_index = table_cell_list[targ_cell_id.i][targ_cell_id.j][targ_cell_id.k]->start[key_j];
          if (start_index != -1){
            int end_index = table_cell_list[targ_cell_id.i][targ_cell_id.j][targ_cell_id.k]->end[key_j];
            current->p_mass += (end_index - start_index);
          }
        }
      }
    }
  }
}
#endif
#endif
//-------------------------------------------------------
// allocate memory for color list
//-------------------------------------------------------
void Level_info::Allocate_color_list(communicator &world)
{
  my_int start;
  my_int end;
  if (world.rank() == 0){
    start.i = glbl_cell_start.i;
    start.j = glbl_cell_start.j;
    start.k = glbl_cell_start.k;
    end.i   = glbl_cell_end.i;
    end.j   = glbl_cell_end.j;
    end.k   = glbl_cell_end.k;
  }else{
    start.i = 0;
    start.j = 0;
    start.k = 0;
    end.i   = 1;
    end.j   = 1;
    end.k   = 1;
  }
  color_list.resize
  (boost::extents[mrange(start.i,end.i)][mrange(start.j,end.j)][mrange(start.k,end.k)]);

  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(start.i, end.i, start.j, end.j, start.k, end.k),
           [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
          if (NULL == color_list[i][j][k]){
            color_list[i][j][k] = p_color_listpool->malloc();
          }
          color_list[i][j][k]->Clear_data();
        }
  }, ap);
}
#ifdef _NARROW_BAND_GRAPH_
//-------------------------------------------------------
// Reset narrow band list
//-------------------------------------------------------
void Level_info::Reset_narrow_band_list(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
           [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
            list_is_narrow_band[i][j][k] = 0;
        }
  }, ap);
}

//-------------------------------------------------------
// find cells that belongs to the narrow band
//-------------------------------------------------------
void Level_info::Find_narrow_band_cells()
{
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
           [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
          if (exist_leaf_particle[i][j][k] == 1){
            int is_narrow_band_cell = 0;
            if ((DIM_X == 1 && (i == cell_start.i || i == cell_end.i-1)) ||
                (DIM_Y == 1 && (j == cell_start.j || j == cell_end.j-1)) ||
                (DIM_Z == 1 && (k == cell_start.k || k == cell_end.k-1)) ) is_narrow_band_cell = 1;
            for( int r=-DIM_X; r<=DIM_X; r++)
              for( int s=-DIM_Y; s<=DIM_Y; s++)
                for( int t=-DIM_Z; t<=DIM_Z; t++){
                  int index_i= AMAX1(cell_start.i,AMIN1(i+r,cell_end.i-1));
                  int index_j= AMAX1(cell_start.j,AMIN1(j+s,cell_end.j-1));
                  int index_k= AMAX1(cell_start.k,AMIN1(k+t,cell_end.k-1));
                  if (exist_leaf_particle[index_i][index_j][index_k] == 0){
                    is_narrow_band_cell = 1;
                  }
                }
            if (is_narrow_band_cell == 1){
              list_is_narrow_band[i][j][k] = 1;
            }
          }
        }
  }, ap);
}
#endif
//-------------------------------------------------------
// reset color list
//-------------------------------------------------------
void Level_info::Reset_color_list()
{
  my_int start;
  my_int end;
  start.i = glbl_cell_start.i;
  start.j = glbl_cell_start.j;
  start.k = glbl_cell_start.k;
  end.i   = glbl_cell_end.i;
  end.j   = glbl_cell_end.j;
  end.k   = glbl_cell_end.k;

  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(start.i, end.i, start.j, end.j, start.k, end.k),
           [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
          color_list[i][j][k]->Clear_data();
        }
  }, ap);
}
//-------------------------------------------------------
// Get_local_cell_start_end
//-------------------------------------------------------
void Level_info::Get_local_cell_start_end(SPH *sph)
{
  int dL = level - Lmin;
  int dN = powern(SCALE_RATIO,dL);

  my_set_data  (local_box_l, sph->local_box_l);
  my_set_data  (local_box_r, sph->local_box_r);

  my_real coord_shift_l = my_minus_data (local_box_l, box_l);
  my_real coord_shift_r = my_minus_data (local_box_r, box_l);
  my_real cell_size     = my_multiply_const (dcell, dN);
  my_int  pos_l = get_cell_id (coord_shift_l, cell_size, cell_start, cell_end);
  my_int  pos_r = get_cell_id (coord_shift_r, cell_size, cell_start, cell_end);

  int ii = 0;
  int jj = 0;
  int kk = 0;

#if PERI_DIM != 0
  ii = PERI_DIM_X == 1? 1:0; 
  jj = PERI_DIM_Y == 1? 1:0;
  kk = PERI_DIM_Z == 1? 1:0;
#endif
#if SYM_DIM != 0
  ii = SYM_DIM_X == 1? 1:0;
  jj = SYM_DIM_Y == 1? 1:0;
  kk = SYM_DIM_Z == 1? 1:0;
#endif

#if SYM_DIM != 0
  pos_l.i = DIM_X == 1? AMAX1(pos_l.i-2, 0-ii):0;
  pos_l.j = DIM_Y == 1? AMAX1(pos_l.j-2, 0-jj):0;
  pos_l.k = DIM_Z == 1? AMAX1(pos_l.k-2, 0-kk):0;
  pos_r.i = DIM_X == 1? AMIN1(pos_r.i+3, num_cell.i/dN+ii):1;
  pos_r.j = DIM_Y == 1? AMIN1(pos_r.j+3, num_cell.j/dN+jj):1;
  pos_r.k = DIM_Z == 1? AMIN1(pos_r.k+3, num_cell.k/dN+kk):1;
#endif

#if PERI_DIM != 0
  #if P_DIM_X == 0 && PERI_DIM_X == 1
    pos_l.i = DIM_X == 1? AMAX1(pos_l.i-2, 0-ii):0;
    pos_r.i = DIM_X == 1? AMIN1(pos_r.i+3, num_cell.i/dN+ii):1;
  #elif P_DIM_X == 1 && PERI_DIM_X == 1
    pos_l.i = DIM_X == 1? pos_l.i-2:0;
    pos_r.i = DIM_X == 1? pos_r.i+3:1;
  #endif

  #if P_DIM_Y == 0 && PERI_DIM_Y == 1
    pos_l.j = DIM_Y == 1? AMAX1(pos_l.j-2, 0-jj):0;
    pos_r.j = DIM_Y == 1? AMIN1(pos_r.j+3, num_cell.j/dN+jj):1;
  #elif P_DIM_Y == 1 && PERI_DIM_Y == 1
    pos_l.j = DIM_Y == 1? pos_l.j-2:0;
    pos_r.j = DIM_Y == 1? pos_r.j+3:1;
  #endif

  #if P_DIM_Z == 0 && PERI_DIM_Z == 1
    pos_l.k = DIM_Z == 1? AMAX1(pos_l.k-2, 0-kk):0;
    pos_r.k = DIM_Z == 1? AMIN1(pos_r.k+3, num_cell.k/dN+kk):1;
  #elif P_DIM_Z == 1 && PERI_DIM_Z == 1
    pos_l.k = DIM_Z == 1? pos_l.k-2:0;
    pos_r.k = DIM_Z == 1? pos_r.k+3:1;
  #endif
#endif

#if PERI_DIM == 0 && SYM_DIM == 0
  pos_l.i = DIM_X == 1? AMAX1(pos_l.i-2, 0):0;
  pos_l.j = DIM_Y == 1? AMAX1(pos_l.j-2, 0):0;
  pos_l.k = DIM_Z == 1? AMAX1(pos_l.k-2, 0):0;
  pos_r.i = DIM_X == 1? AMIN1(pos_r.i+3, num_cell.i/dN):1;
  pos_r.j = DIM_Y == 1? AMIN1(pos_r.j+3, num_cell.j/dN):1;
  pos_r.k = DIM_Z == 1? AMIN1(pos_r.k+3, num_cell.k/dN):1;
#endif

  pos_l = my_multiply_const (pos_l, dN);
  pos_r = my_multiply_const (pos_r, dN);

  cell_start.i = DIM_X == 1? AMAX1(pos_l.i,cell_start.i):cell_start.i;
  cell_start.j = DIM_Y == 1? AMAX1(pos_l.j,cell_start.j):cell_start.j;
  cell_start.k = DIM_Z == 1? AMAX1(pos_l.k,cell_start.k):cell_start.k;
  cell_end.i   = DIM_X == 1? AMIN1(pos_r.i,cell_end.i):cell_end.i;
  cell_end.j   = DIM_Y == 1? AMIN1(pos_r.j,cell_end.j):cell_end.j;
  cell_end.k   = DIM_Z == 1? AMIN1(pos_r.k,cell_end.k):cell_end.k;
}
//-------------------------------------------------------
// reset cell list info
//-------------------------------------------------------
void Level_info::Initialize_color_infor()
{
  // set color list in master node
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
                [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
#ifndef _NARROW_BAND_GRAPH_
          if(exist_leaf_particle[i][j][k] == 1){
#else
          if(exist_leaf_particle[i][j][k] == 1 && list_is_narrow_band[i][j][k] == 1){
#endif            
            std:: pair <int,int> pair_t;
            pair_t.first = 0;
            pair_t.second = level - Lmin;

            int t, s, m; t = s = m = 0;
            Shift_cell_index(t, s, m, i, j, k);
            color_list[t][s][m]->Add_color(pair_t);
          }
        }
  }, ap);
}
//-------------------------------------------------------
// Shift cell index in periodical tree
//-------------------------------------------------------
void Level_info::Shift_cell_index(int &t, int &s, int &m, int i, int j, int k)
{
#if PERI_DIM != 0
  t = i; s = j; m = k;

  #if P_DIM_X == 1 && PERI_DIM_X == 1
    if (i >= num_cell.i) t = i - num_cell.i;
    else if (i < 0) t = i + num_cell.i;
    else t = i;
  #endif

  #if P_DIM_Y == 1 && PERI_DIM_Y == 1
    if (j >= num_cell.j) s = j - num_cell.j;
    else if (j < 0) s = j + num_cell.j;
    else s = j;
  #endif

  #if P_DIM_Z == 1 && PERI_DIM_Z == 1
    if (k >= num_cell.k) m = k - num_cell.k;
    else if (k < 0) m = k + num_cell.k;
    else m = k;
  #endif

#elif PERI_DIM == 0
  t = i; s = j; m = k;
#endif
}
//-------------------------------------------------------
// Return th start & end index of local_tree
// especially for periodical BCs
//-------------------------------------------------------
void Level_info::Calculate_local_tree_start_end(my_int &start, my_int &end)
{
#if PERI_DIM == 0
  start.i = AMAX1(0, cell_start.i);
  start.j = AMAX1(0, cell_start.j);
  start.k = AMAX1(0, cell_start.k);
  end.i   = AMIN1(num_cell.i, cell_end.i);
  end.j   = AMIN1(num_cell.j, cell_end.j);
  end.k   = AMIN1(num_cell.k, cell_end.k);
#elif PERI_DIM != 0
  start.i = 0; start.j = 0; start.k = 0;
  end.i   = 1; end.j   = 1; end.k   = 1;
  #if P_DIM_X == 0 && PERI_DIM_X == 1
    start.i = AMAX1(0, cell_start.i);
    end.i   = AMIN1(num_cell.i, cell_end.i);
  #elif P_DIM_X == 1 && PERI_DIM_X == 1
    start.i = cell_start.i;
    end.i   = cell_end.i;
  #endif

  #if P_DIM_Y == 0 && PERI_DIM_Y == 1
    start.j = AMAX1(0, cell_start.j);
    end.j   = AMIN1(num_cell.j, cell_end.j);
  #elif P_DIM_Y == 1 && PERI_DIM_Y == 1
    start.j = cell_start.j;
    end.j   = cell_end.j;
  #endif

  #if P_DIM_Z == 0 && PERI_DIM_Z == 1
    start.k = AMAX1(0, cell_start.k);
    end.k   = AMIN1(num_cell.k, cell_end.k);
  #elif P_DIM_Z == 1 && PERI_DIM_Z == 1
    start.k = cell_start.k;
    end.k   = cell_end.k;
  #endif
#endif
}
//-------------------------------------------------------
// initialize the exchange tag system for current level
//-------------------------------------------------------
void Level_info::Init_exchange_status()
{
  my_int start, end;

  Calculate_local_tree_start_end(start, end);

  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>
              (start.i, end.i, start.j, end.j, start.k, end.k),
           [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++)
          if (exist_leaf_particle[i][j][k] == 1){
            for( int r=-DIM_X; r<=DIM_X; r++)
              for( int s=-DIM_Y; s<=DIM_Y; s++)
                for( int t=-DIM_Z; t<=DIM_Z; t++){
                  int index_i= AMAX1(cell_start.i,AMIN1(i+r,cell_end.i-1));
                  int index_j= AMAX1(cell_start.j,AMIN1(j+s,cell_end.j-1));
                  int index_k= AMAX1(cell_start.k,AMIN1(k+t,cell_end.k-1));
                  
                  exist_exchange_particle[index_i][index_j][index_k] = 1;
                }
          }
  }, ap);
}
//-------------------------------------------------------
// update the exist tag system for current level
//-------------------------------------------------------
void Level_info::Update_exchange_status(p_Level_info *level_infos, int parent_level)
{
  array_tag parent_exist_exchange_particle = level_infos[parent_level]->exist_exchange_particle;

  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
                [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
          int ii = i/SCALE_RATIO;
          int jj = j/SCALE_RATIO;
          int kk = k/SCALE_RATIO;

          if(parent_exist_exchange_particle[ii][jj][kk] == 1)
            exist_exchange_particle[i][j][k] = 1;
        }
  }, ap);
}
//-------------------------------------------------------
// initialize the exchange tag system for current level
//-------------------------------------------------------
void Level_info::Merge_exchange_status()
{
  my_int start, end;

  Calculate_local_tree_start_end(start, end);

  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>
              (start.i, end.i, start.j, end.j, start.k, end.k),
           [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++)
          if (exist_cell_list[i][j][k] == 1){
            for( int r=-DIM_X; r<=DIM_X; r++)
              for( int s=-DIM_Y; s<=DIM_Y; s++)
                for( int t=-DIM_Z; t<=DIM_Z; t++){
                  int index_i= AMAX1(cell_start.i,AMIN1(i+r,cell_end.i-1));
                  int index_j= AMAX1(cell_start.j,AMIN1(j+s,cell_end.j-1));
                  int index_k= AMAX1(cell_start.k,AMIN1(k+t,cell_end.k-1));
                  
                  exist_exchange_particle[index_i][index_j][index_k] = 1;
                }
          }
  }, ap);
}
//-------------------------------------------------------
// update every level color info for cell list date structure 
// according to the chirldren infomation
//-------------------------------------------------------
void Level_info::Update_every_level_color_info(p_Level_info *level_infos, int child_level)
{
  array_color_list chirldren_color_list       = level_infos[child_level]->color_list;

  int i_c = DIM_X==1 ? SCALE_RATIO : 1;
  int j_c = DIM_Y==1 ? SCALE_RATIO : 1;
  int k_c = DIM_Z==1 ? SCALE_RATIO : 1;

  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(glbl_cell_start.i, glbl_cell_end.i, glbl_cell_start.j, glbl_cell_end.j, glbl_cell_start.k, glbl_cell_end.k),
                [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++){
      for(int j = r.rows().begin(); j < r.rows().end(); j++){
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
          int ii = i*SCALE_RATIO;
          int jj = j*SCALE_RATIO;
          int kk = k*SCALE_RATIO;
          for(int t=0; t<i_c; t++){
            for(int s=0; s<j_c; s++){
              for(int m=0; m<k_c; m++){
                if (chirldren_color_list[ii+t][jj+s][kk+m]->color_list.size() > 0){
                  int num_color_in_chirldren 
                      = int(chirldren_color_list[ii+t][jj+s][kk+m]->color_list.size());
                  for(int icyc=0; icyc<num_color_in_chirldren; icyc++){
                    std:: pair <int,int> pair_t;
                    pair_t = chirldren_color_list[ii+t][jj+s][kk+m]->color_list[icyc];
                    color_list[i][j][k]->Add_color(pair_t);
                  }
                }
              }
            }
          }
        }
      }
    }
  }, ap);
}
#if PERI_DIM != 0
//-------------------------------------------------------
// update every level color info in PBC area
//-------------------------------------------------------
void Level_info::Update_PBC_color_info()
{
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(glbl_cell_start.i, glbl_cell_end.i,
                                     glbl_cell_start.j, glbl_cell_end.j,
                                     glbl_cell_start.k, glbl_cell_end.k),
                [&](const blocked_range3d<int>& r){
    int t, s, m; t = s = m = 0; int tag = 0;
    my_real pos; pos.i = pos.j = pos.k = 0.;
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
#if PERI_DIM_X == 1
  #if P_DIM_X == 1
          if (i < -0.5 || i > num_cell.i-0.5)
            {t = i - i/abs(i)*num_cell.i; tag = 1; pos.i = i/abs(i)*domain.i;}
          else { t = i; pos.i = 0.;}
  #elif P_DIM_X == 0
          t = i; pos.i = 0.;
  #endif
#endif
#if PERI_DIM_Y == 1
  #if P_DIM_Y == 1
          if (j < -0.5 || j > num_cell.j-0.5)
            {s = j - j/abs(j)*num_cell.j; tag = 1; pos.j = j/abs(j)*domain.j;}
          else {s = j; pos.j = 0.;}
  #elif P_DIM_Y == 0
          s = j; pos.j = 0.;
  #endif
#endif
#if PERI_DIM_Z == 1
  #if P_DIM_Z == 1
          if (k < -0.5 || k > num_cell.k-0.5)
            {m = k - k/abs(k)*num_cell.k; tag = 1; pos.k = k/abs(k)*domain.k;}
          else {m = k; pos.k = 0.;}
  #elif P_DIM_Z == 0
          m = k; pos.k = 0.;
  #endif
#endif
          int p_list_size = color_list[t][s][m]->color_list.size();
          if (tag == 1 && p_list_size != 0){
            for(int icyc=0; icyc<p_list_size; icyc++){
              std::pair <int,int> pair_t = color_list[t][s][m]->color_list[icyc];
              color_list[i][j][k]->Add_color(pair_t);
            }
          }
          tag = 0;
        }
  }, ap);
}
//-------------------------------------------------------
// Update the cell_list in current level
//-------------------------------------------------------
void Level_info::Update_cell_list_and_shift_coord(Particle *current_particle, communicator &world)
{
  int i_c = (DIM_X==1 && PERI_DIM_X == 1 && P_DIM_X == 1) ? 1 : 0;
  int j_c = (DIM_Y==1 && PERI_DIM_Y == 1 && P_DIM_Y == 1) ? 1 : 0;
  int k_c = (DIM_Z==1 && PERI_DIM_Z == 1 && P_DIM_Z == 1) ? 1 : 0;
  my_real pos, p_coord;
  my_int shift;
  shift.i = shift.j = shift.k = 0;
  p_coord.i = p_coord.j = p_coord.k = 0.;
  pos.i = current_particle->coord.i;
  pos.j = current_particle->coord.j;
  pos.k = current_particle->coord.k;

  int count = 0;
  for (int rr = -1*i_c; rr <= i_c; rr++){
    p_coord.i = pos.i + rr*domain.i - box_l.i;
    shift.i = DIM_X==1 ? int(floor(p_coord.i/dcell.i)) : 0;
    if (shift.i >= cell_start.i && shift.i < cell_end.i){
      for (int ss = -1*j_c; ss <= j_c; ss++){
        p_coord.j = pos.j + ss*domain.j - box_l.j;
        shift.j = DIM_Y==1 ? int(floor(p_coord.j/dcell.j)) : 0;
        if (shift.j >= cell_start.j && shift.j < cell_end.j){
          for (int tt = -1*k_c; tt <= k_c; tt++){
            p_coord.k = pos.k + tt*domain.k - box_l.k;
            shift.k = DIM_Z==1 ? int(floor(p_coord.k/dcell.k)) : 0;
            if (shift.k >= cell_start.k && shift.k < cell_end.k){
              if (count == 0){
                current_particle->coord = my_add_data (p_coord, box_l);
#ifdef _SCLL_
                int key = 0;
                my_real coord_shift_sub_cell = my_minus_data (p_coord, my_multiply_data(dcell,shift));

                my_int  pos_sub = get_subcell_id (coord_shift_sub_cell, dsubcell, nsubcell);
                key = (pos_sub.k*DIM_Z*nsubcell.j + pos_sub.j*DIM_Y)*nsubcell.i + pos_sub.i*DIM_X;

                std::pair <int,p_Particle> current_particle_pair;
                current_particle_pair.first = key;
                current_particle_pair.second = current_particle;
                table_cell_list[shift.i][shift.j][shift.k]->Add_particle_and_key(current_particle_pair);
#else
                table_cell_list[shift.i][shift.j][shift.k]->Add_particle(current_particle);
#endif
                current_particle->level = level;
                count++;
              }
            }
          }
        }
      }
    }
  }
}
#endif
//-------------------------------------------------------
// traverse the tree and establish global topology graph
//-------------------------------------------------------
void Level_info::Traverse_tree_and_construct_graph(Graphcls *sph_graph)
{
  int i_c = DIM_X==1 ? 2 : 0;
  int j_c = DIM_Y==1 ? 2 : 0;
  int k_c = DIM_Z==1 ? 2 : 0;
  
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(0, num_cell.i, 0, num_cell.j, 0, num_cell.k),
                [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
          int c_list_size = color_list[i][j][k]->color_list.size();
          if( c_list_size> 0){
            for(int icyc=0; icyc<c_list_size; icyc++){
              std::pair <int,int> current_pair = color_list[i][j][k]->color_list[icyc];
              if(current_pair.second == level - Lmin){
                for(int t=AMAX1((i-i_c),glbl_cell_start.i); t<=AMIN1((i+i_c),(glbl_cell_end.i-1)); t++)
                  for(int s=AMAX1((j-j_c),glbl_cell_start.j); s<=AMIN1((j+j_c),(glbl_cell_end.j-1)); s++)
                    for(int m=AMAX1((k-k_c),glbl_cell_start.k); m<=AMIN1((k+k_c),(glbl_cell_end.k-1)); m++){
                      int cj_list_size = color_list[t][s][m]->color_list.size();
                      if (cj_list_size > 0){
                        for(int jcyc=0; jcyc<cj_list_size; jcyc++){
                          std::pair <int,int> neighbor_pair = color_list[t][s][m]->color_list[jcyc];
                          int color_i = current_pair.first;
                          int color_j = neighbor_pair.first;

                          if(color_i != color_j)
                            sph_graph->Add_edge(color_i, color_j);
                        }
                      }
                    }
              }
            }
          }
        }
  }, ap);
}
//-------------------------------------------------------
// traverse the tree and establish global topology graph
//-------------------------------------------------------
void Level_info::Traverse_tree_and_construct_graph_modified(Graphcls *sph_graph)
{
  int i_c = DIM_X==1 ? 2 : 0;
  int j_c = DIM_Y==1 ? 2 : 0;
  int k_c = DIM_Z==1 ? 2 : 0;
  
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(0, num_cell.i, 0, num_cell.j, 0, num_cell.k),
                [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
          int color_list_size_i = color_list[i][j][k]->color_list.size();
          if( color_list_size_i > 0){
            for(int icyc=0; icyc<color_list_size_i; icyc++){
              std::pair <int,int> current_pair = color_list[i][j][k]->color_list[icyc];
              if(current_pair.second == level - Lmin){
                for(int t=AMAX1((i-i_c),glbl_cell_start.i); t<=AMIN1((i+i_c),(glbl_cell_end.i-1)); t++)
                  for(int s=AMAX1((j-j_c),glbl_cell_start.j); s<=AMIN1((j+j_c),(glbl_cell_end.j-1)); s++)
                    for(int m=AMAX1((k-k_c),glbl_cell_start.k); m<=AMIN1((k+k_c),(glbl_cell_end.k-1)); m++){
                      int color_list_size_j = color_list[t][s][m]->color_list.size();
                      if ( color_list_size_j > 0){
                        for(int jcyc=0; jcyc<color_list_size_j; jcyc++){
                          std::pair <int,int> neighbor_pair = color_list[t][s][m]->color_list[jcyc];
                          int color_i = current_pair.first;
                          int color_j = neighbor_pair.first;

                          if(color_i != color_j){
                            sph_graph->graph_matrix[color_i][color_j] = 1;
                          }
                        }
                      }
                    }
              }
            }
          }
        }
  }, ap);
}
//-------------------------------------------------------
// exchange the graph topology information
//-------------------------------------------------------
void Level_info::Exchange_graph_topology(communicator &world, serialization_vector <unsigned int> &exchange_vector)
{
  exchange_vector.Vector.push_back((unsigned int)cell_start.i);
  exchange_vector.Vector.push_back((unsigned int)cell_start.j);
  exchange_vector.Vector.push_back((unsigned int)cell_start.k);
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
                [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++){
#ifndef _NARROW_BAND_GRAPH_
          if(exist_leaf_particle[i][j][k] == 1){
#else
          if(exist_leaf_particle[i][j][k] == 1 && list_is_narrow_band[i][j][k] == 1){
#endif            
            unsigned int data = compress(i-cell_start.i, j-cell_start.j, k-cell_start.k, 0);
            if (data == 0xffffffff) {
              std::cout<< "compress error" << endl; world.abort(-1);
            }
            exchange_vector.Vector.push_back(data);
          }
        }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Exchange_graph_topology in "<<level<<" is finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// exchange the buffer particle between pairs
//-------------------------------------------------------
void Level_info::Exchange_buffer_between_pairs(communicator &world,int target_processor, concurrent_vector  <p_Particle> &particle_ghost, int flag)
{
  if(target_processor != -1){

    serialization_vector <unsigned int> exchange_vector_in;
    serialization_vector <unsigned int> exchange_vector_out;
    exchange_vector_in.Vector.clear();
    exchange_vector_out.Vector.clear();

    // reset matrix
    static affinity_partitioner ap;
    parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
             [&](const blocked_range3d<int>& r){
      for(int i = r.pages().begin(); i < r.pages().end(); i++)
        for(int j = r.rows().begin(); j < r.rows().end(); j++)
          for(int k = r.cols().begin(); k < r.cols().end(); k++)
              exchange_exist_cell_list[i][j][k] = 0;
    }, ap);

    // fill the output data structure
    exchange_vector_out.Vector.push_back((unsigned int)cell_start.i);
    exchange_vector_out.Vector.push_back((unsigned int)cell_start.j);
    exchange_vector_out.Vector.push_back((unsigned int)cell_start.k);

    parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
             [&](const blocked_range3d<int>& r){
      for(int i = r.pages().begin(); i < r.pages().end(); i++)
        for(int j = r.rows().begin(); j < r.rows().end(); j++)
          for(int k = r.cols().begin(); k < r.cols().end(); k++)
            if ( exist_exchange_particle[i][j][k] == 1 ){
              unsigned int data = compress(i-cell_start.i,j-cell_start.j,k-cell_start.k,0);
              if (data == 0xffffffff) {
                std::cout<< "compress error" << endl; world.abort(-1);
              }
              exchange_vector_out.Vector.push_back(data);
            }
    }, ap);
    
    exchange_vector_out.mem_size = int(exchange_vector_out.Vector.size());
    exchange_vector_out.tag = 1;

    // nonblocking communication
    mpi::request reqs[2];
    reqs[0] = world.isend(target_processor, level-Lmin, exchange_vector_out);
    reqs[1] = world.irecv(target_processor, level-Lmin, exchange_vector_in);
    mpi::wait_all(reqs, reqs + 2);

    my_int st;
    my_set_const(st, 0);
    st.i = int(exchange_vector_in.Vector[0]);
    st.j = int(exchange_vector_in.Vector[1]);
    st.k = int(exchange_vector_in.Vector[2]);

#if PERI_DIM != 0
    int i_c = (DIM_X==1 && PERI_DIM_X == 1 && P_DIM_X == 1) ? 1 : 0;
    int j_c = (DIM_Y==1 && PERI_DIM_Y == 1 && P_DIM_Y == 1) ? 1 : 0;
    int k_c = (DIM_Z==1 && PERI_DIM_Z == 1 && P_DIM_Z == 1) ? 1 : 0;
#endif

    parallel_for( blocked_range<int>(3, int(exchange_vector_in.Vector.size())),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        int i_t, j_t, k_t, level_temp;
        unsigned int data;
        data = exchange_vector_in.Vector[i];

        if (data == 0xffffffff) {
          std::cout<< "compress error" << endl; world.abort(-1);
        }
        // uncompress
        uncompress(data, i_t, j_t, k_t, level_temp);
	i_t += st.i;
	j_t += st.j;
	k_t += st.k;

        int i_tt, j_tt, k_tt; i_tt = j_tt = k_tt = 0;
        Shift_cell_index(i_tt, j_tt, k_tt, i_t, j_t, k_t);
#if PERI_DIM == 0
	if ( i_tt >= cell_start.i && i_tt < cell_end.i
          && j_tt >= cell_start.j && j_tt < cell_end.j
          && k_tt >= cell_start.k && k_tt < cell_end.k){
          exchange_exist_cell_list[i_tt][j_tt][k_tt] = 1;
        }
#elif PERI_DIM != 0
        for (int rr = -1*i_c*num_cell.i; rr <= i_c*num_cell.i; rr+=num_cell.i){
          int id_i = i_tt+rr;
          if (id_i >= cell_start.i && id_i < cell_end.i){
            for (int ss = -1*j_c*num_cell.j; ss <= j_c*num_cell.j; ss+=num_cell.j){
              int id_j = j_tt+ss;
              if (id_j >= cell_start.j && id_j < cell_end.j){
                for (int tt = -1*k_c*num_cell.k; tt <= k_c*num_cell.k; tt+=num_cell.k){
                  int id_k = k_tt+tt;
                  if (id_k >= cell_start.k && id_k < cell_end.k){
                    exchange_exist_cell_list[id_i][id_j][id_k] = 1;
                  }
                }
              }
            }
          }
        }
#endif
      }
    }, ap);

    // prepare the data for transfering

    parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
             [&](const blocked_range3d<int>& r){
      for(int i = r.pages().begin(); i < r.pages().end(); i++)
        for(int j = r.rows().begin(); j < r.rows().end(); j++)
          for(int k = r.cols().begin(); k < r.cols().end(); k++)
            if ( exist_leaf_particle[i][j][k] == 1 && exchange_exist_cell_list[i][j][k] == 1){
              int p_list_size = table_cell_list[i][j][k]->particle_list.size();
              for(int icyc=0; icyc<p_list_size; icyc++){
#ifndef _SCLL_
                Particle *current_particle = table_cell_list[i][j][k]->particle_list[icyc];
#else
                Particle *current_particle = table_cell_list[i][j][k]->particle_list[icyc].second;
#endif
                if(current_particle->level == level && current_particle->color == world.rank()){
                  current_particle->tag = flag;
                  particle_ghost.push_back(current_particle);
                }
              }
            }
    }, ap);
    exchange_vector_in.Vector.clear();
    exchange_vector_out.Vector.clear();
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Exchange_buffer_between"<<world.rank()<<" & "<<target_processor<<" in level "<<level<<" is finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// output current level 
//--------------------------------------------------
void Level_info::Output_color_list(communicator &world, char *filename, int n)
{
  ofstream out(filename, ios::app);

  //defining header for tecplot(plot software)
  out<<"\n";
  out<<"title='View'"<<"\n";
  out<<"variables=";
#if DIM_X
  out<<"x, ";
#endif
#if DIM_Y
  out<<"y, ";
#endif
#if DIM_Z
  out<<"z, ";
#endif

  int LEN;
  LEN = 3;
  out<<"level, color, np"<<"\n";

  my_int start;
  my_int end;
 
  start.i = glbl_cell_start.i;
  start.j = glbl_cell_start.j;
  start.k = glbl_cell_start.k;
  end.i   = glbl_cell_end.i;
  end.j   = glbl_cell_end.j;
  end.k   = glbl_cell_end.k;

  out<<"zone t='sub_filed_"<<n<<"_"<<world.rank()<<"_"<<level-Lmin<<"'"
     <<"  i= "<<end.i-start.i+1
     <<"  j= "<<end.j-start.j+1
     <<"  k= "<<end.k-start.k+1
     <<"  DATAPACKING=BLOCK, VARLOCATION=([";
  int pos_s = DIM_X+DIM_Y+DIM_Z+1;
  out<<pos_s<<"-";
  out<<2*pos_s -1 + LEN-1<<"]=CELLCENTERED) SOLUTIONTIME="<<n<<"\n";

#if DIM_X
  for(int k=start.k; k<=end.k; k++){
    for(int j=start.j; j<=end.j; j++){
      for(int i=start.i; i<=end.i; i++){
        out<<dcell.i*i+box_l.i<<" ";
      }
      out<<"\n";
    }
  }
#endif
#if DIM_Y
  for(int k=start.k; k<=end.k; k++){
    for(int j=start.j; j<=end.j; j++){
      for(int i=start.i; i<=end.i; i++){
        out<<dcell.j*j+box_l.j<<" ";
      }
      out<<"\n";
    }
  }
#endif
#if DIM_Z
  for(int k=start.k; k<=end.k; k++){
    for(int j=start.j; j<=end.j; j++){
      for(int i=start.i; i<=end.i; i++){
        out<<dcell.k*k+box_l.k<<" ";
      }
      out<<"\n";
    }
  }
#endif

  // level
  for(int k=start.k; k<end.k; k++){
    for(int j=start.j; j<end.j; j++){
      for(int i=start.i; i<end.i; i++){
        out<<level<<" ";
      }
      out<<"\n";
    }
  }

  // color
  for(int k=start.k; k<end.k; k++){
    for(int j=start.j; j<end.j; j++){
      for(int i=start.i; i<end.i; i++){
        if (int(color_list[i][j][k]->color_list.size()) > 0){
          int m = 0;
          for (int l=0; l < color_list[i][j][k]->color_list.size(); l ++){
            if (color_list[i][j][k]->color_list[0].first != color_list[i][j][k]->color_list[l].first) m = 1;
          }
          if (m == 0)
            out<<color_list[i][j][k]->color_list[0].first<<" ";
          else if (m == 1)
            out<<"999 ";
        }else out<<"-1 ";
      }
      out<<"\n";
    }
  }

  // np
  for(int k=start.k; k<end.k; k++){
    for(int j=start.j; j<end.j; j++){
      for(int i=start.i; i<end.i; i++){
        out<<color_list[i][j][k]->color_list.size()<<" ";
      }
      out<<"\n";
    }
  }
  out.close();

}
#endif
#ifdef _MEM_CHECK_
//-------------------------------------------------------
// calculate memory consumption in current level
// NOTE: if you are using concurrent_vector, please be 
// aware that push_back operation in function
// Refresh_neighbor_info will allocate large chunk of
// memory, with cannot be accounted in current function.
// Just be carefull of using safeguard factor!
//-------------------------------------------------------
void Level_info::Check_memory_consumption(communicator &world, long long &mem_cell_temp, long long &mem_list_temp)
{
  long long count1 = 0;
#ifdef _MPI_
  long long count2 = 0;
#endif
  for (int i = cell_start.i; i < cell_end.i; i++){
    for (int j = cell_start.j; j < cell_end.j; j++){
      for (int k = cell_start.k; k < cell_end.k; k++){
        count1 += (long long)(table_cell_list[i][j][k]->particle_list.capacity());
#ifdef _MPI_
        count2 += (long long)(color_list[i][j][k]->color_list.capacity());
#endif
      }
    }
  }
  mem_cell_temp += (long long)(count1*sizeof(p_Particle));
#ifdef _MPI_
  mem_cell_temp += (long long)(count2*sizeof(std::pair <int,int>));
#endif
  long long num_cell = (long long)(cell_end.i-cell_start.i)*
                      (long long)(cell_end.j-cell_start.j)*
                      (long long)(cell_end.k-cell_start.k);

  mem_list_temp += (long long)(num_cell*sizeof(DTAG)*2);
#ifdef _MPI_
  mem_list_temp += (long long)(num_cell*sizeof(DTAG)*2);
#endif
}
#endif
