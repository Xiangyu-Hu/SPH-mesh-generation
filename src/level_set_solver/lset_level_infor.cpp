#include "glbfunc.h"
#include "level_set.h"
#include "lset_level_infor.h"

typedef boost::multi_array_types::extent_range mrange;

/***************************************************/
/*                                                 */
/* Functions defined in class "Levelset_levelinfo" */
/*                                                 */
/***************************************************/

//----------------------------------------------------------------
// initialze all the necessary parameters for 
// Levelset_levelinfo
//---------------------------------------------------------------
void Levelset_levelinfo::Initialize(int i, Levelset *level_set, communicator &world)
{
  level = i;
  Lmin = level_set->Lmin;
  Lmax = level_set->Lmax;
  
  p_lset_pkg_pool = &(level_set->lset_pkg_pool);
  
  Real multi = powern(Real(SCALE_RATIO),Lmax-Lmin);
  
  my_set_data  (domain, level_set->domain);
  my_set_data  (box_l, level_set->box_l);
  my_set_data  (box_r, level_set->box_r);
  my_set_const (num_pkg, 1);
  my_set_data  (num_pkg, my_multiply_const(level_set->ini_num_cell, multi));
  my_self_multiply (num_pkg, glbl_total_num_pkg);
  dpkg = my_devide_data (domain, num_pkg);
  dcell.i = dpkg.i / Real (ICPX);
  dcell.j = dpkg.j / Real (ICPY);
  dcell.k = dpkg.k / Real (ICPZ);
  
#if DIM == 2
  if (fabs(dpkg.i - dpkg.j) >= 1.e-12){
    if (world.rank() == 0) cout<<"ERROR: The package size is wrong!!!!"<<endl;
  }
#elif DIM == 3
  if (fabs(dpkg.i - dpkg.k) >= 1.e-12){
    if (world.rank() == 0) cout<<"ERROR: The package size is wrong!!!!"<<endl;
  }
#endif

  scale = dpkg.i;
  dl = 1.e20;
  if (DIM_X) dl = AMIN1(dcell.i, dl);
  if (DIM_Y) dl = AMIN1(dcell.j, dl);
  if (DIM_Z) dl = AMIN1(dcell.k, dl);

  pkg_start.i = pkg_start.j = pkg_start.k = 0;
  pkg_end.i = pkg_end.j = pkg_end.k = 1;
  my_set_data (pkg_end, num_pkg);

  int dL = Lmax - Lmin;
  int dN = powern(SCALE_RATIO,dL);
#if DIM_X
  pkg_start.i -= dN;
  pkg_end.i += dN;
  num_pkg_margin.i = dN;
#endif
#if DIM_Y
  pkg_start.j -= dN;
  pkg_end.j += dN;
  num_pkg_margin.j = dN;
#endif
#if DIM_Z
  pkg_start.k -= dN;
  pkg_end.k += dN;
  num_pkg_margin.k = dN;
#endif
  glbl_pkg_start.i = pkg_start.i;
  glbl_pkg_start.j = pkg_start.j;
  glbl_pkg_start.k = pkg_start.k;
  glbl_pkg_end.i = pkg_end.i;
  glbl_pkg_end.j = pkg_end.j;
  glbl_pkg_end.k = pkg_end.k;
  my_set_const (glbl_num_pkg, 1);
  glbl_num_pkg.i = glbl_pkg_end.i - glbl_pkg_start.i;
  glbl_num_pkg.j = glbl_pkg_end.j - glbl_pkg_start.j;
  glbl_num_pkg.k = glbl_pkg_end.k - glbl_pkg_start.k;
  
  my_self_multiply (glbl_num_pkg, glbl_total_num_pkg);
  
  Allocate_memory();

  Initial_pkg_list_and_reset_tags(level_set);
  
  #if !defined(_READ_TARGET_FIELD_) && !defined(_READ_TARGET_FIELD_UNCHANGE_) && !defined(_READ_TARGET_FIELD_FOR_POST_)
    #ifdef _READ_SDF_
    Load_SDF_file(level_set, world);
    #else
    Define_levelset_full_field(level_set);
    #endif
  #endif
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Level set Level     : "<<level<<" is initialized\n";
    cout<<"<<<<< Total number of pkg : "<<right<<setw(10)<<glbl_total_num_pkg<<"\n";
    cout<<"<<<<< Scale size          : "<<right<<setw(10)<<scale             <<"\n";
    cout<<"<<<<< Cell size           : "<<right<<setw(10)<<dcell.i           <<" | "<<right<<setw(10)<<dcell.j         <<" | "<<right<<setw(10)<<dcell.k<<"\n";
    cout<<"<<<<< pkg start           : "<<right<<setw(10)<<glbl_pkg_start.i  <<" | "<<right<<setw(10)<<glbl_pkg_start.j<<" | "<<right<<setw(10)<<glbl_pkg_start.k<<"\n";
    cout<<"<<<<< pkg end             : "<<right<<setw(10)<<glbl_pkg_end.i    <<" | "<<right<<setw(10)<<glbl_pkg_end.j  <<" | "<<right<<setw(10)<<glbl_pkg_end.k<<"\n";
    cout<<"**********************************************************\n";
  }
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Levelset level "<<level<<" is initialized\n";
  out.close();
#endif
}
//-------------------------------------------------------
// Set_cell_topology
//-------------------------------------------------------
void Levelset_levelinfo::Set_cell_topology (Levelset *level_set, communicator &world)
{
  static affinity_partitioner ap;

  parallel_for( blocked_range3d<int>(0, num_pkg.i, 0, num_pkg.j, 0, num_pkg.k),
                [&](const blocked_range3d<int>& r){
    for(int i=r.pages().begin(); i!=r.pages().end(); ++i){
      for(int j=r.rows().begin(); j!=r.rows().end(); ++j){
        for(int k=r.cols().begin(); k!=r.cols().end(); ++k){
          table_lset_pkg_list[i][j][k]->Set_cell_topology(this);
        }
      }
    }
  }, ap);
}
//-------------------------------------------------------
// Find_interface_packages
//-------------------------------------------------------
void Levelset_levelinfo::Find_interface_packages (Levelset *level_set, communicator &world)
{
  static affinity_partitioner ap;
  tbb::mutex gMutex;
  parallel_for( blocked_range3d<int>(0, num_pkg.i, 0, num_pkg.j, 0, num_pkg.k),
                [&](const blocked_range3d<int>& r){
    for(int i=r.pages().begin(); i!=r.pages().end(); ++i){
      for(int j=r.rows().begin(); j!=r.rows().end(); ++j){
        for(int k=r.cols().begin(); k!=r.cols().end(); ++k){
          if (table_lset_pkg_list[i][j][k]->Return_package_interface_tag(this) > 0){
            tbb::mutex::scoped_lock lock(gMutex);
            level_set->interface_pkg.push_back (table_lset_pkg_list[i][j][k]);
          }
        }
      }
    }
  }, ap);
}
//-------------------------------------------------------
// Get_extended_cell_tags
//-------------------------------------------------------
void Levelset_levelinfo::Get_extended_cell_tags (Levelset *level_set, communicator &world)
{
  std::vector<std::pair <p_Levelset_cell, my_real>> seg_cells; seg_cells.clear();
  std::vector<std::pair <p_Levelset_cell, my_real>> sing_cells; sing_cells.clear();
  static affinity_partitioner ap;
  tbb::mutex eMutex;
  tbb::mutex iMutex;
  parallel_for( blocked_range3d<int>(0, num_pkg.i, 0, num_pkg.j, 0, num_pkg.k),
                [&](const blocked_range3d<int>& r){
    for(int i=r.pages().begin(); i!=r.pages().end(); ++i){
      for(int j=r.rows().begin(); j!=r.rows().end(); ++j){
        for(int k=r.cols().begin(); k!=r.cols().end(); ++k){
          p_Levelset_package pkg_i = table_lset_pkg_list[i][j][k];
          my_int index = pkg_i->index;
          for (int i = 0; i < ICPX; i++)
            for (int j = 0; j < ICPY; j++)
              for (int k = 0; k < ICPZ; k++){
                p_Levelset_cell cell_i = pkg_i->p_cell[i][j][k];
                if (cell_i->tag_characteristic == SINGULARITY_CELL){
                  std::pair <p_Levelset_cell, my_real> tmp;
                  my_real coord = pkg_i->get_cell_position (index.i, index.j, index.k, i, j, k, dpkg, dcell, box_l);
                  tmp.first  = cell_i;
                  tmp.second = coord;
                  tbb::mutex::scoped_lock lock(iMutex);
                  sing_cells.push_back (tmp);
                } else if (cell_i->tag_characteristic == SEGMENT_CELL ){
                  std::pair <p_Levelset_cell, my_real> tmp;
                  my_real coord = pkg_i->get_cell_position (index.i, index.j, index.k, i, j, k, dpkg, dcell, box_l);
                  tmp.first  = cell_i;
                  tmp.second = coord;                  
                  tbb::mutex::scoped_lock lock(eMutex);
                  seg_cells.push_back (tmp);
                }
              }
        }
      }
    }
  }, ap);

  int seg_cell_size = seg_cells.size();
  for (int i = 0; i < seg_cell_size; ++i)
  {
    p_Levelset_cell cell_i = seg_cells[i].first;

    Real rad = 0.5*level_set->Get_scale(cell_i->curv, cell_i->phi);
    my_real center = seg_cells[i].second;
    my_int id_pkg = get_cell_id (center, dpkg, glbl_pkg_start, glbl_pkg_end);

    p_Levelset_package pkg_i = table_lset_pkg_list[id_pkg.i][id_pkg.j][id_pkg.k];

    if (cell_i->tag_characteristic != SEGMENT_CELL){
      cout<<"<<<<< ERROR!!! in SEGMENT_CELL for Get_extended_cell_tags"<<endl;
    }
    int npkg = int(ceil(rad/scale));

    int im = DIM_X==1 ? -npkg   : 0;
    int ip = DIM_X==1 ?  npkg+1 : 1;
    int jm = DIM_Y==1 ? -npkg   : 0;
    int jp = DIM_Y==1 ?  npkg+1 : 1;
    int km = DIM_Z==1 ? -npkg   : 0;
    int kp = DIM_Z==1 ?  npkg+1 : 1;

    for (int i = im; i < ip; ++i)
      for (int j = jm; j < jp; ++j)
        for (int k = km; k < kp; ++k)
        {
          my_int pkg_j;
          pkg_j.i = AMAX1 (AMIN1(id_pkg.i + i, num_pkg.i-1), 0);
          pkg_j.j = AMAX1 (AMIN1(id_pkg.j + j, num_pkg.j-1), 0);
          pkg_j.k = AMAX1 (AMIN1(id_pkg.k + k, num_pkg.k-1), 0);

          p_Levelset_package current_pkg = table_lset_pkg_list[pkg_j.i][pkg_j.j][pkg_j.k];

          for (int ii = 0; ii < ICPX; ii++)
            for (int jj = 0; jj < ICPY; jj++)
              for (int kk = 0; kk < ICPZ; kk++){
                p_Levelset_cell current_cell = current_pkg->p_cell[ii][jj][kk];
                  my_real cell_center = current_pkg->get_cell_position (pkg_j.i, pkg_j.j, pkg_j.k, ii, jj, kk, dpkg, dcell, box_l);
                  Real dist = get_distance_2p (center, cell_center);

                  if (dist < rad){
                    if (current_cell->tag_characteristic == NORMAL_CELL){
                      current_cell->tag_characteristic = SEGMENT_CELL;
                      current_cell->idx_characteristic = cell_i->idx_characteristic;
                    }
                  }
              }
        }
  }

  int sing_cell_size = sing_cells.size();
  for (int i = 0; i < sing_cell_size; ++i)
  {
    p_Levelset_cell cell_i = sing_cells[i].first;

    Real rad = 0.5*level_set->Get_scale(cell_i->curv, cell_i->phi);
    my_real center = sing_cells[i].second;
    my_int id_pkg = get_cell_id (center, dpkg, glbl_pkg_start, glbl_pkg_end);

    p_Levelset_package pkg_i = table_lset_pkg_list[id_pkg.i][id_pkg.j][id_pkg.k];

    if (cell_i->tag_characteristic != SINGULARITY_CELL){
      cout<<"<<<<< ERROR!!! in SINGULARITY_CELL for Get_extended_cell_tags"<<endl;
    }
    int npkg = int(ceil(rad/scale));

    int im = DIM_X==1 ? -npkg   : 0;
    int ip = DIM_X==1 ?  npkg+1 : 1;
    int jm = DIM_Y==1 ? -npkg   : 0;
    int jp = DIM_Y==1 ?  npkg+1 : 1;
    int km = DIM_Z==1 ? -npkg   : 0;
    int kp = DIM_Z==1 ?  npkg+1 : 1;

    for (int i = im; i < ip; ++i)
      for (int j = jm; j < jp; ++j)
        for (int k = km; k < kp; ++k)
        {
          my_int pkg_j;
          pkg_j.i = AMAX1 (AMIN1(id_pkg.i + i, num_pkg.i-1), 0);
          pkg_j.j = AMAX1 (AMIN1(id_pkg.j + j, num_pkg.j-1), 0);
          pkg_j.k = AMAX1 (AMIN1(id_pkg.k + k, num_pkg.k-1), 0);

          p_Levelset_package current_pkg = table_lset_pkg_list[pkg_j.i][pkg_j.j][pkg_j.k];

          for (int ii = 0; ii < ICPX; ii++)
            for (int jj = 0; jj < ICPY; jj++)
              for (int kk = 0; kk < ICPZ; kk++){
                p_Levelset_cell current_cell = current_pkg->p_cell[ii][jj][kk];
                  my_real cell_center = current_pkg->get_cell_position (pkg_j.i, pkg_j.j, pkg_j.k, ii, jj, kk, dpkg, dcell, box_l);
                  Real dist = get_distance_2p (center, cell_center);

                  if (dist < rad){
                    if (current_cell->tag_characteristic != SINGULARITY_CELL){
                      current_cell->tag_characteristic = SINGULARITY_CELL;
                      current_cell->idx_characteristic = cell_i->idx_characteristic;
                    }
                  }
              }
        }
  }
}
//-------------------------------------------------------
// Allocate_memory
//-------------------------------------------------------
void Levelset_levelinfo::Allocate_memory()
{
  table_lset_pkg_list.resize
    (boost::extents[mrange(pkg_start.i,pkg_end.i)][mrange(pkg_start.j,pkg_end.j)][mrange(pkg_start.k,pkg_end.k)]);
}
//-------------------------------------------------------
// reset tags to 0 for initialization
//-------------------------------------------------------
void Levelset_levelinfo::Initial_pkg_list_and_reset_tags(Levelset *level_set)
{
  for(int i=pkg_start.i; i!=pkg_end.i; ++i){
    for(int j=pkg_start.j; j!=pkg_end.j; ++j){
      for(int k=pkg_start.k; k!=pkg_end.k; ++k){
        if (NULL == table_lset_pkg_list[i][j][k]){
          table_lset_pkg_list[i][j][k]   = p_lset_pkg_pool->malloc();
          if (i>=0 && i<num_pkg.i &&
              j>=0 && j<num_pkg.j &&
              k>=0 && k<num_pkg.k) level_set->lset_pkg.push_back(table_lset_pkg_list[i][j][k]);
          else level_set->lset_bc_pkg.push_back(table_lset_pkg_list[i][j][k]);
        }
        table_lset_pkg_list[i][j][k]->Initialize(this, i, j, k);
      }
    }
  }
}
//-------------------------------------------------------
// Define_levelset_full_field
//-------------------------------------------------------
void Levelset_levelinfo::Define_levelset_full_field(Levelset *level_set)
{
  static affinity_partitioner ap;

  parallel_for( blocked_range3d<int>(pkg_start.i, pkg_end.i, pkg_start.j, pkg_end.j, pkg_start.k, pkg_end.k),
                [&](const blocked_range3d<int>& r){
    for(int i=r.pages().begin(); i!=r.pages().end(); ++i){
      for(int j=r.rows().begin(); j!=r.rows().end(); ++j){
        for(int k=r.cols().begin(); k!=r.cols().end(); ++k){
          table_lset_pkg_list[i][j][k]->Define_levelset_full_field(i, j, k, dpkg, dcell, level_set);
        }
      }
    }
  }, ap);
}
//-------------------------------------------------------
// Get_id_pkg_cell
//-------------------------------------------------------
void Levelset_levelinfo::Get_id_pkg_cell(my_real pos, my_int &id_pkg, my_int &id_cell)
{
  
  id_pkg = get_cell_id (pos, dpkg, glbl_pkg_start, glbl_pkg_end);
  
  id_cell.i = id_cell.j = id_cell.k = 0.;
  
  my_real pos_shift_pkg = my_minus_data (pos, my_multiply_data (dpkg, id_pkg));
  
  id_cell.i = DIM_X==1 ? int(floor(pos_shift_pkg.i/dcell.i)) : 0;
  id_cell.j = DIM_Y==1 ? int(floor(pos_shift_pkg.j/dcell.j)) : 0;
  id_cell.k = DIM_Z==1 ? int(floor(pos_shift_pkg.k/dcell.k)) : 0;

  id_cell.i = AMAX1(0,AMIN1(id_cell.i,ICPX-1));
  id_cell.j = AMAX1(0,AMIN1(id_cell.j,ICPY-1));
  id_cell.k = AMAX1(0,AMIN1(id_cell.k,ICPZ-1));
}
//-------------------------------------------------------
// Search_for_char_cell_within_dx
//-------------------------------------------------------
void Levelset_levelinfo::Search_for_char_cell_within_dx(my_real coord, Real rad, int &tag_interface, int &tag_characteristic, int &idx_characteristic)
{
  
  my_int id_pkg = get_cell_id (coord, dpkg, glbl_pkg_start, glbl_pkg_end);
  
  my_int id_cell; id_cell.i = id_cell.j = id_cell.k = 0.;
  
  my_real pos_shift_pkg = my_minus_data (coord, my_multiply_data (dpkg, id_pkg));
  
  id_cell.i = DIM_X==1 ? int(floor(pos_shift_pkg.i/dcell.i)) : 0;
  id_cell.j = DIM_Y==1 ? int(floor(pos_shift_pkg.j/dcell.j)) : 0;
  id_cell.k = DIM_Z==1 ? int(floor(pos_shift_pkg.k/dcell.k)) : 0;

  id_cell.i = AMAX1(0,AMIN1(id_cell.i,ICPX-1));
  id_cell.j = AMAX1(0,AMIN1(id_cell.j,ICPY-1));
  id_cell.k = AMAX1(0,AMIN1(id_cell.k,ICPZ-1));

  p_Levelset_package pkg_i = table_lset_pkg_list[id_pkg.i][id_pkg.j][id_pkg.k];

  if (pkg_i->p_cell[id_cell.i][id_cell.j][id_cell.k]->tag_interface != NORMAL_CELL){

    if (pkg_i->p_cell[id_cell.i][id_cell.j][id_cell.k]->tag_interface == CUT_CELL)
      tag_interface = CUT_CELL;

    int npkg = int(ceil(rad/scale));

    int im = DIM_X==1 ? -npkg   : 0;
    int ip = DIM_X==1 ?  npkg+1 : 1;
    int jm = DIM_Y==1 ? -npkg   : 0;
    int jp = DIM_Y==1 ?  npkg+1 : 1;
    int km = DIM_Z==1 ? -npkg   : 0;
    int kp = DIM_Z==1 ?  npkg+1 : 1;

    for (int i = im; i < ip; ++i)
      for (int j = jm; j < jp; ++j)
        for (int k = km; k < kp; ++k)
        {
          my_int pkg_j;
          pkg_j.i = AMAX1 (AMIN1(id_pkg.i + i, pkg_end.i), pkg_start.i);
          pkg_j.j = AMAX1 (AMIN1(id_pkg.j + j, pkg_end.j), pkg_start.j);
          pkg_j.k = AMAX1 (AMIN1(id_pkg.k + k, pkg_end.k), pkg_start.k);

          p_Levelset_package current_pkg = table_lset_pkg_list[pkg_j.i][pkg_j.j][pkg_j.k];

          for (int ii = 0; ii < ICPX; ii++)
            for (int jj = 0; jj < ICPY; jj++)
              for (int kk = 0; kk < ICPZ; kk++){
                p_Levelset_cell current_cell = current_pkg->p_cell[ii][jj][kk];
                if (current_cell->tag_interface != NORMAL_CELL){
                  my_real cell_center = current_pkg->get_cell_position (pkg_j.i, pkg_j.j, pkg_j.k, ii, jj, kk, dpkg, dcell, box_l);
                  Real dist = get_distance_2p (coord, cell_center);

                  if (dist < rad){
                    if (current_cell->tag_characteristic == SEGMENT_CELL && tag_characteristic != SINGULARITY_CELL){
                      tag_characteristic = SEGMENT_CELL;
                      idx_characteristic = current_cell->idx_characteristic;
                    }else if (current_cell->tag_characteristic == SINGULARITY_CELL){
                      tag_characteristic = SINGULARITY_CELL;
                      idx_characteristic = current_cell->idx_characteristic;
                      break;
                    }
                  }
                }
              }
        }
  }
}
//-------------------------------------------------------
// Get_unique_cell_id
//-------------------------------------------------------
int Levelset_levelinfo::Get_unique_cell_id(my_int id_pkg, my_int id_cell)
{
  int unique_index = 0;
  
  my_int idd_cell;
  
  idd_cell.i = DIM_X == 1 ? ICPX*id_pkg.i + id_cell.i : 0;
  idd_cell.j = DIM_Y == 1 ? ICPY*id_pkg.j + id_cell.j : 0;
  idd_cell.k = DIM_Z == 1 ? ICPZ*id_pkg.k + id_cell.k : 0;
  
  unique_index = (idd_cell.i*glbl_num_pkg.j*ICPY + idd_cell.j)*glbl_num_pkg.k*ICPZ + idd_cell.k;
  
  return unique_index;
}
//------------------------------------------------------------
// Output_level_set
//------------------------------------------------------------
void Levelset_levelinfo::Output_level_set(Levelset *level_set, communicator &world, char *filename, int n)
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
  LEN = 11;
  out<<"phi, n_x, n_y, n_z, curv, psi, interface, characteristic, h, vol, area"<<"\n";
  
  out<<"zone t='sub_filed_levelset'  i= "<<num_pkg.i*ICPX+DIM_X<<"  j= "<<num_pkg.j*ICPY+DIM_Y<<"  k= "<<num_pkg.k*ICPZ+DIM_Z<<"  DATAPACKING=BLOCK, VARLOCATION=([";
  int pos_s = DIM_X+DIM_Y+DIM_Z+1;
  out<<pos_s<<"-";
  out<<2*pos_s -1 + LEN-1<<"]=CELLCENTERED) SOLUTIONTIME="<<n<<"\n";

#if DIM_X
  for(int k=AMAX1(pkg_start.k, 0); k<=AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ+DIM_Z-1; k++){
    for(int j=AMAX1(pkg_start.j, 0); j<=AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY+DIM_Y-1; j++){
      for(int i=AMAX1(pkg_start.i, 0); i<=AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX+DIM_X-1; i++){
        out<<dcell.i*i+box_l.i<<" ";
      }
      out<<"\n";
    }
  }
#endif
#if DIM_Y
  for(int k=AMAX1(pkg_start.k, 0); k<=AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ+DIM_Z-1; k++){
    for(int j=AMAX1(pkg_start.j, 0); j<=AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY+DIM_Y-1; j++){
      for(int i=AMAX1(pkg_start.i, 0); i<=AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX+DIM_X-1; i++){
        out<<dcell.j*j+box_l.j<<" ";
      }
      out<<"\n";
    }
  }
#endif
#if DIM_Z
  for(int k=AMAX1(pkg_start.k, 0); k<=AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ+DIM_Z-1; k++){
    for(int j=AMAX1(pkg_start.j, 0); j<=AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY+DIM_Y-1; j++){
      for(int i=AMAX1(pkg_start.i, 0); i<=AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX+DIM_X-1; i++){
        out<<dcell.k*k+box_l.k<<" ";
      }
      out<<"\n";
    }
  }
#endif

  // phi
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].phi<<" ";
      }
      out<<"\n";
    }
  }

  // n_x
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_x<<" ";
      }
      out<<"\n";
    }
  }

  // n_y
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_y<<" ";
      }
      out<<"\n";
    }
  }

  // n_z
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_z<<" ";
      }
      out<<"\n";
    }
  }

  // curv
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].curv<<" ";
      }
      out<<"\n";
    }
  }

  // psi
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].psi<<" ";
      }
      out<<"\n";
    }
  }

  // interface
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        out<<int(table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].tag_interface)<<" ";
      }
      out<<"\n";
    }
  }
  
  // interface
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        out<<int(table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].tag_characteristic)<<" ";
      }
      out<<"\n";
    }
  }
  
  // scale;
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].scale<<" ";
      }
      out<<"\n";
    }
  }

  // volume;
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].vol<<" ";
      }
      out<<"\n";
    }
  }

  // area;
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].area<<" ";
      }
      out<<"\n";
    }
  }
  out.close();
}
//------------------------------------------------------------
// Write_level_set_rstfile
//------------------------------------------------------------
void Levelset_levelinfo::Write_level_set_rstfile(Levelset *level_set, communicator &world, char *filename)
{
  ofstream out(filename, ios::app | ios::binary);
  
  out.flush();
  
  // out<<"//<<<<restart parameters: glbl_total_num_pkg | scale | dl | num_pkg | glbl_num_pkg | glbl_pkg_start | glbl_pkg_end | pkg_start | pkg_end | num_pkg_margin | dpkg | dcell \n";
       
  // out<<glbl_total_num_pkg<<" "<<scale<<" "<<dl<<" "
  //      <<num_pkg.i<<" "<<num_pkg.j<<" "<<num_pkg.k<<" "
  //      <<glbl_num_pkg.i<<" "<<glbl_num_pkg.j<<" "<<glbl_num_pkg.k<<" "
  //      <<glbl_pkg_start.i<<" "<<glbl_pkg_start.j<<" "<<glbl_pkg_start.k<<" "
  //      <<glbl_pkg_end.i<<" "<<glbl_pkg_end.j<<" "<<glbl_pkg_end.k<<" "
  //      <<pkg_start.i<<" "<<pkg_start.j<<" "<<pkg_start.k<<" "
  //      <<pkg_end.i<<" "<<pkg_end.j<<" "<<pkg_end.k<<" "
  //      <<num_pkg_margin.i<<" "<<num_pkg_margin.j<<" "<<num_pkg_margin.k<<" "
  //      <<dpkg.i<<" "<<dpkg.j<<" "<<dpkg.k<<" "
  //      <<dcell.i<<" "<<dcell.j<<" "<<dcell.k<<"\n";

  out.write(reinterpret_cast<char*>(&glbl_total_num_pkg), sizeof(int ));
  out.write(reinterpret_cast<char*>(&scale)             , sizeof(Real));
  out.write(reinterpret_cast<char*>(&dl)                , sizeof(Real));
  out.write(reinterpret_cast<char*>(&num_pkg.i)         , sizeof(int ));
  out.write(reinterpret_cast<char*>(&num_pkg.j)         , sizeof(int ));
  out.write(reinterpret_cast<char*>(&num_pkg.k)         , sizeof(int ));
  out.write(reinterpret_cast<char*>(&glbl_num_pkg.i)    , sizeof(int ));
  out.write(reinterpret_cast<char*>(&glbl_num_pkg.j)    , sizeof(int ));
  out.write(reinterpret_cast<char*>(&glbl_num_pkg.k)    , sizeof(int ));
  out.write(reinterpret_cast<char*>(&glbl_pkg_start.i)  , sizeof(int ));
  out.write(reinterpret_cast<char*>(&glbl_pkg_start.j)  , sizeof(int ));
  out.write(reinterpret_cast<char*>(&glbl_pkg_start.k)  , sizeof(int ));
  out.write(reinterpret_cast<char*>(&glbl_pkg_end.i)    , sizeof(int ));
  out.write(reinterpret_cast<char*>(&glbl_pkg_end.j)    , sizeof(int ));
  out.write(reinterpret_cast<char*>(&glbl_pkg_end.k)    , sizeof(int ));
  out.write(reinterpret_cast<char*>(&pkg_start.i)       , sizeof(int ));
  out.write(reinterpret_cast<char*>(&pkg_start.j)       , sizeof(int ));
  out.write(reinterpret_cast<char*>(&pkg_start.k)       , sizeof(int ));
  out.write(reinterpret_cast<char*>(&pkg_end.i)         , sizeof(int ));
  out.write(reinterpret_cast<char*>(&pkg_end.j)         , sizeof(int ));
  out.write(reinterpret_cast<char*>(&pkg_end.k)         , sizeof(int ));
  out.write(reinterpret_cast<char*>(&num_pkg_margin.i)  , sizeof(int ));
  out.write(reinterpret_cast<char*>(&num_pkg_margin.j)  , sizeof(int ));
  out.write(reinterpret_cast<char*>(&num_pkg_margin.k)  , sizeof(int ));
  out.write(reinterpret_cast<char*>(&dpkg.i)            , sizeof(Real));
  out.write(reinterpret_cast<char*>(&dpkg.j)            , sizeof(Real));
  out.write(reinterpret_cast<char*>(&dpkg.k)            , sizeof(Real));
  out.write(reinterpret_cast<char*>(&dcell.i)           , sizeof(Real));
  out.write(reinterpret_cast<char*>(&dcell.j)           , sizeof(Real));
  out.write(reinterpret_cast<char*>(&dcell.k)           , sizeof(Real));

  out.flush();
  
  // out<<"//<<<<restart parameters: phi | n_x | n_y | n_z | curv | psi | interface | characteristic | idx_characteristic | h\n";
  // phi
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        // out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].phi<<" ";
        Real var = table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].phi; 
        out.write(reinterpret_cast<char*>(&var), sizeof(var));
      }
      // out<<"\n";
    }
  }

  // n_x
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        // out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_x<<" ";
        Real var = table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_x; 
        out.write(reinterpret_cast<char*>(&var), sizeof(var));
      }
      // out<<"\n";
    }
  }

  // n_y
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        // out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_y<<" ";
        Real var = table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_y; 
        out.write(reinterpret_cast<char*>(&var), sizeof(var));
      }
      // out<<"\n";
    }
  }

  // n_z
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        // out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_z<<" ";
        Real var = table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_z; 
        out.write(reinterpret_cast<char*>(&var), sizeof(var));
      }
      // out<<"\n";
    }
  }

  // curv
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        // out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].curv<<" ";
        Real var = table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].curv; 
        out.write(reinterpret_cast<char*>(&var), sizeof(var));
      }
      // out<<"\n";
    }
  }

  // psi
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        // out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].psi<<" ";
        Real var = table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].psi; 
        out.write(reinterpret_cast<char*>(&var), sizeof(var));
      }
      // out<<"\n";
    }
  }

  // interface
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        // out<<int(table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].tag_interface)<<" ";
        int var = table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].tag_interface; 
        out.write(reinterpret_cast<char*>(&var), sizeof(var));
      }
      // out<<"\n";
    }
  }
  
  // interface
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        // out<<int(table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].tag_characteristic)<<" ";
        int var = table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].tag_characteristic; 
        out.write(reinterpret_cast<char*>(&var), sizeof(var));
      }
      // out<<"\n";
    }
  }

  // interface
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        // out<<int(table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].idx_characteristic)<<" ";
        int var = table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].idx_characteristic; 
        out.write(reinterpret_cast<char*>(&var), sizeof(var));
      }
      // out<<"\n";
    }
  }
  
  // scale;
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        // out<<table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].scale<<" ";
        Real var = table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].scale; 
        out.write(reinterpret_cast<char*>(&var), sizeof(var));
      }
      // out<<"\n";
    }
  }
  out.close();
}
//------------------------------------------------------------
// Write_level_set_rstfile
//------------------------------------------------------------
void Levelset_levelinfo::Load_level_set_rstfile(ifstream &load, Levelset *level_set, communicator &world)
{
  // char buffer[256];
  // load.getline (buffer,256);
  // load.getline (buffer,256);
  
  // load>>glbl_total_num_pkg>>scale>>dl
  //       >>num_pkg.i>>num_pkg.j>>num_pkg.k
  //       >>glbl_num_pkg.i>>glbl_num_pkg.j>>glbl_num_pkg.k
  //       >>glbl_pkg_start.i>>glbl_pkg_start.j>>glbl_pkg_start.k
  //       >>glbl_pkg_end.i>>glbl_pkg_end.j>>glbl_pkg_end.k
  //       >>pkg_start.i>>pkg_start.j>>pkg_start.k
  //       >>pkg_end.i>>pkg_end.j>>pkg_end.k
  //       >>num_pkg_margin.i>>num_pkg_margin.j>>num_pkg_margin.k
  //       >>dpkg.i>>dpkg.j>>dpkg.k
  //       >>dcell.i>>dcell.j>>dcell.k;
  
  // load.getline (buffer,256);
  // load.getline (buffer,256);
  
  load.read(reinterpret_cast<char*>(&glbl_total_num_pkg), sizeof(int ));
  load.read(reinterpret_cast<char*>(&scale)             , sizeof(Real));
  load.read(reinterpret_cast<char*>(&dl)                , sizeof(Real));
  load.read(reinterpret_cast<char*>(&num_pkg.i)         , sizeof(int ));
  load.read(reinterpret_cast<char*>(&num_pkg.j)         , sizeof(int ));
  load.read(reinterpret_cast<char*>(&num_pkg.k)         , sizeof(int ));
  load.read(reinterpret_cast<char*>(&glbl_num_pkg.i)    , sizeof(int ));
  load.read(reinterpret_cast<char*>(&glbl_num_pkg.j)    , sizeof(int ));
  load.read(reinterpret_cast<char*>(&glbl_num_pkg.k)    , sizeof(int ));
  load.read(reinterpret_cast<char*>(&glbl_pkg_start.i)  , sizeof(int ));
  load.read(reinterpret_cast<char*>(&glbl_pkg_start.j)  , sizeof(int ));
  load.read(reinterpret_cast<char*>(&glbl_pkg_start.k)  , sizeof(int ));
  load.read(reinterpret_cast<char*>(&glbl_pkg_end.i)    , sizeof(int ));
  load.read(reinterpret_cast<char*>(&glbl_pkg_end.j)    , sizeof(int ));
  load.read(reinterpret_cast<char*>(&glbl_pkg_end.k)    , sizeof(int ));
  load.read(reinterpret_cast<char*>(&pkg_start.i)       , sizeof(int ));
  load.read(reinterpret_cast<char*>(&pkg_start.j)       , sizeof(int ));
  load.read(reinterpret_cast<char*>(&pkg_start.k)       , sizeof(int ));
  load.read(reinterpret_cast<char*>(&pkg_end.i)         , sizeof(int ));
  load.read(reinterpret_cast<char*>(&pkg_end.j)         , sizeof(int ));
  load.read(reinterpret_cast<char*>(&pkg_end.k)         , sizeof(int ));
  load.read(reinterpret_cast<char*>(&num_pkg_margin.i)  , sizeof(int ));
  load.read(reinterpret_cast<char*>(&num_pkg_margin.j)  , sizeof(int ));
  load.read(reinterpret_cast<char*>(&num_pkg_margin.k)  , sizeof(int ));
  load.read(reinterpret_cast<char*>(&dpkg.i)            , sizeof(Real));
  load.read(reinterpret_cast<char*>(&dpkg.j)            , sizeof(Real));
  load.read(reinterpret_cast<char*>(&dpkg.k)            , sizeof(Real));
  load.read(reinterpret_cast<char*>(&dcell.i)           , sizeof(Real));
  load.read(reinterpret_cast<char*>(&dcell.j)           , sizeof(Real));
  load.read(reinterpret_cast<char*>(&dcell.k)           , sizeof(Real));

    // phi
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        Real var = 0.;
        load.read(reinterpret_cast<char*>(&var), sizeof(Real));
        table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].phi = var;
        // load>>table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].phi;
      }
    }
  }

  // n_x
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        Real var = 0.;
        load.read(reinterpret_cast<char*>(&var), sizeof(Real));
        table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_x = var;
        // load>>table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_x;
      }
    }
  }

  // n_y
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        Real var = 0.;
        load.read(reinterpret_cast<char*>(&var), sizeof(Real));
        table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_y = var;
        // load>>table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_y;
      }
    }
  }

  // n_z
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        Real var = 0.;
        load.read(reinterpret_cast<char*>(&var), sizeof(Real));
        table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_z = var;
        // load>>table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].n_z;
      }
    }
  }

  // curv
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        Real var = 0.;
        load.read(reinterpret_cast<char*>(&var), sizeof(Real));
        table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].curv = var;
        // load>>table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].curv;
      }
    }
  }

  // psi
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        Real var = 0.;
        load.read(reinterpret_cast<char*>(&var), sizeof(Real));
        table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].psi = var;
        // load>>table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].psi;
      }
    }
  }

  // interface
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        int var = 0;
        load.read(reinterpret_cast<char*>(&var), sizeof(int));
        table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].tag_interface = var;
       // load>>table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].tag_interface;
      }
    }
  }
  
  // interface
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        int var = 0;
        load.read(reinterpret_cast<char*>(&var), sizeof(int));
        table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].tag_characteristic = var;
        // load>>table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].tag_characteristic;
      }
    }
  }

  // interface
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        int var = 0;
        load.read(reinterpret_cast<char*>(&var), sizeof(int));
        table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].idx_characteristic = var;
        // load>>table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].idx_characteristic;
      }
    }
  }
  
  // scale;
  for(int k=AMAX1(pkg_start.k, 0)*ICPZ; k<AMIN1(glbl_num_pkg.k, num_pkg.k)*ICPZ; k++){
    for(int j=AMAX1(pkg_start.j, 0)*ICPY; j<AMIN1(glbl_num_pkg.j, num_pkg.j)*ICPY; j++){
      for(int i=AMAX1(pkg_start.i, 0)*ICPX; i<AMIN1(glbl_num_pkg.i, num_pkg.i)*ICPX; i++){
        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;
        Real var = 0.;
        load.read(reinterpret_cast<char*>(&var), sizeof(Real));
        table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].scale = var;
        // load>>table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].scale;
      }
    }
  }
}


#ifdef _READ_SDF_
//------------------------------------------------------------
// Load_SDF_file
//------------------------------------------------------------
void Levelset_levelinfo::Load_SDF_file(Levelset *level_set, communicator &world)
{
  char sdf_file[256];
  size_t _size = sizeof (level_set->sdf_file);
  snprintf(sdf_file, _size, "%s", level_set->sdf_file);
  my_int cells_to_cut = level_set->cells_to_cut;

  if (DIM == 2){
    cout<<"<<<<< ["<<world.rank()<<"]->ERROR!!! Unable to read 2D SDF files"<<endl;
    world.abort(-1);
  }

  ifstream load(sdf_file);
  
  if (! load.is_open())  { 
    cout << "<<<<< ["<<world.rank()<<"]->Error opening .sdf file."<<endl; 
    world.abort (-1); 
  }
  cout<<"<<<<< ["<<world.rank()<<"]-> loading SDF file....."<<endl;

  my_int grid_size;
  Real dx = 0.;
  Real tmp = 0.;

  load>>grid_size.i>>grid_size.j>>grid_size.k;
  load>>tmp        >>tmp        >>tmp;
  load>>dx;

  if (fabs(dx-dl)> 1.e-15 ){
    cout<<"<<<<< ["<<world.rank()<<"]->ERROR!!! dx ("<<dx<<") is different from dcell ("<<dl<<")."<<endl;
    world.abort(-1);
  }

  // array_real phis;
  // phis.resize (boost::extents[mrange(0,grid_size.i)][mrange(0,grid_size.j)][mrange(0,grid_size.k)]);

  // for (int k = 0; k < grid_size.k; k++){
  //   for (int j = 0; j < grid_size.j; j++)
  //     for (int i = 0; i < grid_size.i; i++)
  //       load>>phis[i][j][k];
  //     }

  // for(int k=0; k<num_pkg.k*ICPZ; k++){
  //   for(int j=0; j<num_pkg.j*ICPY; j++){
  //     for(int i=0; i<num_pkg.i*ICPX; i++){
  for (int k = 0; k < grid_size.k; k++){
    for (int j = 0; j < grid_size.j; j++){
      for (int i = 0; i < grid_size.i; i++){

        int pkg_i = int(i/ICPX);
        int pkg_j = int(j/ICPY);
        int pkg_k = int(k/ICPZ);
        int cell_i = i - pkg_i*ICPX;
        int cell_j = j - pkg_j*ICPY;
        int cell_k = k - pkg_k*ICPZ;

        Real phi_ = 0.;
        load>>phi_;

        if (i < num_pkg.i*ICPX && j < num_pkg.j*ICPY && k < num_pkg.k*ICPZ)
          table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].Set_phi(phi_);
      }
    }
  }
  cout<<"<<<<< ["<<world.rank()<<"]-> loading SDF file finished."<<endl;

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Load_SDF_file finished\n";
  out.close();
#endif
}
#endif