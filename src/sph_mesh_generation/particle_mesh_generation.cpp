#include "glbfunc.h"
#include "voronoi_construction.h"
#include "graph_mesh.h"
#include "sph_mesh_generation.h"
#include "particle_mesh_generation.h"
#include "level_set.h"
#include "lset_level_infor.h"
#include "lset_package.h"
#include "lset_cell.h"

/***************************************************/
/*                                                 */
/*      Functions defined in class "Particle"      */
/*                                                 */
/***************************************************/

//-----------------------------------------------------
// constructor
//-----------------------------------------------------
Particle_mesh_generation::Particle_mesh_generation()
{
  type       = REAL_PARTICLE;
  level      = 0;
  color      = 0;
  P          = 0.;
  p_mass     = 0.;
  sum        = 0.;
  rho        = 0.;
  mass       = 0.;
  vol        = 0.;
  h          = 0.;
  scale      = 0.;
  id         = 0;
  local_id   = 0;
  safe_guard = 0.;
  h_nominal  = 0.;
  phi        = 0.;
  curv       = 0.;
  nu         = 0.;

  my_set_const (coord, 0.);
  my_set_const (v, 0.0);
  my_set_const (a, 0.0);
  my_set_const (norm, 0.0);
  my_set_const (coord_corr, 0.);
#if PERI_DIM != 0 || SYM_DIM != 0
  copy.clear();
#endif
}
//-----------------------------------------------------------------
// Initialize
//-----------------------------------------------------------------
void Particle_mesh_generation::Initialize (SOLVER *sph, int id_, int type_, int rank_)
{
  id         = id_;
  color      = rank_;
  scale      = sph->max_scale;
  safe_guard = sph->safe_guard;
  h          = scale / CUT_OFF;
  mass       = 1.;
  P          = 1.;
  type       = type_;
}
//-----------------------------------------------------------------
// Set_local_index
//-----------------------------------------------------------------
void Particle_mesh_generation::Set_local_index (int index_)
{
  local_id   = index_;
}
//-----------------------------------------------------------------
// clear particle info
//-----------------------------------------------------------------
void Particle_mesh_generation::Clearup(){
  neighbor.clear();

#if PERI_DIM != 0 || SYM_DIM != 0
  copy.clear();
#endif
}
//-----------------------------------------------------
// reset particle neighbor information
//-----------------------------------------------------
void Particle_mesh_generation::Reset_neighbor_info()
{
  neighbor.clear();

}
//-----------------------------------------------------
// add neighbor for particle
//-----------------------------------------------------
void Particle_mesh_generation::Add_neighbor(Particle *current_particle)
{
  neighbor.push_back(current_particle);
}
//-----------------------------------------------------
// Refresh neighbor infor
//-----------------------------------------------------
void Particle_mesh_generation::Refresh_neighbor_info(SPH *sph, int flag)
{
  p_Level_info    current_level = sph->level_info[level-sph->Lmin];
  my_real coord_shift = my_minus_data (coord, sph->box_l);
#ifndef _SCLL_
  my_int  pos = get_cell_id (coord_shift, current_level->dcell, current_level->cell_start, current_level->cell_end);
#else
  my_int  pos        = get_cell_id (coord_shift, current_level->dsubcell, current_level->subcell_start, current_level->subcell_end);
#endif
  current_level->Refresh_neighbor_info(sph, this, pos, flag);
}
//-----------------------------------------------------------------
// set time step
//-----------------------------------------------------------------
void Particle_mesh_generation::Set_timestep()
{
  Real    U = get_distance(v);
  
  timestep  = 0.25*sqrt(h/(get_distance(a)+1.e-10));
  
  timestep  = AMIN1(timestep, h/40./(U+1.e-10));
  
  nu        = 0.5*h*U;
  timestep  = AMIN1(timestep, 0.125*h*h/(nu+1.e-10));
}
//-----------------------------------------------------------------
// set time step
//-----------------------------------------------------------------
void Particle_mesh_generation::Set_timestep(int stage, SOLVER *sph)
{
  Real   U   = get_distance(v);
  
  timestep   = 0.25*sqrt(h/(get_distance(a)+1.e-15));
  
  timestep   = AMIN1(timestep, h/40./(U+1.e-15));
  
  // if (stage  == 2 && sph->current_stage != FILL){
  if (stage == 2 ){
    Real nu_coeff = my_get_ramping_factor (sph->nu_coeff_start, 
                                           sph->nu_coeff_end, 
                                           sph->run_time, 
                                           sph->nu_ramping_start, 
                                           sph->nu_ramping_end);
    nu       = nu_coeff*h*U;
    timestep = AMIN1(timestep, 0.125*h*h/(nu+1.e-15));
  }
}
//-----------------------------------------------------------------
// Reset particle force
//-----------------------------------------------------------------
void Particle_mesh_generation::Reset_particle_force(SOLVER *sph, int stage)
{
  if (stage == 1) sum = 0.;

  my_set_const (a, 0.0);
  my_set_const (F1, 0.0);
  my_set_const (F2, 0.0);
  my_set_const (coord_corr, 0.);
  
  // if (stage == 1){
  //   if (   (sph->current_stage == FILL && int(sph->run_time)%10 == 0)
  //        || (sph->current_stage != FILL )){
  //     my_set_const (v, 0.0);
  //   }
  // }
  // if (stage == 1 && int(sph->run_time)%40 == 0){
  if (stage == 1){
    if (sph->run_time < sph->v_reini_change && int(sph->run_time)%sph->v_reini_start == 0)
      my_set_const (v, 0.0);
    else if (sph->run_time >= sph->v_reini_change && int(sph->run_time)%sph->v_reini_end == 0)
      my_set_const (v, 0.0);
  }
}
//-----------------------------------------------------------------
// Calculate_particle_infor
//-----------------------------------------------------------------
void Particle_mesh_generation::Calculate_particle_infor(SOLVER *sph)
{
  p_Levelset_levelinfo current_level = sph->level_set.lset_level_info[0];
  my_real lset_domain = sph->level_set.domain;
  
  // my_real coord_shift = my_minus_data (coord, sph->box_l);

  my_real coord_shift;
  coord_shift.i = DIM_X == 1 ? AMIN1(AMAX1(coord.i-sph->box_l.i, 0.), lset_domain.i-1.e-8) : 0.;
  coord_shift.j = DIM_Y == 1 ? AMIN1(AMAX1(coord.j-sph->box_l.j, 0.), lset_domain.j-1.e-8) : 0.;
  coord_shift.k = DIM_Z == 1 ? AMIN1(AMAX1(coord.k-sph->box_l.k, 0.), lset_domain.k-1.e-8) : 0.;
  
  my_int  pos = get_cell_id (coord_shift, current_level->dpkg, current_level->glbl_pkg_start, current_level->glbl_pkg_end);
  
  p_Levelset_package current_pkg = current_level->table_lset_pkg_list[pos.i][pos.j][pos.k];

  my_real rr; rr.i = rr.j = rr.k = 0.;
  my_real dcell = current_level->dcell;
  my_real dpkg = current_level->dpkg;
  
  my_real coord_shift_pkg = my_minus_data (coord_shift, my_multiply_data (dpkg, pos));
  
  rr.i = DIM_X==1 ? coord_shift_pkg.i/dcell.i - 0.5 : 0.;
  rr.j = DIM_Y==1 ? coord_shift_pkg.j/dcell.j - 0.5 : 0.;
  rr.k = DIM_Z==1 ? coord_shift_pkg.k/dcell.k - 0.5 : 0.;
  
  rr.i = AMAX1(-0.5,AMIN1(rr.i,Real(ICPX)-0.5));
  rr.j = AMAX1(-0.5,AMIN1(rr.j,Real(ICPY)-0.5));
  rr.k = AMAX1(-0.5,AMIN1(rr.k,Real(ICPZ)-0.5));

  int     im = 0;
  int     jm = 0;
  int     km = 0;

  Real ii = 0.;
  Real jj = 0.;
  Real kk = 0.;

  #if DIM_X
    im = floor (rr.i);
    ii =  rr.i - Real(im);
  #endif
  #if DIM_Y
    jm = floor (rr.j);
    jj = rr.j - Real(jm);
  #endif
  #if DIM_Z
    km = floor (rr.k);
    kk =  rr.k - Real(km);
  #endif

  Real Target_smooth_infor_scale[1+DIM_X][1+DIM_Y][1+DIM_Z];
  Real Target_nx                [1+DIM_X][1+DIM_Y][1+DIM_Z];
  Real Target_ny                [1+DIM_X][1+DIM_Y][1+DIM_Z];
  Real Target_nz                [1+DIM_X][1+DIM_Y][1+DIM_Z];
  Real Target_curv              [1+DIM_X][1+DIM_Y][1+DIM_Z];
  Real coeffecient              [1+DIM_X][1+DIM_Y][1+DIM_Z];
  
  for(int r=0; r<=DIM_X; r++)
    for(int s=0; s<=DIM_Y; s++)
      for(int t=0; t<=DIM_Z; t++){
        Target_smooth_infor_scale[r][s][t] = current_pkg->p_cell[r+im][s+jm][t+km]->phi;
        Target_nx                [r][s][t] = current_pkg->p_cell[r+im][s+jm][t+km]->n_x;
        Target_ny                [r][s][t] = current_pkg->p_cell[r+im][s+jm][t+km]->n_y;
        Target_nz                [r][s][t] = current_pkg->p_cell[r+im][s+jm][t+km]->n_z;
        Target_curv              [r][s][t] = current_pkg->p_cell[r+im][s+jm][t+km]->curv;

        coeffecient[r][s][t] = (DIM_X==1 ? fabs(ii-(1.0-r)):1.0)*(DIM_Y==1 ? fabs(jj-(1.0-s)):1.0)*(DIM_Z==1 ? fabs(kk-(1.0-t)):1.0);
      }

  phi = 0. ;
  norm.i = 0.;
  norm.j = 0.;
  norm.k = 0.;
  curv = 0. ;

  for(int r=0; r<=DIM_X; r++)
    for(int s=0; s<=DIM_Y; s++)
      for(int t=0; t<=DIM_Z; t++){
        phi    += Target_smooth_infor_scale[r][s][t]*coeffecient[r][s][t];
        norm.i += Target_nx[r][s][t]*coeffecient[r][s][t];
        norm.j += Target_ny[r][s][t]*coeffecient[r][s][t];
        norm.k += Target_nz[r][s][t]*coeffecient[r][s][t];
        curv   += Target_curv[r][s][t]*coeffecient[r][s][t];
      }

  // renormalize the normal direction
  Real normal = get_distance (norm) + 1.e-20;
  my_multiply_const (norm, 1./normal);
  
}
//-----------------------------------------------------------------
// Calculate_particle_infor
//-----------------------------------------------------------------
void Particle_mesh_generation::Calculate_particle_infor(Levelset *level_set)
{
  p_Levelset_levelinfo current_level = level_set->lset_level_info[0];
  my_real lset_domain = level_set->domain;
  
  my_real coord_shift;
  coord_shift.i = DIM_X == 1 ? AMIN1(AMAX1(coord.i-level_set->box_l.i, 0.), lset_domain.i-1.e-8) : 0.;
  coord_shift.j = DIM_Y == 1 ? AMIN1(AMAX1(coord.j-level_set->box_l.j, 0.), lset_domain.j-1.e-8) : 0.;
  coord_shift.k = DIM_Z == 1 ? AMIN1(AMAX1(coord.k-level_set->box_l.k, 0.), lset_domain.k-1.e-8) : 0.;

  my_int  pos = get_cell_id (coord_shift, current_level->dpkg, current_level->glbl_pkg_start, current_level->glbl_pkg_end);
  
  p_Levelset_package current_pkg = current_level->table_lset_pkg_list[pos.i][pos.j][pos.k];

  my_real rr; rr.i = rr.j = rr.k = 0.;
  my_real dcell = current_level->dcell;
  my_real dpkg = current_level->dpkg;
  
  my_real coord_shift_pkg = my_minus_data (coord_shift, my_multiply_data (dpkg, pos));
  
  rr.i = DIM_X==1 ? coord_shift_pkg.i/dcell.i - 0.5 : 0.;
  rr.j = DIM_Y==1 ? coord_shift_pkg.j/dcell.j - 0.5 : 0.;
  rr.k = DIM_Z==1 ? coord_shift_pkg.k/dcell.k - 0.5 : 0.;
  
  rr.i = AMAX1(-0.5,AMIN1(rr.i,Real(ICPX)-0.5));
  rr.j = AMAX1(-0.5,AMIN1(rr.j,Real(ICPY)-0.5));
  rr.k = AMAX1(-0.5,AMIN1(rr.k,Real(ICPZ)-0.5));

  int     im = 0;
  int     jm = 0;
  int     km = 0;

  Real ii = 0.;
  Real jj = 0.;
  Real kk = 0.;

  #if DIM_X
    im = floor (rr.i);
    ii =  rr.i - Real(im);
  #endif
  #if DIM_Y
    jm = floor (rr.j);
    jj = rr.j - Real(jm);
  #endif
  #if DIM_Z
    km = floor (rr.k);
    kk =  rr.k - Real(km);
  #endif

  Real Target_smooth_infor_scale[1+DIM_X][1+DIM_Y][1+DIM_Z];
  Real Target_nx                [1+DIM_X][1+DIM_Y][1+DIM_Z];
  Real Target_ny                [1+DIM_X][1+DIM_Y][1+DIM_Z];
  Real Target_nz                [1+DIM_X][1+DIM_Y][1+DIM_Z];
  Real Target_curv              [1+DIM_X][1+DIM_Y][1+DIM_Z];
  Real coeffecient              [1+DIM_X][1+DIM_Y][1+DIM_Z];
  
  for(int r=0; r<=DIM_X; r++)
    for(int s=0; s<=DIM_Y; s++)
      for(int t=0; t<=DIM_Z; t++){
        Target_smooth_infor_scale[r][s][t] = current_pkg->p_cell[r+im][s+jm][t+km]->phi;
        Target_nx                [r][s][t] = current_pkg->p_cell[r+im][s+jm][t+km]->n_x;
        Target_ny                [r][s][t] = current_pkg->p_cell[r+im][s+jm][t+km]->n_y;
        Target_nz                [r][s][t] = current_pkg->p_cell[r+im][s+jm][t+km]->n_z;
        Target_curv              [r][s][t] = current_pkg->p_cell[r+im][s+jm][t+km]->curv;

        coeffecient[r][s][t] = (DIM_X==1 ? fabs(ii-(1.0-r)):1.0)*(DIM_Y==1 ? fabs(jj-(1.0-s)):1.0)*(DIM_Z==1 ? fabs(kk-(1.0-t)):1.0);
      }

  phi = 0. ;
  norm.i = 0.;
  norm.j = 0.;
  norm.k = 0.;
  curv = 0. ;

  for(int r=0; r<=DIM_X; r++)
    for(int s=0; s<=DIM_Y; s++)
      for(int t=0; t<=DIM_Z; t++){
        phi    += Target_smooth_infor_scale[r][s][t]*coeffecient[r][s][t];
        norm.i += Target_nx[r][s][t]*coeffecient[r][s][t];
        norm.j += Target_ny[r][s][t]*coeffecient[r][s][t];
        norm.k += Target_nz[r][s][t]*coeffecient[r][s][t];
        curv   += Target_curv[r][s][t]*coeffecient[r][s][t];
      }

  // renormalize the normal direction
  Real normal = get_distance (norm) + 1.e-20;
  my_multiply_const (norm, 1./normal);
  
}
//-----------------------------------------------------------------
// Set_particle_scale
//-----------------------------------------------------------------
void Particle_mesh_generation::Set_particle_info(SOLVER *sph)
{
  Real h_tmp = sph->level_set.Get_scale(curv, phi);
  // h = h_tmp;
  // vol = pow(h, DIM);

  h = h_tmp*CELL_RATIO/CUT_OFF;
  vol = pow(h_tmp, DIM);

  // JZ20181214 :: play around with using different h vs. dx here. 
  //               Trying to improve convergence

  // if      (type == REAL_PARTICLE)        h = h_tmp*CELL_RATIO/CUT_OFF;
  // else if (type == SURFACE_PARTICLE)     h = h_tmp*CELL_RATIO/CUT_OFF*1.2;
  // else if (type == SEGMENT_PARTICLE)     h = h_tmp*CELL_RATIO/CUT_OFF*1.2;
  // else if (type == SINGULARITY_PARTICLE) h = h_tmp*CELL_RATIO/CUT_OFF*1.2;

  // vol = pow(h_tmp, DIM);

  rho = mass/vol;
}
//-----------------------------------------------------------------
// Get_dispersed_position
//-----------------------------------------------------------------
void Particle_mesh_generation::Get_dispersed_position(Real dx)
{
  Real     n_d = 0.;
  Real     div = 1.0/(get_distance(norm) + 1.0e-15);
  my_real norm_ = my_multiply_const (norm, div);
  
#if (DIM_X)
  n_d     = phi < 0.0 ? norm_.i : -norm_.i;
  coord.i = coord.i + n_d * dx;
#endif
#if (DIM_Y)
  n_d     = phi < 0.0 ? norm_.j : -norm_.j;
  coord.j = coord.j + n_d * dx;
#endif
#if (DIM_Z)
  n_d     = phi < 0.0 ? norm_.k : -norm_.k;
  coord.k = coord.k + n_d * dx;
#endif

}
//-----------------------------------------------------------------
// Find_and_map_particle
//-----------------------------------------------------------------
void Particle_mesh_generation::Find_and_map_particle(SOLVER *sph)
{
  int tag_interface = NORMAL_CELL;
  int tag_characteristic = NORMAL_CELL;
  int idx_characteristic = -1;

  // Get_extended_char_cell_tag_and_idx (sph, tag_interface, tag_characteristic, idx_characteristic);
  Get_level_set_char_cell_tag_and_idx (sph, tag_interface, tag_characteristic, idx_characteristic);

  if ( (tag_interface == CUT_CELL  && type == REAL_PARTICLE && tag_characteristic == NORMAL_CELL) ||
       (tag_interface == CUT_CELL  && type == SURFACE_PARTICLE && tag_characteristic == NORMAL_CELL)){
      type = SURFACE_PARTICLE;

      Real n_d  = 0.;
      Real dist = fabs(phi);
      Real div  = 1.0/(get_distance(norm) + 1.0e-15);
      
      my_real norm_      = my_multiply_const (norm, div);
      
    #if (DIM_X)                                             
      n_d     = phi < 0.0 ? norm_.i : -norm_.i;
      coord.i = coord.i + n_d * dist;
    #endif
    #if (DIM_Y)
      n_d     = phi < 0.0 ? norm_.j : -norm_.j;
      coord.j = coord.j + n_d * dist;
    #endif
    #if (DIM_Z)
      n_d     = phi < 0.0 ? norm_.k : -norm_.k;
      coord.k = coord.k + n_d * dist;
    #endif
  }else if (tag_characteristic == SEGMENT_CELL ){
    type = SEGMENT_PARTICLE;

    if (idx_characteristic < 0) 
      cout<<"<<<<< ERROR!!!! wrong returned idx_characteristic"<<endl;
    else{
      std::pair <my_real, my_real> seg_i = sph->level_set.segment[idx_characteristic];

      my_real Pij = my_minus_data     (seg_i.second, seg_i.first);
      my_real Nij = my_minus_data     (       coord, seg_i.first);
      my_real eij = my_multiply_const (         Pij, 1./AMAX1(get_distance (Pij), 1.e-15));
      Real   dist = get_dot_product   (         Nij, eij);
      my_real nij = my_multiply_const (         eij, dist);
      my_real  cn = my_add_data       ( seg_i.first, nij);

      my_set_data (coord, cn);
    }
  }else if (tag_characteristic == SINGULARITY_CELL){
    type = SINGULARITY_PARTICLE;

    if (idx_characteristic < 0) 
      cout<<"<<<<< ERROR!!!! wrong returned idx_characteristic"<<endl;
    else{
      my_set_data (coord, sph->level_set.singularity[idx_characteristic]);
    }
  } else {
    type = REAL_PARTICLE;
  }
}
//-----------------------------------------------------------------
// Find_and_map_surface_particle
//-----------------------------------------------------------------
void Particle_mesh_generation::Find_and_map_surface_particle(SOLVER *sph)
{
  int tag_interface = NORMAL_CELL;
  int tag_characteristic = NORMAL_CELL;
  
  Get_current_particle_level_set_cell_tag (sph, tag_interface, tag_characteristic);
  
  if ( tag_interface == CUT_CELL  && type == REAL_PARTICLE){
      type = SURFACE_PARTICLE;
  }
  
  if (type == SURFACE_PARTICLE){
    Real n_d  = 0.;
    Real dist = fabs(phi);
    Real div  = 1.0/(get_distance(norm) + 1.0e-15);
    
    my_real norm_      = my_multiply_const (norm, div);
    
  #if (DIM_X)                                             
    n_d     = phi < 0.0 ? norm_.i : -norm_.i;
    coord.i = coord.i + n_d * dist;
  #endif
  #if (DIM_Y)
    n_d     = phi < 0.0 ? norm_.j : -norm_.j;
    coord.j = coord.j + n_d * dist;
  #endif
  #if (DIM_Z)
    n_d     = phi < 0.0 ? norm_.k : -norm_.k;
    coord.k = coord.k + n_d * dist;
  #endif
  }
}
//-----------------------------------------------------------------
// Find_and_map_segment_particle
//-----------------------------------------------------------------
void Particle_mesh_generation::Find_and_map_segment_particle(SOLVER *sph)
{
  int tag_interface = NORMAL_CELL;
  int tag_characteristic = NORMAL_CELL;
  int idx_characteristic = -1;

  // Get_level_set_char_cell_tag_and_idx (sph, tag_interface, tag_characteristic, idx_characteristic);

  Get_extended_char_cell_tag_and_idx (sph, tag_interface, tag_characteristic, idx_characteristic);

  if (tag_characteristic == SEGMENT_CELL){
    type = SEGMENT_PARTICLE;

    if (idx_characteristic < 0) 
      cout<<"<<<<< ERROR!!!! wrong returned idx_characteristic"<<endl;
    else{
      std::pair <my_real, my_real> seg_i = sph->level_set.segment[idx_characteristic];

      my_real Pij = my_minus_data     (seg_i.second, seg_i.first);
      my_real Nij = my_minus_data     (       coord, seg_i.first);
      my_real eij = my_multiply_const (         Pij, 1./AMAX1(get_distance (Pij), 1.e-15));
      Real   dist = get_dot_product   (         Nij, eij);
      my_real nij = my_multiply_const (         eij, dist);
      my_real  cn = my_add_data       ( seg_i.first, nij);

      my_set_data (coord, cn);
    }
  }
}
//-----------------------------------------------------------------
// Get_extended_char_cell_tag_and_idx
//-----------------------------------------------------------------
void Particle_mesh_generation::Get_extended_char_cell_tag_and_idx(SOLVER *sph, int &tag_interface, int &tag_characteristic, int &idx_characteristic)
{
  p_Levelset_levelinfo current_level = sph->level_set.lset_level_info[0];
  
  my_real coord_shift = my_minus_data (coord, sph->box_l);

  // current_level->Search_for_char_cell_within_dx (coord_shift, current_level->dl*2.0, tag_interface, tag_characteristic, idx_characteristic);

  Real h_search = AMAX1 (current_level->dl*2.0, 0.5*h*CUT_OFF/CELL_RATIO);
  current_level->Search_for_char_cell_within_dx (coord_shift, h_search, tag_interface, tag_characteristic, idx_characteristic);
}
//-----------------------------------------------------------------
// Get_level_set_char_cell_tag_and_idx
//-----------------------------------------------------------------
void Particle_mesh_generation::Get_level_set_char_cell_tag_and_idx(SOLVER *sph, int &tag_interface, int &tag_characteristic, int &idx_characteristic)
{
  p_Levelset_levelinfo current_level = sph->level_set.lset_level_info[0];
  
  my_int  pos_pkg;
  my_int  pos_cell;

  my_real lset_domain = sph->level_set.domain;
  
  my_real coord_shift;
  coord_shift.i = DIM_X == 1 ? AMIN1(AMAX1(coord.i-sph->box_l.i, 0.), lset_domain.i-1.e-8) : 0.;
  coord_shift.j = DIM_Y == 1 ? AMIN1(AMAX1(coord.j-sph->box_l.j, 0.), lset_domain.j-1.e-8) : 0.;
  coord_shift.k = DIM_Z == 1 ? AMIN1(AMAX1(coord.k-sph->box_l.k, 0.), lset_domain.k-1.e-8) : 0.;

  current_level->Get_id_pkg_cell (coord_shift, pos_pkg, pos_cell);
  
  p_Levelset_package current_pkg = current_level->table_lset_pkg_list[pos_pkg.i][pos_pkg.j][pos_pkg.k];

  tag_interface      = current_pkg->p_cell[pos_cell.i][pos_cell.j][pos_cell.k]->tag_interface;
  tag_characteristic = current_pkg->p_cell[pos_cell.i][pos_cell.j][pos_cell.k]->tag_characteristic;
  idx_characteristic = current_pkg->p_cell[pos_cell.i][pos_cell.j][pos_cell.k]->idx_characteristic;
}
//-----------------------------------------------------------------
// Get_current_particle_level_set_cell_tag
//-----------------------------------------------------------------
void Particle_mesh_generation::Get_current_particle_level_set_cell_tag(SOLVER *sph, int &tag_interface, int &tag_characteristic)
{
  p_Levelset_levelinfo current_level = sph->level_set.lset_level_info[0];
  
  my_int  pos_pkg;
  my_int  pos_cell;
  
  my_real lset_domain = sph->level_set.domain;
  
  my_real coord_shift;
  coord_shift.i = DIM_X == 1 ? AMIN1(AMAX1(coord.i-sph->box_l.i, 0.), lset_domain.i-1.e-8) : 0.;
  coord_shift.j = DIM_Y == 1 ? AMIN1(AMAX1(coord.j-sph->box_l.j, 0.), lset_domain.j-1.e-8) : 0.;
  coord_shift.k = DIM_Z == 1 ? AMIN1(AMAX1(coord.k-sph->box_l.k, 0.), lset_domain.k-1.e-8) : 0.;
  current_level->Get_id_pkg_cell (coord_shift, pos_pkg, pos_cell);
  
  p_Levelset_package current_pkg = current_level->table_lset_pkg_list[pos_pkg.i][pos_pkg.j][pos_pkg.k];

  tag_interface = current_pkg->p_cell[pos_cell.i][pos_cell.j][pos_cell.k]->tag_interface;
  tag_characteristic = current_pkg->p_cell[pos_cell.i][pos_cell.j][pos_cell.k]->tag_characteristic;
}
//-----------------------------------------------------------------
// Find_and_map_characteristic_particle
//-----------------------------------------------------------------
void Particle_mesh_generation::Find_characteristic_particle(SOLVER *sph)
{
  int tag_interface = NORMAL_CELL;
  int tag_characteristic = NORMAL_CELL;
  
  Get_current_particle_level_set_cell_tag (sph, tag_interface, tag_characteristic);
  
  if (      tag_characteristic == SINGULARITY_CELL && type != SINGULARITY_PARTICLE){
    type = SINGULARITY_PARTICLE;
  }
}
//-----------------------------------------------------------------
// Find_and_map_characteristic_particle
//-----------------------------------------------------------------
void Particle_mesh_generation::Map_characteristic_particle(SOLVER *sph)
{
  if (type == SINGULARITY_PARTICLE){
    Real n_d = 0.;
    Real dist = fabs(phi);
    Real div = 1.0/(get_distance(norm) + 1.0e-15);
    
    my_real norm_ = my_multiply_const (norm, div);
    
  #if (DIM_X)
    n_d = phi < 0.0 ? norm_.i : -norm_.i;
    coord.i = coord.i + n_d * dist;
  #endif
  #if (DIM_Y)
    n_d = phi < 0.0 ? norm_.j : -norm_.j;
    coord.j = coord.j + n_d * dist;
  #endif
  #if (DIM_Z)
    n_d = phi < 0.0 ? norm_.k : -norm_.k;
    coord.k = coord.k + n_d * dist;
  #endif
  }
}
//-----------------------------------------------------------------
// Find_and_map_surface_particle
//-----------------------------------------------------------------
void Particle_mesh_generation::Pop_surface_particle(SOLVER *sph)
{
//   Real n_d = 0.;
//   Real dist = 1.1*h;
//   Real div = 1.0/(get_distance(norm) + 1.0e-15);
//     
//   norm = my_multiply_const (norm, div);
//     
//   #if (DIM_X)                                             
//   n_d = phi < 0.0 ? norm.i : -norm.i;
//   coord.i = coord.i + norm.i * dist;
//   #endif
//   #if (DIM_Y)
//   n_d = phi < 0.0 ? norm.j : -norm.j;
//   coord.j = coord.j + norm.j * dist;
//   #endif
//   #if (DIM_Z)
//   n_d = phi < 0.0 ? norm.k : -norm.k;
//   coord.k = coord.k + norm.k * dist;
//   #endif

  srand((unsigned)time(NULL)); 
  my_real position;
  
      bool positive_phase = false;
      while (!positive_phase ){
    #if (DIM_X)
        position.i = rand()/ double(RAND_MAX+1.0);
        coord.i = position.i*sph->particle_box.i + rand()/ double(RAND_MAX+1.0)*1.e-3 + sph->particle_box_l.i;
        coord.i = AMAX1(sph->particle_box_l.i+1.e-20,AMIN1(coord.i, sph->particle_box_r.i-1.e-20));
    #endif
    #if (DIM_Y)
        position.j = rand()/ double(RAND_MAX+1.0);
        coord.j = position.j*sph->particle_box.j + rand()/ double(RAND_MAX+1.0)*1.e-3 + sph->particle_box_l.j;
        coord.j = AMAX1(sph->particle_box_l.j+1.e-20,AMIN1(coord.j, sph->particle_box_r.j-1.e-20)); 
    #endif
    #if (DIM_Z)
        position.k = rand()/ double(RAND_MAX+1.0);
        coord.k = position.k*sph->particle_box.k + rand()/ double(RAND_MAX+1.0)*1.e-3 + sph->particle_box_l.k;
        coord.k = AMAX1(sph->particle_box_l.k+1.e-20,AMIN1(coord.k, sph->particle_box_r.k-1.e-20));
    #endif
        Calculate_particle_infor(sph);
        if (phi > h) positive_phase = true;
      }
  
  type = REAL_PARTICLE;
}
//-----------------------------------------------------------------
// Get_kernel_summation
//-----------------------------------------------------------------
void Particle_mesh_generation::Get_kernel_summation(SOLVER *sph)
{
  sum = 0.;
  Real angle_thres = sph->angle_thres;

  int n_neighbor = neighbor.size();
  int ncount = 0;
  for (int i = 0; i < n_neighbor; i++){
    Particle *neighbor_particle = neighbor[i];
    my_real neighbor_coord = neighbor_particle->coord;
    my_real neighbor_norm  = neighbor_particle->norm;
    int     neighbor_type  = neighbor_particle->type;
    
    if (type == REAL_PARTICLE){
      my_real dr = my_minus_data(coord, neighbor_coord);
      Real  dist = get_distance(dr);
      Real    w  = Kernel_function(dist, h);
            sum += w;
      ncount++;
    }else if (type == SURFACE_PARTICLE){
      // JZ20190131::project particle to tanj. plane
      my_real vec_ij  = my_minus_data (coord, neighbor_coord);
      Real dist_norm  = get_dot_product (vec_ij, norm);
      my_real norm_ij = my_multiply_const (norm, dist_norm);
      my_real      dr = my_minus_data (vec_ij, norm_ij);
      Real       dist = get_distance(dr);

      if (neighbor_type == SURFACE_PARTICLE){
        Real cos_alpha = get_dot_product (norm, neighbor_norm);
        if (fabs(cos_alpha) > cos(angle_thres)){
          #if DIM == 3
          Real w = Kernel_function_2D(dist, h);
          #elif DIM == 2
          Real w = Kernel_function_1D(dist, h);
          #endif
          sum += w;
          ncount++;
        }
      }else if (neighbor_type == SEGMENT_PARTICLE || neighbor_type == SINGULARITY_PARTICLE){
        #if DIM == 3
        Real w = Kernel_function_2D(dist, h);
        #elif DIM == 2
        Real w = Kernel_function_1D(dist, h);
        #endif
        sum += w;
        ncount++;
      }
    }else if (type == SEGMENT_PARTICLE){
      // JZ20190131::project particle to tanj. plane
      my_real vec_ij  = my_minus_data (coord, neighbor_coord);
      Real dist_norm  = get_dot_product (vec_ij, norm);
      my_real norm_ij = my_multiply_const (norm, dist_norm);
      my_real      dr = my_minus_data (vec_ij, norm_ij);
      Real       dist = get_distance(dr);

      if (neighbor_type == SEGMENT_PARTICLE){
        Real cos_alpha = get_dot_product (norm, neighbor_norm);
        if (fabs(cos_alpha) > cos(angle_thres)){
          #if DIM == 3
          Real w = Kernel_function_1D(dist, h);
          sum += w;
          ncount++;
          #endif
        }
      }else if (neighbor_type == SINGULARITY_PARTICLE){
        #if DIM == 3
        Real w = Kernel_function_1D(dist, h);
        sum += w;
        ncount++;
        #endif
      }
    }
  }
  if (ncount > 1)
    sum = 1./(sum + 1.e-15);
}
//-----------------------------------------------------------------
// Accumulate_particle_force for segment particle
//-----------------------------------------------------------------
void Particle_mesh_generation::Accumulate_particle_force_seg(SOLVER *sph, int stage)
{
  int n_neighbor = neighbor.size();
  int current_stage = sph->current_stage;
  Real angle_thres = sph->angle_thres;

  Real sigma_ = 1.;
  #ifdef _MPI_
  sigma_ = my_get_ramping_factor (sph->sigma_coeff_start, 
				  sph->sigma_coeff_end, 
				  sph->run_time, 
				  sph->sigma_ramping_start, 
				  sph->sigma_ramping_end);
  #endif
  
  for (int i = 0; i < n_neighbor; i++){
    Particle *neighbor_particle = neighbor[i];

    int type_j = neighbor_particle->type;
    if (type_j == REAL_PARTICLE || type_j == SURFACE_PARTICLE) continue;
    
    if (id != neighbor_particle->id){
      // load particle j data
      Real h_j     = neighbor_particle->h;
      Real mass_j  = neighbor_particle->mass;
      Real P_j     = neighbor_particle->P;
      Real rho_j   = neighbor_particle->rho;
      Real nu_j    = neighbor_particle->nu;
      int  color_j = neighbor_particle->color;
      
      my_real coord_j = neighbor_particle->coord;
      my_real v_j     = neighbor_particle->v;
      my_real norm_j  = neighbor_particle->norm;

      // my_real      dr = my_minus_data(coord, coord_j);
      // Real       dist = get_distance(dr);

      // JZ20190131::project particle to tanj. plane
      my_real vec_ij  = my_minus_data (coord, coord_j);
      Real dist_norm  = get_dot_product (vec_ij, norm);
      my_real norm_ij = my_multiply_const (norm, dist_norm);
      my_real      dr = my_minus_data (vec_ij, norm_ij);
      Real       dist = get_distance(dr);

      Real      vol_i = P/powern(rho,2);
      Real      vol_j = P_j/powern(rho_j,2);

      Real         aw = 0.;
      my_real    awax;

      bool skip = false;
      if (type_j == SEGMENT_PARTICLE){
        Real cos_alpha = get_dot_product (norm, norm_j);
        if (fabs(cos_alpha) <= cos(angle_thres)){
          skip = true;
        }
      }

      if (!skip){
        #if DIM == 3
        aw = Derivative_kernel_function_1D(dist, dr, 0.5*(h+h_j));
        #endif
        awax = my_multiply_const (dr, aw);
        
        // pressure force
        if (stage == 1 ){
      	  // JZ20181230::surface tension force for partitioning
      	  Real sigma = 1.;
      	  if (color_j != color) sigma = sigma_;
	  
          if (type_j == SEGMENT_PARTICLE){
            Real P_ij = sigma * mass_j * sum * (vol_i + vol_j);
            my_real dFP = my_multiply_const (awax, -P_ij);
            F1 = my_add_data(F1, dFP);              
          }else if (type_j == SINGULARITY_PARTICLE){
            Real P_ij = sigma * mass_j * sum * (vol_i + vol_j);
            my_real dFP = my_multiply_const (awax, -P_ij);
            F1 = my_add_data(F1, dFP); 
          }
        }

        // if (stage == 2 && sph->current_stage != FILL){
        if (stage == 2 ){
          Real      eta_i = nu*rho;
          Real      eta_j = nu_j*rho_j;
          Real     eta_ij = 2.*eta_j*eta_i/(eta_i+eta_j+1.e-15);
          my_real    v_ij = my_minus_data (v, v_j);
          Real      temp0 = aw*eta_ij*mass_j*(1./rho/rho + 1./rho_j/rho_j);
          my_real  dF_v  = my_multiply_const (v_ij, temp0);
                   F2  = my_add_data(F2, dF_v);
        }
      }
    }
  }
  a  = my_add_data(F1, F2);

  if (stage == 2){
    // if (sph->current_stage == FULL){
      my_real damping;

      Real damping_fac = my_get_ramping_factor (sph->damp_coeff_start, 
                                                sph->damp_coeff_end, 
                                                sph->run_time, 
                                                sph->damp_ramping_start, 
                                                sph->damp_ramping_end);
      damping = my_multiply_const (v, damping_fac);

      a = my_add_data(a, damping);
    // }
  }
}
//-----------------------------------------------------------------
// Accumulate_particle_force for surface particle
//-----------------------------------------------------------------
void Particle_mesh_generation::Accumulate_particle_force_surf(SOLVER *sph, int stage)
{
  int n_neighbor = neighbor.size();
  int current_stage = sph->current_stage;
  Real angle_thres = sph->angle_thres;

  Real sigma_ = 1.;
  #ifdef _MPI_
  sigma_ = my_get_ramping_factor (sph->sigma_coeff_start, 
				  sph->sigma_coeff_end, 
				  sph->run_time, 
				  sph->sigma_ramping_start, 
				  sph->sigma_ramping_end);
  #endif
  
  for (int i = 0; i < n_neighbor; i++){
    Particle *neighbor_particle = neighbor[i];

    int type_j = neighbor_particle->type;
    if (type_j == REAL_PARTICLE) continue;
    
    if (id != neighbor_particle->id){
      // load particle j data
      Real h_j    = neighbor_particle->h;
      Real mass_j = neighbor_particle->mass;
      Real P_j    = neighbor_particle->P;
      Real rho_j  = neighbor_particle->rho;
      Real nu_j   = neighbor_particle->nu;
      int color_j = neighbor_particle->color;
      
      my_real coord_j = neighbor_particle->coord;

      my_real v_j     = neighbor_particle->v;
      my_real norm_j  = neighbor_particle->norm;

      // my_real      dr = my_minus_data(coord, coord_j);
      // Real       dist = get_distance(dr);

      // JZ20190131::project particle to tanj. plane
      my_real vec_ij  = my_minus_data (coord, coord_j);
      Real dist_norm  = get_dot_product (vec_ij, norm);
      my_real norm_ij = my_multiply_const (norm, dist_norm);
      my_real      dr = my_minus_data (vec_ij, norm_ij);
      Real       dist = get_distance(dr);

      Real      vol_i = P/powern(rho,2);
      Real      vol_j = P_j/powern(rho_j,2);

      Real         aw = 0.;
      my_real    awax;

      bool skip = false;
      if (type_j == SURFACE_PARTICLE){
        Real cos_alpha = get_dot_product (norm, norm_j);
        if (fabs(cos_alpha) <= cos(angle_thres)){
          skip = true;
        }
      }

      if (!skip){
        #if DIM == 3
        aw = Derivative_kernel_function_2D(dist, dr, 0.5*(h+h_j));
        #elif DIM ==2
        aw = Derivative_kernel_function_1D(dist, dr, 0.5*(h+h_j));
        #endif
        awax = my_multiply_const (dr, aw);
        
        // pressure force
        // if (stage == 1 || (stage == 2 && sph->current_stage == FILL)){
        if (stage == 1 ){
          // JZ20181230::surface tension force for partitioning
          Real sigma = 1.;
          if (color_j != color) sigma = sigma_;
	  
          Real P_ij = sigma * mass_j * sum * (vol_i + vol_j);
          my_real dFP = my_multiply_const (awax, -P_ij);
          F1 = my_add_data(F1, dFP);
        }

        // if (stage == 2 && sph->current_stage != FILL){
        if (stage == 2 ){
          Real      eta_i = nu*rho;
          Real      eta_j = nu_j*rho_j;
          Real     eta_ij = 2.*eta_j*eta_i/(eta_i+eta_j+1.e-15);
          my_real    v_ij = my_minus_data (v, v_j);
          Real      temp0 = aw*eta_ij*mass_j*(1./rho/rho + 1./rho_j/rho_j);
          my_real  dF_v  = my_multiply_const (v_ij, temp0);
                   F2  = my_add_data(F2, dF_v);
        }
      }
    }
  }
  a  = my_add_data(F1, F2);

  if (stage == 2){
    // if (sph->current_stage == FULL){
      my_real damping;

      Real damping_fac = my_get_ramping_factor (sph->damp_coeff_start, 
                                                sph->damp_coeff_end, 
                                                sph->run_time, 
                                                sph->damp_ramping_start, 
                                                sph->damp_ramping_end);
      damping = my_multiply_const (v, damping_fac);

      a = my_add_data(a, damping);
    // }
  }
}
//-----------------------------------------------------------------
// Accumulate_particle_force for real particle
//-----------------------------------------------------------------
void Particle_mesh_generation::Accumulate_particle_force_real(SOLVER *sph, int stage)
{
  int n_neighbor = neighbor.size();
  int current_stage = sph->current_stage;
  Real angle_thres = sph->angle_thres;

  Real sigma_ = 1.;
  #ifdef _MPI_
  sigma_ = my_get_ramping_factor (sph->sigma_coeff_start, 
				  sph->sigma_coeff_end, 
				  sph->run_time, 
				  sph->sigma_ramping_start, 
				  sph->sigma_ramping_end);
  #endif

  #ifdef _APD_
  Real         r_avg = 0.;
  int         ncount = 0;
  #endif
  for (int i = 0; i < n_neighbor; i++){
    Particle *neighbor_particle = neighbor[i];

    int type_j = neighbor_particle->type;
    
    if (id != neighbor_particle->id){
      // load particle j data
      Real h_j    = neighbor_particle->h;
      Real mass_j = neighbor_particle->mass;
      Real P_j    = neighbor_particle->P;
      Real rho_j  = neighbor_particle->rho;
      Real nu_j   = neighbor_particle->nu;
      int color_j = neighbor_particle->color;
      
      my_real coord_j = neighbor_particle->coord;
      my_real v_j     = neighbor_particle->v;
      if (type_j != REAL_PARTICLE)
        v_j = my_multiply_const (v, -1.0);

      my_real            dr = my_minus_data(coord, coord_j);
      Real             dist = get_distance(dr);
      my_real           dr_ = my_multiply_const(dr,1./(dist+1.e-20));
      Real             h_ij = 0.5*(h+h_j);

      Real             vol_i = P/powern(rho,2);
      Real             vol_j = P_j/powern(rho_j,2);
      
      Real                w = Kernel_function(dist, h);
      Real               aw = Derivative_kernel_function(dist, dr, 0.5*(h+h_j));
      my_real          awax = my_multiply_const (dr, aw);
      
      // if (stage == 1 || (stage == 2 && sph->current_stage == FILL)){
      if (stage == 1 ){
         // pressure force
	  // JZ20181230::surface tension force for partitioning
	  Real sigma = 1.;
	  if (color_j != color) sigma = sigma_;
	  
          Real P_ij = sigma * mass_j * sum * (vol_i + vol_j);
          my_real dFP = my_multiply_const (awax, -P_ij);
          F1 = my_add_data(F1, dFP);

        // if (type_j == REAL_PARTICLE){
        //   Real P_ij = mass_j * sum * (vol_i + vol_j);
        //   my_real dFP = my_multiply_const (awax, -P_ij);
        //   F1 = my_add_data(F1, dFP);
        // }else{
        //         // v_j = my_multiply_const (v, -1);
        //   my_real norm_j = neighbor_particle->norm;
        //   Real P_ij = mass_j / h_j * sum * w *(vol_i + vol_j);
        //   my_real dFP = my_multiply_const (norm_j, P_ij);
        //   F1 = my_add_data(F1, dFP);
        // }
      }
      // if (stage == 2 && sph->current_stage != FILL){
      if (stage == 2 ){
        Real      eta_i = nu*rho;
        Real      eta_j = nu_j*rho_j;
        Real     eta_ij = 2.*eta_j*eta_i/(eta_i+eta_j+1.e-15);
        my_real    v_ij = my_minus_data (v, v_j);
        Real      temp0 = aw*eta_ij*mass_j*(1./rho/rho + 1./rho_j/rho_j);
        my_real  dF_v  = my_multiply_const (v_ij, temp0);
                 F2  = my_add_data(F2, dF_v);
      }
      #ifdef _WVT_
      if (stage == 1 || (stage == 2 && sph->current_stage == FILL)){
        Real F = pow (h_ij/(dist+pow(0.1*h_ij, DIM)), 2.);
        if (dist > CUT_OFF*h_ij) F = 0.;
        if (dist < 1.e-15) F = 0.;
        my_real ddx = my_multiply_const (dr_, F*sph->WVT_coeff*h);
        coord_corr = my_add_data(coord_corr, ddx);
      }
      #endif
      #ifdef _APD_
      if (stage == 1 || (stage == 2 && sph->current_stage == FULL)){
        r_avg += dist; ncount++;
        coord_corr = my_add_data (coord_corr, my_multiply_const(dr, 1./(pow(dist, 3) + 1.e-15)));
      }
      #endif
    }
  }
  
  a  = my_add_data(F1, F2);
  
  #ifdef _APD_
  if (stage == 1 || (stage == 2 && sph->current_stage == FULL)){
    r_avg /= (Real(ncount)+1.e-15);
    coord_corr = my_multiply_const (coord_corr, sph->APD_coeff*r_avg*r_avg*sph->glbl_max_v*sph->glbl_timestep);
  }
  #endif
  if (stage == 2){
    // if (sph->current_stage == FULL){
      my_real damping;

      Real damping_fac = my_get_ramping_factor (sph->damp_coeff_start, 
                                                sph->damp_coeff_end, 
                                                sph->run_time, 
                                                sph->damp_ramping_start, 
                                                sph->damp_ramping_end);
      damping = my_multiply_const (v, damping_fac);

      a = my_add_data(a, damping);
    // }
  }
}
//-----------------------------------------------------------------
// Map_surface_particle_force
//-----------------------------------------------------------------
void Particle_mesh_generation::Map_surface_particle_force(SOLVER *sph)
{
  Real normal_force = get_dot_product (a, norm);
  a = my_minus_data (a, my_multiply_const (norm, normal_force));
}
//-----------------------------------------------------------------
// Update_velocity_half_step
//-----------------------------------------------------------------
void Particle_mesh_generation::Update_velocity_half_step(SOLVER *sph)
{
  Real        dt = sph->glbl_timestep;
  
  v = my_add_data (v, my_multiply_const(a, 0.5*dt));
}
//-----------------------------------------------------------------
// Update_coord_full_step
//-----------------------------------------------------------------
void Particle_mesh_generation::Update_coord_full_step(SOLVER *sph)
{
  Real        dt = sph->glbl_timestep;
  my_real  box_l = sph->box_l;
  my_real  box_r = sph->box_r;
  my_real domain = sph->domain;
  
  coord = my_add_data (coord, my_multiply_const(v, dt));

  #if defined (_WVT_) || defined (_APD_)
  coord = my_add_data (coord, coord_corr);
  #endif
#ifndef _MPI_
  #if PERI_DIM_X
  if ( coord.i > box_r.i + 1.e-10) coord.i -= domain.i;
  else if ( coord.i < box_l.i - 1.e-10) coord.i += domain.i;
  #endif
  #if PERI_DIM_Y
  if ( coord.j > box_r.j + 1.e-10) coord.j -= domain.j;
  else if ( coord.j < box_l.j - 1.e-10) coord.j += domain.j;
  #endif
  #if PERI_DIM_Z
  if ( coord.k > box_r.k + 1.e-10) coord.k -= domain.k;
  else if ( coord.k < box_l.k - 1.e-10) coord.k += domain.k;
  #endif
#endif
#ifdef _MPI_
  #if P_DIM_X == 0 && PERI_DIM_X == 1
  if ( coord.i > box_r.i + 1.e-10) coord.i -= domain.i;
  else if ( coord.i < box_l.i - 1.e-10) coord.i += domain.i;
  #endif
  #if P_DIM_Y == 0 && PERI_DIM_Y == 1
  if ( coord.j > box_r.j + 1.e-10) coord.j -= domain.j;
  else if ( coord.j < box_l.j - 1.e-10) coord.j += domain.j;
  #endif
  #if P_DIM_Z == 0 && PERI_DIM_Z == 1
  if ( coord.k > box_r.k + 1.e-10) coord.k -= domain.k;
  else if ( coord.k < box_l.k - 1.e-10) coord.k += domain.k;
  #endif
#endif
}
//-----------------------------------------------------------------
// construct local voronoi diagram
//-----------------------------------------------------------------
void Particle_mesh_generation::Local_voronoi_diagram(SOLVER *sph)
{
  Voronoi_construction voro_cons;
  voro_cons.Initialize();

  voro_cons.Local_VD_generation (this, &sph->mesh);
}
//-----------------------------------------------------------------
// Reconstruce_dummy_particles
//-----------------------------------------------------------------
bool Particle_mesh_generation::Reconstruce_dummy_particles(SOLVER *sph, p_Particle dummy_particle)
{
  Real     n_d = 0.;
  Real     div = 1.0/(get_distance(norm) + 1.0e-15);
  my_real  box_l = sph->level_set.box_l;
  my_real  box_r = sph->level_set.box_r;
  Real     padding = sph->level_set.lset_level_info[0]->scale;

  srand((unsigned)time(NULL)); 
  Real   rr =  rand()/ double(RAND_MAX+1.0);
  Real      dx = (3.0 + 0.01*rr) * h;

  // Real      dx = 3.0 * h;

  my_real norm_ = my_multiply_const (norm, div);
  my_real coord_dummy;
  coord_dummy.i = 0.; coord_dummy.j = 0.; coord_dummy.k = 0.;
  Real phi_dummy = 0.;

  // while(phi_dummy > -0.8*dx){
    #if (DIM_X)
    coord_dummy.i = coord.i - norm_.i * dx;
    coord_dummy.i = AMIN1(AMAX1(coord_dummy.i, box_l.i+padding), box_r.i-padding);
    #endif
    #if (DIM_Y)
    coord_dummy.j = coord.j - norm_.j * dx;
    coord_dummy.j = AMIN1(AMAX1(coord_dummy.j, box_l.j+padding), box_r.j-padding);
    #endif
    #if (DIM_Z)
    coord_dummy.k = coord.k - norm_.k * dx;
    coord_dummy.k = AMIN1(AMAX1(coord_dummy.k, box_l.k+padding), box_r.k-padding);
    #endif

    my_set_data (dummy_particle->coord, coord_dummy);

    dummy_particle->Calculate_particle_infor(sph);

    phi_dummy = dummy_particle->phi;

    if (phi_dummy <= -0.99*dx) return true;

    return false;

  //   dx *= 0.9;
  // }
}
#if SYM_DIM != 0
//-----------------------------------------------------
// set symmetric particle information
//-----------------------------------------------------
void Particle_mesh_generation::Set_bc_particle_info(p_Particle copy_particle, my_real shift, my_real pos)
{
  color    = copy_particle->color;
  id       = copy_particle->id;
  local_id = copy_particle->local_id;
  sum      = copy_particle->sum;
  h        = copy_particle->h;
  scale    = copy_particle->scale;
  mass     = copy_particle->mass;
  type     = copy_particle->type;
  level    = copy_particle->level;
  rho      = copy_particle->rho;
  P        = copy_particle->P;

  copy.clear();
  copy.push_back(copy_particle);
  Real slip = 0.;
#ifdef _NON_SLIP_
  slip = -1.;
#endif
#ifdef _FREE_SLIP_
  slip = 1.;
#endif

  v.i = pos.i == 0? slip*copy_particle->v.i : (-1.*copy_particle->v.i);
  v.j = pos.j == 0? slip*copy_particle->v.j : (-1.*copy_particle->v.j);
  v.k = pos.k == 0? slip*copy_particle->v.k : (-1.*copy_particle->v.k);

  my_set_data   (coord, shift);
  my_set_data   ( norm, copy_particle->norm);
  my_set_const  (    a, 0);
  my_set_const  (   F1, 0.0);
  my_set_const  (   F2, 0.0);
}
//-----------------------------------------------------
// refresh symmetric particle information
//-----------------------------------------------------
void Particle_mesh_generation::Refresh_sbc_particle_info(p_Particle copy_particle, int flag, SPH *sph)
{
  if (flag == 0){

  }else if (flag == 1){

  }else if (flag == 2){
    
  }
}
#endif
#if PERI_DIM != 0
//-----------------------------------------------------
// set periodical particle information
//-----------------------------------------------------
void Particle_mesh_generation::Set_bc_particle_info(p_Particle copy_particle, my_real shift)
{
  color    = copy_particle->color;
  id       = copy_particle->id;
  local_id = copy_particle->local_id;
  sum      = copy_particle->sum;
  h        = copy_particle->h;
  scale    = copy_particle->scale;
  mass     = copy_particle->mass;
  type     = copy_particle->type;
  level    = copy_particle->level;
  rho      = copy_particle->rho;
  P        = copy_particle->P;
  coord    = my_add_data(copy_particle->coord, shift);
  copy.clear();

  copy.push_back(copy_particle);
  my_set_data   (   v, copy_particle->v);
  my_set_data   (norm, copy_particle->norm);
  my_set_const  (   a, 0);
  my_set_const  (  F1, 0.0);
  my_set_const  (  F2, 0.0);
}
//-----------------------------------------------------
// refresh periodical particle information
//-----------------------------------------------------
void Particle_mesh_generation::Refresh_pbc_particle_info(p_Particle copy_particle, int flag)
{
  if (flag == 0){

  }else if (flag == 1){

  }else if (flag == 2){

  }
}
#endif
#ifdef _MPI_
//-----------------------------------------------------
// Set buffer particle information from incoming particle
//-----------------------------------------------------
void Particle_mesh_generation::Set_buffer_particle_info(p_Particle temp, int flag)
{
  if (flag == 1){
    sum = temp->sum;
  }else if (flag == 2){
    nu  = temp->nu;
    my_set_data (v, temp->v);
  }else if (flag == 4){

  }
}
//-----------------------------------------------------
// Set particle information from incoming particle
//-----------------------------------------------------
void Particle_mesh_generation::Set_particle_info(p_Particle income_particle)
{
  color      = income_particle->color;
  id         = income_particle->id;
  type       = income_particle->type;
  p_mass     = income_particle->p_mass;
  scale      = income_particle->scale;
  h          = income_particle->h;
  level      = income_particle->level;
  rho        = income_particle->rho;
  vol        = income_particle->vol;
  mass       = income_particle->mass;
  P          = income_particle->P;
  safe_guard = income_particle->safe_guard;
  
  my_set_data(coord, income_particle->coord);
}
#endif
