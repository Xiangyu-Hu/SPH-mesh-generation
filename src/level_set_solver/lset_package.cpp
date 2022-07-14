#include "glbfunc.h"
#include "level_set.h"
#include "lset_level_infor.h"
#include "lset_package.h"

typedef boost::multi_array_types::extent_range mrange;

/***************************************************/
/*                                                 */
/*   Functions defined in class "Levelset_package" */
/*                                                 */
/***************************************************/
//------------------------------------------------------------
// initialze all the necessary parameters 
// for Levelset_package
//------------------------------------------------------------
void Levelset_package::Initialize(p_Levelset_levelinfo lset_level_info, int pkg_i, int pkg_j, int pkg_k)
{
  level = lset_level_info->level;
  index.i = pkg_i;
  index.j = pkg_j;
  index.k = pkg_k;
  
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        lset_cell[i][j][k].Initialize(this);
      }

  total_mass = 0.;
  total_mass_surface = 0.;
  total_mass_segment = 0.;
  total_volume = 0.;
  total_area = 0.;
  total_length = 0.;
  maximum_phi = 0.;
  maximum_dl = 0.;
  minimum_dl = 0.;
  middel_dl = 0.;
  maximum_curv = 0.;
  minimum_curv = 0.;
  maximum_psi = 0.;
  tag_convergence = 0;
  n_interface_cell = 0;

#ifdef _MPI_
  color = -1;
#endif
}
//------------------------------------------------------------
// initialze levelset phi
//------------------------------------------------------------
void Levelset_package::Define_levelset_full_field(int pkg_i, int pkg_j, int pkg_k, my_real dpkg, my_real dcell, Levelset *level_set)
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        my_real pos = get_cell_position (pkg_i, pkg_j, pkg_k, i, j, k, dpkg, dcell, level_set->box_l);
        Real phi = level_set->F_phi (pos);
        
        lset_cell[i][j][k].Set_phi (phi);
      }
}
//------------------------------------------------------------
// get_cell_position
//------------------------------------------------------------
my_real Levelset_package::get_cell_position (int pkg_i, int pkg_j, int pkg_k, int cell_i, int cell_j, int cell_k, my_real dpkg, my_real dcell, my_real box_l)
{
  my_real pos_;
  
  pos_.i = DIM_X==1 ? pkg_i*dpkg.i + ( cell_i + 0.5)*dcell.i + box_l.i : 0.;
  pos_.j = DIM_Y==1 ? pkg_j*dpkg.j + ( cell_j + 0.5)*dcell.j + box_l.j : 0.;
  pos_.k = DIM_Z==1 ? pkg_k*dpkg.k + ( cell_k + 0.5)*dcell.k + box_l.k : 0.;
  
  return pos_;
}
//------------------------------------------------------------
// get_cell_position corner point
//------------------------------------------------------------
my_real Levelset_package::get_cell_position_corner (int pkg_i, int pkg_j, int pkg_k, int cell_i, int cell_j, int cell_k, my_real dpkg, my_real dcell, my_real box_l)
{
  my_real pos_;
  
  pos_.i = DIM_X==1 ? pkg_i*dpkg.i + cell_i*dcell.i + box_l.i : 0.;
  pos_.j = DIM_Y==1 ? pkg_j*dpkg.j + cell_j*dcell.j + box_l.j : 0.;
  pos_.k = DIM_Z==1 ? pkg_k*dpkg.k + cell_k*dcell.k + box_l.k : 0.;
  
  return pos_;
}
//------------------------------------------------------------
// get_pkg_position_corner corner point
//------------------------------------------------------------
my_real Levelset_package::get_pkg_position_corner (int pkg_i, int pkg_j, int pkg_k, my_real dpkg, my_real box_l)
{
  my_real pos_;
  
  pos_.i = DIM_X==1 ? pkg_i*dpkg.i + box_l.i : 0.;
  pos_.j = DIM_Y==1 ? pkg_j*dpkg.j + box_l.j : 0.;
  pos_.k = DIM_Z==1 ? pkg_k*dpkg.k + box_l.k : 0.;
  
  return pos_;
}
//------------------------------------------------------------
// get_pkg_position_corner corner point
//------------------------------------------------------------
my_real Levelset_package::get_pkg_position (int pkg_i, int pkg_j, int pkg_k, my_real dpkg, my_real box_l)
{
  my_real pos_;
  
  pos_.i = DIM_X==1 ? (Real(pkg_i)+0.5)*dpkg.i + box_l.i : 0.;
  pos_.j = DIM_Y==1 ? (Real(pkg_j)+0.5)*dpkg.j + box_l.j : 0.;
  pos_.k = DIM_Z==1 ? (Real(pkg_k)+0.5)*dpkg.k + box_l.k : 0.;
  
  return pos_;
}
//------------------------------------------------------------
// Copy_phi_value
//------------------------------------------------------------
void Levelset_package::Copy_phi_value()
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        lset_cell[i][j][k].Set_phi_1();
      }
}
//------------------------------------------------------------
// Set_cell_topology
//------------------------------------------------------------
void Levelset_package::Set_cell_topology(p_Levelset_levelinfo lset_level_info)
{
  p_cell.resize(boost::extents[mrange(-BwidthX,ICPX+BwidthX)]
                              [mrange(-BwidthY,ICPY+BwidthY)]
                              [mrange(-BwidthZ,ICPZ+BwidthZ)]);
                           
  int     i_c = DIM_X==1 ? 1 : 0;
  int     j_c = DIM_Y==1 ? 1 : 0;
  int     k_c = DIM_Z==1 ? 1 : 0;
  

  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        p_cell[i][j][k] = & (lset_cell[i][j][k]);
      }
  
  my_int pkg_start = lset_level_info->pkg_start;
  my_int pkg_end = lset_level_info->pkg_end;
#if DIM_X
  for(int t=AMAX1((index.i-i_c),pkg_start.i); t<=AMIN1((index.i+i_c),(pkg_end.i-1)); t+=2){
    p_Levelset_package pkg_j = lset_level_info->table_lset_pkg_list[t][index.j][index.k];
    int tt = t - index.i;
    int offset = (tt == 1) ? ICPX-1 : 0;
    
    for (int i = 0; i < BwidthX; i++){
      for (int j = 0; j < ICPY; j++){
        for (int k = 0; k < ICPZ; k++){
          int i_ = tt*(i+1)+offset;
          int i__ = i_ - ICPX*tt;
          
          p_cell[i_][j][k] = &pkg_j->lset_cell[i__][j][k];
        }
      }
    }
  }
#endif
#if DIM_Y
  for(int s=AMAX1((index.j-j_c),pkg_start.j); s<=AMIN1((index.j+j_c),(pkg_end.j-1)); s+=2){
    p_Levelset_package pkg_j = lset_level_info->table_lset_pkg_list[index.i][s][index.k];
    int ss = s - index.j;
    int offset = (ss == 1) ? ICPY-1 : 0;
    
    for (int i = 0; i < ICPX; i++){
      for (int j = 0; j < BwidthY; j++){
        for (int k = 0; k < ICPZ; k++){
          int j_ = ss*(j+1)+offset;
          int j__ = j_ - ICPY*ss;
          
          p_cell[i][j_][k] = &pkg_j->lset_cell[i][j__][k];
        }
      }
    }
  }
#endif
#if DIM_Z
  for(int m=AMAX1((index.k-k_c),pkg_start.k);  m<=AMIN1((index.k+k_c),(pkg_end.k-1)); m+=2){
    p_Levelset_package pkg_j = lset_level_info->table_lset_pkg_list[index.i][index.j][m];
    int mm = m - index.k;
    int offset = (mm == 1) ? ICPZ-1 : 0;
    
    for (int i = 0; i < ICPX; i++){
      for (int j = 0; j < ICPY; j++){
        for (int k = 0; k < BwidthZ; k++){
          int k_ = mm*(k+1)+offset;
          int k__ = k_ - ICPZ*mm;
          
          p_cell[i][j][k_] = &pkg_j->lset_cell[i][j][k__];
        }
      }
    }
  }
#endif
#if DIM_X*DIM_Y
  for(int t=AMAX1((index.i-i_c),pkg_start.i); t<=AMIN1((index.i+i_c),(pkg_end.i-1)); t+=2){
    for(int s=AMAX1((index.j-j_c),pkg_start.j); s<=AMIN1((index.j+j_c),(pkg_end.j-1)); s+=2){
      p_Levelset_package pkg_j = lset_level_info->table_lset_pkg_list[t][s][index.k];
      int tt = t - index.i;
      int offset_x = (tt == 1) ? ICPX-1 : 0;
      
      int ss = s - index.j;
      int offset_y = (ss == 1) ? ICPY-1 : 0;
    
      for (int i = 0; i < BwidthX; i++){
        for (int j = 0; j < BwidthY; j++){
          for (int k = 0; k < ICPZ; k++){
            int i_ = tt*(i+1)+offset_x;
            int i__ = i_ - ICPX*tt;
            int j_ = ss*(j+1)+offset_y;
            int j__ = j_ - ICPY*ss;
          
            p_cell[i_][j_][k] = &pkg_j->lset_cell[i__][j__][k];
          }
        }
      }
    }
  }
#endif
#if DIM_X*DIM_Z
  for(int t=AMAX1((index.i-i_c),pkg_start.i); t<=AMIN1((index.i+i_c),(pkg_end.i-1)); t+=2){
    for(int m=AMAX1((index.k-k_c),pkg_start.k);  m<=AMIN1((index.k+k_c),(pkg_end.k-1)); m+=2){
      p_Levelset_package pkg_j = lset_level_info->table_lset_pkg_list[t][index.j][m];
      int tt = t - index.i;
      int offset_x = (tt == 1) ? ICPX-1 : 0;
      
      int mm = m - index.k;
      int offset_z = (mm == 1) ? ICPZ-1 : 0;
    
      for (int i = 0; i < BwidthX; i++){
        for (int j = 0; j < ICPY; j++){
          for (int k = 0; k < BwidthZ; k++){
            int i_ = tt*(i+1)+offset_x;
            int i__ = i_ - ICPX*tt;
            int k_ = mm*(k+1)+offset_z;
            int k__ = k_ - ICPZ*mm;
          
            p_cell[i_][j][k_] = &pkg_j->lset_cell[i__][j][k__];
          }
        }
      }
    }
  }
#endif
#if DIM_Y*DIM_Z
  for(int s=AMAX1((index.j-j_c),pkg_start.j); s<=AMIN1((index.j+j_c),(pkg_end.j-1)); s+=2){
    for(int m=AMAX1((index.k-k_c),pkg_start.k);  m<=AMIN1((index.k+k_c),(pkg_end.k-1)); m+=2){
      p_Levelset_package pkg_j = lset_level_info->table_lset_pkg_list[index.i][s][m];
      int ss = s - index.j;
      int offset_y = (ss == 1) ? ICPY-1 : 0;
      
      int mm = m - index.k;
      int offset_z = (mm == 1) ? ICPZ-1 : 0;
    
      for (int i = 0; i < ICPX; i++){
        for (int j = 0; j < BwidthY; j++){
          for (int k = 0; k < BwidthZ; k++){
            int j_ = ss*(j+1)+offset_y;
            int j__ = j_ - ICPY*ss;
            int k_ = mm*(k+1)+offset_z;
            int k__ = k_ - ICPZ*mm;
          
            p_cell[i][j_][k_] = &pkg_j->lset_cell[i][j__][k__];
          }
        }
      }
    }
  }
#endif
#if DIM_X*DIM_Y*DIM_Z
  for(int t=AMAX1((index.i-i_c),pkg_start.i); t<=AMIN1((index.i+i_c),(pkg_end.i-1)); t+=2){
    for(int s=AMAX1((index.j-j_c),pkg_start.j); s<=AMIN1((index.j+j_c),(pkg_end.j-1)); s+=2){
      for(int m=AMAX1((index.k-k_c),pkg_start.k);  m<=AMIN1((index.k+k_c),(pkg_end.k-1)); m+=2){
        p_Levelset_package pkg_j = lset_level_info->table_lset_pkg_list[t][s][m];
        int tt = t - index.i;
        int offset_x = (tt == 1) ? ICPX-1 : 0;
        
        int ss = s - index.j;
        int offset_y = (ss == 1) ? ICPY-1 : 0;
      
        int mm = m - index.k;
        int offset_z = (mm == 1) ? ICPZ-1 : 0;
    
        for (int i = 0; i < BwidthX; i++){
          for (int j = 0; j < BwidthY; j++){
            for (int k = 0; k < BwidthZ; k++){
              int i_ = tt*(i+1)+offset_x;
              int i__ = i_ - ICPX*tt;
              int j_ = ss*(j+1)+offset_y;
              int j__ = j_ - ICPY*ss;
              int k_ = mm*(k+1)+offset_z;
              int k__ = k_ - ICPZ*mm;
          
              p_cell[i_][j_][k_] = &pkg_j->lset_cell[i__][j__][k__];
            }
          }
        }
      }
    }
  }
#endif
}
//------------------------------------------------------------
// Boundary
//------------------------------------------------------------
void Levelset_package::Boundary(p_Levelset_levelinfo lset_level_info)
{
  
}
//------------------------------------------------------------
// Boundary
//------------------------------------------------------------
void Levelset_package::Boundary_curv(p_Levelset_levelinfo lset_level_info)
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        lset_cell[i][j][k].curv = 0.;
      }
}
//------------------------------------------------------------
// Get_normal
//------------------------------------------------------------
void Levelset_package::Get_normal(p_Levelset_levelinfo lset_level_info)
{
  my_real dcell = lset_level_info->dcell;
  // normal
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        my_real dPhidr; my_set_const (dPhidr, 0.);
#if DIM_X
          dPhidr.i = 0.5 * (p_cell[i+1][j][k]->phi - p_cell[i-1][j][k]->phi) / dcell.i;
#endif
#if DIM_Y
          dPhidr.j = 0.5 * (p_cell[i][j+1][k]->phi - p_cell[i][j-1][k]->phi) / dcell.j;
#endif
#if DIM_Z
          dPhidr.k = 0.5 * (p_cell[i][j][k+1]->phi - p_cell[i][j][k-1]->phi) / dcell.k;
#endif
        Real dist = get_distance (dPhidr) + 1.e-12;
        dPhidr = my_multiply_const (dPhidr, 1./dist);
        p_cell[i][j][k]->n_x = DIM_X ? dPhidr.i : 0.;
        p_cell[i][j][k]->n_y = DIM_Y ? dPhidr.j : 0.;
        p_cell[i][j][k]->n_z = DIM_Z ? dPhidr.k : 0.;
      }
}
//------------------------------------------------------------
// Get_normal_for_BC_pkg
//------------------------------------------------------------
void Levelset_package::Get_BC_normal(p_Levelset_levelinfo lset_level_info)
{
  my_real dcell = lset_level_info->dcell;
  // normal
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        my_real dPhidr; my_set_const (dPhidr, 0.);
#if DIM_X
          if (i == 0) dPhidr.i = (lset_cell[i+1][j][k].phi - lset_cell[i][j][k].phi) / dcell.i;
          else if (i == ICPX-1) dPhidr.i = (lset_cell[i][j][k].phi - lset_cell[i-1][j][k].phi) / dcell.i;
          else dPhidr.i = 0.5 * (lset_cell[i+1][j][k].phi - lset_cell[i-1][j][k].phi) / dcell.i;
          
#endif
#if DIM_Y
          if (j == 0) dPhidr.j = (lset_cell[i][j+1][k].phi - lset_cell[i][j][k].phi) / dcell.j;
          else if (j == ICPY-1) dPhidr.j = (lset_cell[i][j][k].phi - lset_cell[i][j-1][k].phi) / dcell.j;
          else dPhidr.j = 0.5 * (lset_cell[i][j+1][k].phi - lset_cell[i][j-1][k].phi) / dcell.j;
#endif
#if DIM_Z
          if (k == 0) dPhidr.k = (lset_cell[i][j][k+1].phi - lset_cell[i][j][k].phi) / dcell.k;
          else if (k == ICPZ-1) dPhidr.k = (lset_cell[i][j][k].phi - lset_cell[i][j][k-1].phi) / dcell.k;
          else dPhidr.k = 0.5 * (lset_cell[i][j][k+1].phi - lset_cell[i][j][k-1].phi) / dcell.k;
#endif
        Real dist = get_distance (dPhidr) + 1.e-12;
        dPhidr = my_multiply_const (dPhidr, 1./dist);
        lset_cell[i][j][k].n_x = DIM_X ? dPhidr.i : 0.;
        lset_cell[i][j][k].n_y = DIM_Y ? dPhidr.j : 0.;
        lset_cell[i][j][k].n_z = DIM_Z ? dPhidr.k : 0.;
      }
}
//------------------------------------------------------------
// Get_curvature
//------------------------------------------------------------
void Levelset_package::Get_curvature(p_Levelset_levelinfo lset_level_info)
{
  my_real dcell = lset_level_info->dcell;
  
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        if (p_cell[i][j][k]->tag_interface == CUT_CELL){
          Real phi_x, phi_y, phi_z, phi_xx, phi_yy, phi_zz, phi_xy, phi_yz, phi_xz;
          phi_x = phi_y = phi_z = phi_xx = phi_yy = phi_zz = phi_xy = phi_yz = phi_xz = 0.;
          #if DIM_X
          phi_x = (8.0*(p_cell[i+1][j][k]->phi -   p_cell[i-1][j][k]->phi) - (p_cell[i+2][j][k]->phi - p_cell[i-2][j][k]->phi))/dcell.i/12.0;
          phi_xx =     (p_cell[i+1][j][k]->phi - 2*p_cell[i  ][j][k]->phi  +  p_cell[i-1][j][k]->phi)/dcell.i/dcell.i;
          #endif
          #if DIM_Y
          phi_y = (8.0*(p_cell[i][j+1][k]->phi -   p_cell[i][j-1][k]->phi) - (p_cell[i][j+2][k]->phi - p_cell[i][j-2][k]->phi))/dcell.j/12.0;
          phi_yy =     (p_cell[i][j+1][k]->phi - 2*p_cell[i][j  ][k]->phi  +  p_cell[i][j-1][k]->phi)/dcell.j/dcell.j;
          #endif
          #if DIM_Z
          phi_z = (8.0*(p_cell[i][j][k+1]->phi -   p_cell[i][j][k-1]->phi) - (p_cell[i][j][k+2]->phi - p_cell[i][j][k-2]->phi))/dcell.k/12.0;
          phi_zz =     (p_cell[i][j][k+1]->phi - 2*p_cell[i][j][k  ]->phi  +  p_cell[i][j][k-1]->phi)/dcell.k/dcell.k;
          #endif
          #if DIM_X*DIM_Y
          phi_xy = 0.25*(p_cell[i+1][j+1][k]->phi-p_cell[i-1][j+1][k]->phi-p_cell[i+1][j-1][k]->phi+p_cell[i-1][j-1][k]->phi)/dcell.i/dcell.j;
          #endif
          #if DIM_X*DIM_Z
          phi_xz = 0.25*(p_cell[i+1][j][k+1]->phi-p_cell[i-1][j][k+1]->phi-p_cell[i+1][j][k-1]->phi+p_cell[i-1][j][k-1]->phi)/dcell.i/dcell.k;
          #endif
          #if DIM_Y*DIM_Z
          phi_yz = 0.25*(p_cell[i][j+1][k+1]->phi-p_cell[i][j+1][k-1]->phi-p_cell[i][j-1][k+1]->phi+p_cell[i][j-1][k-1]->phi)/dcell.j/dcell.k;
          #endif
          Real curv = ( phi_x*phi_x*phi_yy  - 2*phi_x*phi_y*phi_xy + phi_y*phi_y*phi_xx
                      + phi_x*phi_x*phi_zz  - 2*phi_x*phi_z*phi_xz + phi_z*phi_z*phi_xx
                      + phi_y*phi_y*phi_zz  - 2*phi_y*phi_z*phi_yz + phi_z*phi_z*phi_yy )
                      / pow(phi_x*phi_x + phi_y*phi_y + phi_z*phi_z, (Real)1.5);

          Real dist = get_distance (dcell);
          p_cell[i][j][k]->curv = fabs(AMIN1(AMAX1(curv, -2./dist), 2./dist));
        }
      }
}
//------------------------------------------------------------
// Clean_curvature
//------------------------------------------------------------
void Levelset_package::Clean_curvature()
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        if (p_cell[i][j][k]->tag_interface == NORMAL_CELL){
          p_cell[i][j][k]->curv = 0.;
        }
      }
}
//------------------------------------------------------------
// tag interface
//------------------------------------------------------------
void Levelset_package::Get_interface_tag(p_Levelset_levelinfo lset_level_info)
{
  my_real dcell = lset_level_info->dcell;
  n_interface_cell = 0;
  
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        p_cell[i][j][k]->tag_interface = NORMAL_CELL;
      }
      
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        if (Is_cut_cell(i, j, k)){
          p_cell[i][j][k]->tag_interface = CUT_CELL;
          n_interface_cell ++;
        }
      }
}
//------------------------------------------------------------
// tag interface
//------------------------------------------------------------
void Levelset_package::Get_extended_interface_tag(p_Levelset_levelinfo lset_level_info)
{
  my_real dcell = lset_level_info->dcell;

  int narrow_band_small_i = DIM_X ? NARROW_BAND_SMALL: 0;
  int narrow_band_small_j = DIM_Y ? NARROW_BAND_SMALL: 0;
  int narrow_band_small_k = DIM_Z ? NARROW_BAND_SMALL: 0;

  int narrow_band_large_i = DIM_X ? NARROW_BAND_LARGE: 0;
  int narrow_band_large_j = DIM_Y ? NARROW_BAND_LARGE: 0;
  int narrow_band_large_k = DIM_Z ? NARROW_BAND_LARGE: 0;


  for (int i = 0; i < ICPX; i++){
    for (int j = 0; j < ICPY; j++){
      for (int k = 0; k < ICPZ; k++){
        for(int r = -narrow_band_small_i; r < narrow_band_small_i+1; ++r){
          for(int s= -narrow_band_small_j; s < narrow_band_small_j+1; ++s){
            for(int t= -narrow_band_small_k; t < narrow_band_small_k+1; ++t){
              if(    p_cell[i+r][j+s][k+t]->tag_interface == CUT_CELL 
                  && p_cell[i][j][k]->tag_interface != CUT_CELL)
                p_cell[i][j][k]->tag_interface = NARROW_BAND_SMALL;
            }
          }
        }
        for(int r = -narrow_band_large_i; r < narrow_band_large_i+1; ++r){
          for(int s = -narrow_band_large_j; s < narrow_band_large_j+1; ++s){
            for(int t= -narrow_band_large_k; t < narrow_band_large_k+1; ++t){
              if(    p_cell[i+r][j+s][k+t]->tag_interface == CUT_CELL 
                  && p_cell[i][j][k]->tag_interface != CUT_CELL 
                  && p_cell[i][j][k]->tag_interface != NARROW_BAND_SMALL)
                p_cell[i][j][k]->tag_interface = NARROW_BAND_LARGE;
            }
          }
        }
      }
    }
  }
}
//------------------------------------------------------------
// Get_singularity_tag
//------------------------------------------------------------
void Levelset_package::Get_singularity_tag(Levelset *level_set, my_real *Pos_singularity)
{
  my_real dcell = level_set->lset_level_info[level]->dcell;
  my_real dpkg = level_set->lset_level_info[level]->dpkg;
  Real num_singularity = level_set->num_singularity;
  Real dl = level_set->lset_level_info[level]->dl;
  
  for (int ic = 0; ic < num_singularity; ic++){
    my_real singularity_i = level_set->singularity[ic];
    my_real coord_pkg = get_pkg_position_corner (index.i, index.j, index.k, dpkg, level_set->box_l);

    if (is_in_current_pkg (singularity_i, coord_pkg, dpkg)){
      std::vector<std::pair<Real, my_int>> dist_id; dist_id.clear();
      for (int i = 0; i < ICPX; i++)
        for (int j = 0; j < ICPY; j++)
          for (int k = 0; k < ICPZ; k++){
            if (p_cell[i][j][k]->tag_interface != NORMAL_CELL){
              std::pair<Real, my_int> dis_id_;
              
              my_real coord = get_cell_position (index.i, index.j, index.k, i, j, k, dpkg, dcell, level_set->box_l);
              Real dist = get_distance_2p (coord, singularity_i);
              my_int ID; 
              ID.i = i; ID.j = j; ID.k = k;
              
              dis_id_.first = dist;
              dis_id_.second = ID;
              dist_id.push_back(dis_id_);
            }
          }

      std::sort( dist_id.begin(), dist_id.end(), pairCompare);
      my_int id_min = dist_id[0].second;
      p_cell[id_min.i][id_min.j][id_min.k]->tag_characteristic = SINGULARITY_CELL;
      p_cell[id_min.i][id_min.j][id_min.k]->idx_characteristic = ic;

      my_real coord = get_cell_position (index.i, index.j, index.k, id_min.i, id_min.j, id_min.k, dpkg, dcell, level_set->box_l);
      my_set_data (Pos_singularity[ic], coord);

      dist_id.clear();
    }
  }
}
//------------------------------------------------------------
// Get_segment_tag
//------------------------------------------------------------
void Levelset_package::Get_segment_tag(Levelset *level_set)
{
  my_real    dcell = level_set->lset_level_info[level]->dcell;
  my_real     dpkg = level_set->lset_level_info[level]->dpkg;
  Real num_segment = level_set->num_segment;
  Real          dl = level_set->lset_level_info[level]->dl;
  std::pair <my_real, my_real> *segment = level_set->segment;

  for (int i = 0; i < ICPX; i++){
    for (int j = 0; j < ICPY; j++){
      for (int k = 0; k < ICPZ; k++){
        if (p_cell[i][j][k]->tag_interface != NORMAL_CELL && p_cell[i][j][k]->tag_characteristic != SINGULARITY_CELL){

          my_real min = get_cell_position_corner (index.i, index.j, index.k, i, j, k, dpkg, dcell, level_set->box_l);
          my_real max = my_add_data (min, dcell);

          for (int iseg= 0; iseg < num_segment; iseg++){
            if (is_Segment_intersec_AABB (segment[iseg].first, segment[iseg].second, min, max)){
              p_cell[i][j][k]->tag_characteristic = SEGMENT_CELL;
              p_cell[i][j][k]->idx_characteristic = iseg;
              // cout<<"find segment cell"<<endl;
              break;
            }
          }
        }
      }
    }
  }
}
//------------------------------------------------------------
// is_Segment_intersec_AABB
//------------------------------------------------------------
bool Levelset_package::is_Segment_intersec_AABB (my_real p1, my_real p2, my_real min, my_real max)
{
  Real EPSILON = 1.e-20;

  my_real d = my_multiply_const (my_minus_data(p2, p1), 0.5);
  my_real e = my_multiply_const (my_minus_data(max, min), 0.5);
  my_real c = my_minus_data(my_add_data(p1, d), my_multiply_const (my_add_data(max, min), 0.5));
  my_real ad;
  ad.i = fabs(d.i);
  ad.j = fabs(d.j);
  ad.k = fabs(d.k);

  if (fabs(c.i) > e.i + ad.i)
      return false;
  if (fabs(c.j) > e.j + ad.j)
      return false;
  if (fabs(c.k) > e.k + ad.k)
      return false;

  if (fabs(d.j * c.k - d.k * c.j) > e.j * ad.k + e.k * ad.j + EPSILON)
      return false;
  if (fabs(d.k * c.i - d.i * c.k) > e.k * ad.i + e.i * ad.k + EPSILON)
      return false;
  if (fabs(d.i * c.j - d.j * c.i) > e.i * ad.j + e.j * ad.i + EPSILON)
      return false;

  return true;
}
//------------------------------------------------------------
// is_in_current_cell
//------------------------------------------------------------
bool Levelset_package::is_in_current_cell (my_real pos, my_real cell_corner, my_real dcell)
{
  my_real dr = my_minus_data (pos, cell_corner);
  bool in_x = DIM_X ? (dr.i < dcell.i+1.e-10 && dr.i > 0.) : true;
  bool in_y = DIM_Y ? (dr.j < dcell.j+1.e-10 && dr.j > 0.) : true;
  bool in_z = DIM_Z ? (dr.k < dcell.k+1.e-10 && dr.k > 0.) : true;
  
  return (in_x && in_y && in_z);
}
//------------------------------------------------------------
// is_in_current_cell
//------------------------------------------------------------
bool Levelset_package::is_in_current_pkg (my_real pos, my_real pkg_corner, my_real dpkg)
{
  my_real dr = my_minus_data (pos, pkg_corner);
  bool in_x = DIM_X ? (dr.i < dpkg.i+1.e-10 && dr.i > 0.) : true;
  bool in_y = DIM_Y ? (dr.j < dpkg.j+1.e-10 && dr.j > 0.) : true;
  bool in_z = DIM_Z ? (dr.k < dpkg.k+1.e-10 && dr.k > 0.) : true;
  
  return (in_x && in_y && in_z);
}
//------------------------------------------------------------
// Return_package_interface_tag
//------------------------------------------------------------
int Levelset_package::Return_package_interface_tag(p_Levelset_levelinfo lset_level_info)
{
  int tmp = NORMAL_CELL;
  
  for (int i = 0; i < ICPX; i++){
    for (int j = 0; j < ICPY; j++){
      for (int k = 0; k < ICPZ; k++){
        tmp = AMAX1 (p_cell[i][j][k]->tag_interface, tmp);
      }
    }
  }
  
  return tmp;
}
//------------------------------------------------------------
// Is_cut_cell
//------------------------------------------------------------
bool Levelset_package::Is_cut_cell (int i, int j, int k)
{
  std::array<Real,8> phi_corner_values;
  for(int c = 0; c < 8; ++c){
      phi_corner_values[c] = 0.0;
  }

  Get_phi_at_corners(i,j,k, phi_corner_values);

  for(int c = 0; c < 8; ++c){
      if(signum(phi_corner_values[c]) != signum(phi_corner_values[0])) return true;
  }

  if(signum(p_cell[i][j][k]->phi) != signum(phi_corner_values[0])) return true;

  return false;
}
//------------------------------------------------------------
// GetPhiAtCorners
//------------------------------------------------------------
void Levelset_package::Get_phi_at_corners(int i, int j, int k, std::array<Real,8> &phi_corner_values)
{
  // Create variables for indices i,j,k +/- 1
  int ip = DIM_X ? i + 1 : 0;
  int im = DIM_X ? i - 1 : 0;
  int jp = DIM_Y ? j + 1 : 0;
  int jm = DIM_Y ? j - 1 : 0;
  int kp = DIM_Z ? k + 1 : 0;
  int km = DIM_Z ? k - 1 : 0;

  // Read variables for the levelset values at cell (i,j,k) and its 26 adjacent cells
  // Values need for 1D
  Real phi_im_j__k = p_cell[im][j][k]->phi;
  Real phi_i__j__k = p_cell[i] [j][k]->phi;
  Real phi_ip_j__k = p_cell[ip][j][k]->phi;

  // Additional values need for 2D
  Real phi_im_jp_k = p_cell[im][jp][k]->phi;
  Real phi_i__jp_k = p_cell[i] [jp][k]->phi;
  Real phi_ip_jp_k = p_cell[ip][jp][k]->phi;
  Real phi_im_jm_k = p_cell[im][jm][k]->phi;
  Real phi_i__jm_k = p_cell[i] [jm][k]->phi;
  Real phi_ip_jm_k = p_cell[ip][jm][k]->phi;

  // Additional values need for 3D
  Real phi_im_j__kp = p_cell[im][j] [kp]->phi;
  Real phi_i__j__kp = p_cell[i] [j] [kp]->phi;
  Real phi_ip_j__kp = p_cell[ip][j] [kp]->phi;
  Real phi_im_jp_kp = p_cell[im][jp][kp]->phi;
  Real phi_i__jp_kp = p_cell[i] [jp][kp]->phi;
  Real phi_ip_jp_kp = p_cell[ip][jp][kp]->phi;
  Real phi_im_jm_kp = p_cell[im][jm][kp]->phi;
  Real phi_i__jm_kp = p_cell[i] [jm][kp]->phi;
  Real phi_ip_jm_kp = p_cell[ip][jm][kp]->phi;

  Real phi_im_j__km = p_cell[im][j] [km]->phi;
  Real phi_i__j__km = p_cell[i] [j] [km]->phi;
  Real phi_ip_j__km = p_cell[ip][j] [km]->phi;
  Real phi_im_jp_km = p_cell[im][jp][km]->phi;
  Real phi_i__jp_km = p_cell[i] [jp][km]->phi;
  Real phi_ip_jp_km = p_cell[ip][jp][km]->phi;
  Real phi_im_jm_km = p_cell[im][jm][km]->phi;
  Real phi_i__jm_km = p_cell[i] [jm][km]->phi;
  Real phi_ip_jm_km = p_cell[ip][jm][km]->phi;

  // Calculate the levelset values at cell corners by averaging the levelset values of the adjacent cells
  Real phi_corner_vectors_[8][8] = {{phi_im_jm_km, phi_i__jm_km, phi_im_j__km, phi_i__j__km,   phi_im_jm_k, phi_i__jm_k, phi_im_j__k, phi_i__j__k },
                                    {phi_im_jm_kp, phi_i__jm_kp, phi_im_j__kp, phi_i__j__kp,   phi_im_jm_k, phi_i__jm_k, phi_im_j__k, phi_i__j__k },
                                    {phi_i__jm_kp, phi_ip_jm_kp, phi_i__j__kp, phi_ip_j__kp,   phi_i__jm_k, phi_ip_jm_k, phi_i__j__k, phi_ip_j__k },
                                    {phi_i__jm_km, phi_ip_jm_km, phi_i__j__km, phi_ip_j__km,   phi_i__jm_k, phi_ip_jm_k, phi_i__j__k, phi_ip_j__k },
                                    {phi_im_jp_km, phi_i__jp_km, phi_im_j__km, phi_i__j__km,   phi_im_jp_k, phi_i__jp_k, phi_im_j__k, phi_i__j__k },
                                    {phi_im_jp_kp, phi_i__jp_kp, phi_im_j__kp, phi_i__j__kp,   phi_im_jp_k, phi_i__jp_k, phi_im_j__k, phi_i__j__k },
                                    {phi_i__jp_kp, phi_ip_jp_kp, phi_i__j__kp, phi_ip_j__kp,   phi_i__jp_k, phi_ip_jp_k, phi_i__j__k, phi_ip_j__k },
                                    {phi_i__jp_km, phi_ip_jp_km, phi_i__j__km, phi_ip_j__km,   phi_i__jp_k, phi_ip_jp_k, phi_i__j__k, phi_ip_j__k }};

  std::array<std::vector<Real>,8> phi_corner_vectors;

  for (int ii = 0; ii < 8; ii++){
    for (int jj = 0; jj < 8; jj++){
      phi_corner_vectors[ii].push_back(phi_corner_vectors_[ii][jj]);
    }
  }

//  phi_corner_vectors[0] = {phi_im_jm_km, phi_i__jm_km, phi_im_j__km, phi_i__j__km,   phi_im_jm_k, phi_i__jm_k, phi_im_j__k, phi_i__j__k }; // - - -
//  phi_corner_vectors[1] = {phi_im_jm_kp, phi_i__jm_kp, phi_im_j__kp, phi_i__j__kp,   phi_im_jm_k, phi_i__jm_k, phi_im_j__k, phi_i__j__k }; // - - +
//  phi_corner_vectors[2] = {phi_i__jm_kp, phi_ip_jm_kp, phi_i__j__kp, phi_ip_j__kp,   phi_i__jm_k, phi_ip_jm_k, phi_i__j__k, phi_ip_j__k }; // + - +
//  phi_corner_vectors[3] = {phi_i__jm_km, phi_ip_jm_km, phi_i__j__km, phi_ip_j__km,   phi_i__jm_k, phi_ip_jm_k, phi_i__j__k, phi_ip_j__k }; // + - -
//  phi_corner_vectors[4] = {phi_im_jp_km, phi_i__jp_km, phi_im_j__km, phi_i__j__km,   phi_im_jp_k, phi_i__jp_k, phi_im_j__k, phi_i__j__k }; // - + -
//  phi_corner_vectors[5] = {phi_im_jp_kp, phi_i__jp_kp, phi_im_j__kp, phi_i__j__kp,   phi_im_jp_k, phi_i__jp_k, phi_im_j__k, phi_i__j__k }; // - + +
//  phi_corner_vectors[6] = {phi_i__jp_kp, phi_ip_jp_kp, phi_i__j__kp, phi_ip_j__kp,   phi_i__jp_k, phi_ip_jp_k, phi_i__j__k, phi_ip_j__k }; // + + +
//  phi_corner_vectors[7] = {phi_i__jp_km, phi_ip_jp_km, phi_i__j__km, phi_ip_j__km,   phi_i__jp_k, phi_ip_jp_k, phi_i__j__k, phi_ip_j__k }; // + + -

  for(int u = 0; u < phi_corner_vectors.size(); ++u) {
      std::sort( phi_corner_vectors[u].begin(), phi_corner_vectors[u].end() );
      phi_corner_values[u] = std::accumulate(phi_corner_vectors[u].begin(), phi_corner_vectors[u].end(), 0.0);
      phi_corner_values[u] *= 0.125;
  }
}
//------------------------------------------------------------
// Get_global_scale
//------------------------------------------------------------
void Levelset_package::Get_global_scale(p_Levelset_levelinfo lset_level_info)
{
  maximum_phi= 0.;
  maximum_curv= 0.;
  minimum_curv= 1.e20;
  
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        p_Levelset_cell Mcell = p_cell[i][j][k];
        if(Mcell->phi > 0.)
        {
          maximum_phi = AMAX1(maximum_phi, fabs(Mcell->phi));
        }
        if(Mcell->tag_interface == CUT_CELL)
        {
          maximum_curv = AMAX1(maximum_curv, fabs(Mcell->curv));
          minimum_curv = AMIN1(minimum_curv, fabs(Mcell->curv));
        }
      }
}
//------------------------------------------------------------
// Extend
//------------------------------------------------------------
void Levelset_package::Extend_psi(Levelset *level_set)
{
  my_real dcell = level_set->lset_level_info[level]->dcell;
  Real dl = level_set->lset_level_info[level]->dl;
  Real dtau = pow(0.5, Real(DIM-1))*dl;//0.5*level_set->lset_level_info[level]->dl;
  Real n_d, dv;
  p_Levelset_cell midcell, cell_mm, cell_pp;

  //extending the variables
  //first order upwind scheme
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        midcell = p_cell[i][j][k];

        midcell->LU[0] = 0.0;

        #if (DIM_X)
        cell_mm = p_cell[i-1][j][k];
        cell_pp = p_cell[i+1][j][k];
        n_d = midcell->phi > 0. ? midcell->n_x : -midcell->n_x;
        n_d /= dcell.i;
        //psi
        dv = n_d >= 0.0 ? midcell->psi - cell_mm->psi : cell_pp->psi - midcell->psi;
        midcell->LU[0] += n_d*dv;
        #endif

        #if (DIM_Y)
        cell_mm = p_cell[i][j-1][k];
        cell_pp = p_cell[i][j+1][k];
        n_d = midcell->phi > 0.0 ? midcell->n_y : -midcell->n_y;
        n_d /= dcell.j;
        //psi
        dv = n_d >= 0.0 ? midcell->psi - cell_mm->psi : cell_pp->psi - midcell->psi;
        midcell->LU[0] += n_d*dv;
        #endif

        #if (DIM_Z)
        cell_mm = p_cell[i][j][k-1];
        cell_pp = p_cell[i][j][k+1];
        n_d = midcell->phi > 0.0 ? midcell->n_z : -midcell->n_z;
        n_d /= dcell.k;
        //psi
        dv = n_d >= 0.0 ? midcell->psi - cell_mm->psi : cell_pp->psi - midcell->psi;
        midcell->LU[0] += n_d*dv;
        #endif

        midcell->LU[0] *= dtau;
        midcell->LU[0] = -midcell->LU[0];
      }
}
//------------------------------------------------------------
// Update_Extend
//------------------------------------------------------------
void Levelset_package::Update_extend_psi(Levelset *level_set)
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        if(p_cell[i][j][k]->tag_characteristic != SINGULARITY_CELL && p_cell[i][j][k]->tag_interface != CUT_CELL){
          p_cell[i][j][k]->Update_extend_psi();
        }
      }
}
//------------------------------------------------------------
// Update_Extend
//------------------------------------------------------------
void Levelset_package::Update_extend_psi_reinitialize(Levelset *level_set)
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        p_cell[i][j][k]->Update_extend_psi();
      }
}
//-------------------------------------------------------
// Reinitialize_psi_backup
//-------------------------------------------------------
void Levelset_package::Reinitialize_psi_backup (Levelset *level_set)
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        p_cell[i][j][k]->psi_1 = p_cell[i][j][k]->psi;
      }
}
//-------------------------------------------------------
// Reinitialize_psi_reset
//-------------------------------------------------------
void Levelset_package::Reinitialize_psi_reset (Levelset *level_set)
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        p_cell[i][j][k]->LU[0] = 0.;
      }
}
//-------------------------------------------------------
// Reinitialize_psi_increment
//-------------------------------------------------------
void Levelset_package::Reinitialize_psi_increment (Levelset *level_set)
{
  Real vrbl, midvalue, inter;
  p_Levelset_cell Tcell;

  Real dx2 = level_set->lset_level_info[level]->dcell.i;
  Real dy2 = level_set->lset_level_info[level]->dcell.j;
  Real dz2 = level_set->lset_level_info[level]->dcell.k;
  Real dl2 = level_set->lset_level_info[level]->dl;
  Real dtau = pow(0.5, Real(DIM-1))*dl2;//0.5*level_set->lset_level_info[level]->dl;

  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        Tcell = p_cell[i][j][k];
        midvalue = Tcell->psi;
        vrbl = Tcell->psi;
      
        if(Tcell->tag_interface == CUT_CELL && Tcell->tag_characteristic != SINGULARITY_CELL){    
          Real s = vrbl;
          Real inter2 = 0.0;

          //x direction
          Real dv_p, dv_m;
          #if (DIM_X==1)
          dv_p = p_cell[i+1][j][k]->psi - midvalue;
          dv_m = midvalue - p_cell[i-1][j][k]->psi;
          Real dv_x = dv_p;
          if(s*dv_p >= 0.0 && s*dv_m >= 0.0) dv_x = dv_m;
          if(s*dv_p <= 0.0 && s*dv_m <= 0.0) dv_x = dv_p;
          if(s*dv_p > 0.0 && s*dv_m < 0.0) dv_x = 0.0;
          if(s*dv_p < 0.0 && s*dv_m > 0.0) {
            if(fabs(dv_p) < fabs(dv_m)) dv_x = dv_m;
          }
          inter2 += dv_x*dv_x/dx2;
          #endif
      
          //y direction
          #if (DIM_Y==1)
          dv_p = p_cell[i][j+1][k]->psi - midvalue;
          dv_m = midvalue - p_cell[i][j-1][k]->psi;
          Real dv_y = dv_p;
          if(s*dv_p >= 0.0 && s*dv_m >= 0.0) dv_y = dv_m;
          if(s*dv_p <= 0.0 && s*dv_m <= 0.0) dv_y = dv_p;
          if(s*dv_p > 0.0 && s*dv_m < 0.0) dv_y = 0.0;
          if(s*dv_p < 0.0 && s*dv_m > 0.0) {
            if(fabs(dv_p) < fabs(dv_m)) dv_y = dv_m;
          }
          inter2 += dv_y*dv_y/dy2;
          #endif

          //z direction
          #if (DIM_Z==1)
          dv_p = p_cell[i][j][k+1]->psi - midvalue;
          dv_m = midvalue - p_cell[i][j][k-1]->psi;
          Real dv_z = dv_p;
          if(s*dv_p >= 0.0 && s*dv_m >= 0.0) dv_z = dv_m;
          if(s*dv_p <= 0.0 && s*dv_m <= 0.0) dv_z = dv_p;
          if(s*dv_p > 0.0 && s*dv_m < 0.0) dv_z = 0.0;
          if(s*dv_p < 0.0 && s*dv_m > 0.0) {
            if(fabs(dv_p) < fabs(dv_m)) dv_z = dv_m;
          }
          inter2 += dv_z*dv_z/dz2;
          #endif

          //time step
          Real inter1 = sqrt(inter2);
          inter = -dtau*(inter1 - 1.0);
        
          Tcell->LU[0]  = inter;
        }
      }
}
//-------------------------------------------------------
// Check_covergence
//-------------------------------------------------------
bool Levelset_package::Check_covergence_psi (Levelset *level_set)
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        if(fabs(p_cell[i][j][k]->psi - p_cell[i][j][k]->psi_1) > level_set->convergence_error){
          return true;
        }
      }
  return false;
}
//-------------------------------------------------------
// Get_psi_max
//-------------------------------------------------------
void Levelset_package::Get_psi_max (Levelset *level_set)
{
  maximum_psi = 0.;
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        if(p_cell[i][j][k]->tag_interface == CUT_CELL || p_cell[i][j][k]->tag_characteristic == SINGULARITY_CELL)
          maximum_psi = AMAX1(maximum_psi, fabs(p_cell[i][j][k]->psi));
      }
}
//-------------------------------------------------------
// Redistribute_curv
//-------------------------------------------------------
void Levelset_package::Redistribute_curv (Levelset *level_set)
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        if(p_cell[i][j][k]->tag_interface == CUT_CELL || p_cell[i][j][k]->tag_characteristic == SINGULARITY_CELL){
          p_cell[i][j][k]->curv = AMAX1(p_cell[i][j][k]->curv, (1. - p_cell[i][j][k]->psi/(level_set->maximum_psi + 1.e-20))*level_set->maximum_curv);

//           if (level_set->num_singularity == 0 ) p_cell[i][j][k]->curv = level_set->maximum_curv;
        }
      }
}
//-------------------------------------------------------
// Extend_curv_backup
//-------------------------------------------------------
void Levelset_package::Extend_curv_backup(p_Levelset_levelinfo lset_level_info)
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        p_cell[i][j][k]->curv_1 = p_cell[i][j][k]->curv;
      }
}
//-------------------------------------------------------
// Extend_curv
//-------------------------------------------------------
void Levelset_package::Extend_curv(p_Levelset_levelinfo lset_level_info)
{
  Real dl = lset_level_info->dl;
  Real dtau = pow(0.5, Real(DIM-1))*pow(dl, 2.0);
//   Real dtau = 0.5*pow(dl, 2.0);
  Real n_d, dv, dv_1, dv_2, dx, dy, dz;
  p_Levelset_cell midcell, cell_mm, cell_pp;
  my_real dcell = lset_level_info->dcell;
  dx = dcell.i;
  dy = dcell.j;
  dz = dcell.k;
  
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++)
      {
        midcell = p_cell[i][j][k];

        midcell->LU[0] = 0.0;

        #if (DIM_X)
        cell_mm = p_cell[i-1][j][k];
        cell_pp = p_cell[i+1][j][k];
        n_d = midcell->n_x*midcell->n_x;
        n_d /= dx*dx;

        //curvature
        dv = cell_pp->curv - 2.* midcell->curv + cell_mm->curv;
        midcell->LU[0] += n_d*dv;
          
          #if (DIM_Y)
          n_d = midcell->n_x*midcell->n_y;
          n_d /= dx*dy;

          //curvature
          dv_1 = (p_cell[i-1][j+1][k]->curv - p_cell[i-1][j-1][k]->curv)/2.;
          dv_2 = (p_cell[i+1][j+1][k]->curv - p_cell[i+1][j-1][k]->curv)/2.;

          dv = (dv_2 - dv_1)/2.;
          midcell->LU[0] += n_d*dv;

          #endif  

          #if (DIM_Z)
          n_d = midcell->n_x*midcell->n_z;
          n_d /= dx*dz;

          //curvature
          dv_1 = (p_cell[i-1][j][k+1]->curv - p_cell[i-1][j][k-1]->curv)/2.;
          dv_2 = (p_cell[i+1][j][k+1]->curv - p_cell[i+1][j][k-1]->curv)/2.;

          dv = (dv_2 - dv_1)/2.;
          midcell->LU[0] += n_d*dv;               
          #endif

        #endif

        #if (DIM_Y)
        cell_mm = p_cell[i][j-1][k];
        cell_pp = p_cell[i][j+1][k];
        n_d = midcell->n_y * midcell->n_y;
        n_d /= dy*dy;

        //curvature
        dv = cell_pp->curv - 2.* midcell->curv + cell_mm->curv;
        midcell->LU[0] += n_d*dv;

          #if (DIM_X)
          n_d = midcell->n_x*midcell->n_y;
          n_d /= dx*dy;

          //curvature
          dv_1 = (p_cell[i+1][j-1][k]->curv - p_cell[i-1][j-1][k]->curv)/2.;
          dv_2 = (p_cell[i+1][j+1][k]->curv - p_cell[i-1][j+1][k]->curv)/2.;

          dv = (dv_2 - dv_1)/2.;
          midcell->LU[0] += n_d*dv;

          #endif  

          #if (DIM_Z)
          n_d = midcell->n_y*midcell->n_z;
          n_d /= dy*dz;

          //curvature
          dv_1 = (p_cell[i][j-1][k+1]->curv - p_cell[i][j-1][k-1]->curv)/2.;
          dv_2 = (p_cell[i][j+1][k+1]->curv - p_cell[i][j+1][k-1]->curv)/2.;

          dv = (dv_2 - dv_1)/2.;
          midcell->LU[0] += n_d*dv;             
          #endif

        #endif

        #if (DIM_Z)
        cell_mm = p_cell[i][j][k-1];
        cell_pp = p_cell[i][j][k+1];
        n_d = midcell->n_z * midcell->n_z;
        n_d /= dz*dz;

        //curvature
        dv = cell_pp->curv - 2.* midcell->curv + cell_mm->curv;
        midcell->LU[0] += n_d*dv;

          #if (DIM_X)
          n_d = midcell->n_x*midcell->n_z;
          n_d /= dx*dz;

          //curvature
          dv_1 = (p_cell[i+1][j][k-1]->curv - p_cell[i-1][j][k-1]->curv)/2.;
          dv_2 = (p_cell[i+1][j][k+1]->curv - p_cell[i-1][j][k+1]->curv)/2.;

          dv = (dv_2 - dv_1)/2.;
          midcell->LU[0] += n_d*dv;

          #endif  

          #if (DIM_Y)
          n_d = midcell->n_y*midcell->n_z;
          n_d /= dy*dz;

          //curvature
          dv_1 = (p_cell[i][j+1][k-1]->curv - p_cell[i][j-1][k-1]->curv)/2.;
          dv_2 = (p_cell[i][j+1][k+1]->curv - p_cell[i][j-1][k+1]->curv)/2.;

          dv = (dv_2 - dv_1)/2.;
          midcell->LU[0] += n_d*dv;
          #endif

        #endif

        midcell->LU[0] *= dtau;
      }
}
//-------------------------------------------------------
// Update_extend_curv
//-------------------------------------------------------
void Levelset_package::Update_extend_curv(Levelset *level_set)
{
  #ifdef _READ_SDF_
  my_real   dpkg = level_set->lset_level_info[0]->dpkg;
  my_real  dcell = level_set->lset_level_info[0]->dcell;
  my_real  box_l = level_set->box_l;
  int n_artifact = level_set->num_artifact_region;
  Real     *radi = level_set->radi_artifact_region;
  my_real  *cent = level_set->center_artifact_region;
  #endif

  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        if(p_cell[i][j][k]->tag_interface != CUT_CELL){
//         if(p_cell[i][j][k]->curv <= 0.2){
          p_cell[i][j][k]->Update_extend_curv();
        #ifdef _READ_SDF_
        }else{
          my_real coord = get_cell_position (index.i, index.j, index.k, i, j, k, dpkg, dcell, box_l);

          for (int ip = 0; ip<n_artifact; ip++){
            Real dist = get_distance_2p (coord, cent[ip]);
            if (dist < radi[ip]){
              p_cell[i][j][k]->Update_extend_curv();
            }
          }
        #endif
        }
      }
}
//------------------------------------------------------------
// Extend_curvature_for_smooth
//------------------------------------------------------------
void Levelset_package::Extend_curvature_for_smooth(Levelset *level_set)
{
  my_real dcell = level_set->lset_level_info[level]->dcell;
  Real dl = level_set->lset_level_info[level]->dl;
  Real dtau = pow(0.5, Real(DIM-1))*dl;//0.5*level_set->lset_level_info[level]->dl;
  Real n_d, dv;
  p_Levelset_cell midcell, cell_mm, cell_pp;

  //extending the variables
  //first order upwind scheme
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        midcell = p_cell[i][j][k];

        midcell->LU[0] = 0.0;

        #if (DIM_X)
        cell_mm = p_cell[i-1][j][k];
        cell_pp = p_cell[i+1][j][k];
        n_d = midcell->phi > 0. ? midcell->n_x : -midcell->n_x;
        n_d /= dcell.i;
        //curv
        dv = n_d >= 0.0 ? midcell->curv - cell_mm->curv : cell_pp->curv - midcell->curv;
        midcell->LU[0] += n_d*dv;
        #endif

        #if (DIM_Y)
        cell_mm = p_cell[i][j-1][k];
        cell_pp = p_cell[i][j+1][k];
        n_d = midcell->phi > 0.0 ? midcell->n_y : -midcell->n_y;
        n_d /= dcell.j;
        //curv
        dv = n_d >= 0.0 ? midcell->curv - cell_mm->curv : cell_pp->curv - midcell->curv;
        midcell->LU[0] += n_d*dv;
        #endif

        #if (DIM_Z)
        cell_mm = p_cell[i][j][k-1];
        cell_pp = p_cell[i][j][k+1];
        n_d = midcell->phi > 0.0 ? midcell->n_z : -midcell->n_z;
        n_d /= dcell.k;
        //curv
        dv = n_d >= 0.0 ? midcell->curv - cell_mm->curv : cell_pp->curv - midcell->curv;
        midcell->LU[0] += n_d*dv;
        #endif

        midcell->LU[0] *= dtau;
        midcell->LU[0] = -midcell->LU[0];
      }
}
//------------------------------------------------------------
// Update_curvature_for_smooth
//------------------------------------------------------------
void Levelset_package::Update_curvature_for_smooth(Levelset *level_set)
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        if(p_cell[i][j][k]->tag_interface != CUT_CELL){
          p_cell[i][j][k]->curv += p_cell[i][j][k]->LU[0];
        }
      }
}
//-------------------------------------------------------
// Update_extend_curv
//-------------------------------------------------------
void Levelset_package::Update_extend_curv_for_smoothing(p_Levelset_levelinfo lset_level_info)
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        if(p_cell[i][j][k]->tag_interface == CUT_CELL){
          p_cell[i][j][k]->Update_extend_curv();
        }
      }
}
//-------------------------------------------------------
// Check_covergence_curv
//-------------------------------------------------------
bool Levelset_package::Check_covergence_curv(Levelset *level_set)
{
  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        if(fabs(p_cell[i][j][k]->curv - p_cell[i][j][k]->curv_1) > level_set->convergence_error){
          return true;
        }
      }
  return false;
}
//-------------------------------------------------------
// Copy_curvature
//-------------------------------------------------------
void Levelset_package::Calculate_volume_mass(Levelset *level_set)
{
  Real dl = level_set->lset_level_info[level]->dl;
  Real dx, dy, dz;
  
  my_real  dpkg = level_set->lset_level_info[level]->dpkg;
  my_real dcell = level_set->lset_level_info[level]->dcell;
  dx = dcell.i;
  dy = dcell.j;
  dz = dcell.k;

  total_mass         = 0.;
  total_mass_surface = 0.;
  total_mass_segment = 0.;
  total_volume       = 0.;
  total_area         = 0.;
  total_length       = 0.;
  
  Real corner[1+2*DIM_X][1+2*DIM_Y][1+2*DIM_Z];
  Real phii0[1+2*DIM_X][1+2*DIM_Y][1+2*DIM_Z];

  int xw = 1+2*DIM_X;
  int yw = 1+2*DIM_Y;
  int zw = 1+2*DIM_Z;

  for (int i = 0; i < ICPX; i++)
    for (int j = 0; j < ICPY; j++)
      for (int k = 0; k < ICPZ; k++){
        p_Levelset_cell Mcell = p_cell[i][j][k];

        // // TODO: more accurate volume calculation
        // if(Mcell->phi > 0. || fabs(Mcell->phi) < dl){
        //   Real volume  = (DIM_X==1 ? dx : 1.)*(DIM_Y==1 ? dy : 1.)*(DIM_Z==1 ? dz : 1.);

        //   Mcell->scale = 0.;
          
        //   Mcell->scale = level_set->Get_scale(Mcell->curv, Mcell->phi);

        //   Real density;

        //   if(DIM_X*DIM_Y*DIM_Z)
        //   density = 1./powern(Mcell->scale*2.0,3);
        //   else
        //   density = 1./powern(Mcell->scale*2.0,2);
        
        //   total_mass +=  density*volume;
        //   total_volume += volume;

        //   if(Mcell->tag_interface == 1)
        //     // TODO: more accurate surface calculation
        //     total_mass_surface +=  pow(dl/Mcell->scale/2.0, DIM-1);
        // }

        // JZ20181122 implementing more accurate volume calculation
        if(Mcell->phi > 0. || fabs(Mcell->phi) < dl){

          Real volume  = 0.;
          Real area    = 0.;
          Real length  = 0.;

          GeometryCalculator gc;

          Mcell->scale = level_set->Get_scale(Mcell->curv, Mcell->phi);

          if(Mcell->tag_interface == CUT_CELL){
            int ii = i-DIM_X;
            int jj = j-DIM_Y;
            int kk = k-DIM_Z;
    
            for(int r=0; r<xw; r++){
              for(int s=0; s<yw; s++){
                for(int t=0; t<zw; t++){
                  phii0 [r][s][t] = p_cell[ii+r][jj+s][kk+t]->phi;
                  corner[r][s][t] = 0.;
                }
              }
            }

            gc.Get_corner(phii0, corner);

            volume = gc.Get_cutcell_volume(corner, dx, dy, dz, dl);

            area = gc.Get_surface_area(corner, dx, dy, dz, dl);
            
            total_area += area;

            total_mass_surface +=  area/powern(Mcell->scale,DIM-1);

          }else{
            volume = (DIM_X==1 ? dx : 1.)*(DIM_Y==1 ? dy : 1.)*(DIM_Z==1 ? dz : 1.);
            area = 0.;
          }

          #if DIM == 3
          if(Mcell->tag_characteristic == SEGMENT_CELL){

            my_real min = get_cell_position_corner (index.i, index.j, index.k, i, j, k, dpkg, dcell, level_set->box_l);
            my_real max = my_add_data (min, dcell);

            length = gc.Get_segment_length (min, max, level_set->segment[Mcell->idx_characteristic].first, level_set->segment[Mcell->idx_characteristic].second);

            total_length += length;
            total_mass_segment += length/Mcell->scale;
          }
          #endif
        
          total_mass   += volume/powern(Mcell->scale,DIM);
          total_volume += volume;

          Mcell->vol  = volume;
          Mcell->area = area;
        }
      }
}
//-------------------------------------------------------
// Copy_curvature
//-------------------------------------------------------
void Levelset_package::Copy_curvature(Levelset *level_set)
{

}
//------------------------------------------------------------
// Output_level_set
//------------------------------------------------------------
void Levelset_package::Output_level_set(Levelset *level_set, communicator &world, char *filename, int n)
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
  LEN = 7;
  out<<"phi, n_x, n_y, n_z, curv, psi, interface"<<"\n";
  
  my_int cell_id;
  cell_id.i = index.i*ICPX;
  cell_id.j = index.j*ICPY;
  cell_id.k = index.k*ICPZ;
  
  my_real box_l = level_set->box_l;
  my_real dcell = level_set->lset_level_info[level]->dcell;
  
  out<<"zone t='sub_filed_"<<cell_id.i<<"_"<<cell_id.j<<"_"<<cell_id.k<<"'  i= "<<ICPX+DIM_X<<"  j= "<<ICPY+DIM_Y<<"  k= "<<ICPZ+DIM_Z<<"  DATAPACKING=BLOCK, VARLOCATION=([";
  int pos_s = DIM_X+DIM_Y+DIM_Z+1;
  out<<pos_s<<"-";
  out<<2*pos_s -1 + LEN-1<<"]=CELLCENTERED) SOLUTIONTIME="<<n<<"\n";

#if DIM_X
  for(int k=0; k<=ICPZ+DIM_Z-1; k++){
    for(int j=0; j<=ICPY+DIM_Y-1; j++){
      for(int i=0; i<=ICPX+DIM_X-1; i++){
        out<<dcell.i*(i+cell_id.i)+box_l.i<<" ";
      }
      out<<"\n";
    }
  }
#endif
#if DIM_Y
  for(int k=0; k<=ICPZ+DIM_Z-1; k++){
    for(int j=0; j<=ICPY+DIM_Y-1; j++){
      for(int i=0; i<=ICPX+DIM_X-1; i++){
        out<<dcell.j*(j+cell_id.j)+box_l.j<<" ";
      }
      out<<"\n";
    }
  }
#endif
#if DIM_Z
  for(int k=0; k<=ICPZ+DIM_Z-1; k++){
    for(int j=0; j<=ICPY+DIM_Y-1; j++){
      for(int i=0; i<=ICPX+DIM_X-1; i++){
        out<<dcell.k*(k+cell_id.k)+box_l.k<<" ";
      }
      out<<"\n";
    }
  }
#endif

  // phi
  for(int k=0; k<ICPZ; k++){
    for(int j=0; j<ICPY; j++){
      for(int i=0; i<ICPX; i++){
        out<<lset_cell[i][j][k].phi<<" ";
      }
      out<<"\n";
    }
  }

  // n_x
  for(int k=0; k<ICPZ; k++){
    for(int j=0; j<ICPY; j++){
      for(int i=0; i<ICPX; i++){
        out<<lset_cell[i][j][k].n_x<<" ";
      }
      out<<"\n";
    }
  }

  // n_y
  for(int k=0; k<ICPZ; k++){
    for(int j=0; j<ICPY; j++){
      for(int i=0; i<ICPX; i++){
        out<<lset_cell[i][j][k].n_y<<" ";
      }
      out<<"\n";
    }
  }

  // n_z
  for(int k=0; k<ICPZ; k++){
    for(int j=0; j<ICPY; j++){
      for(int i=0; i<ICPX; i++){
        out<<lset_cell[i][j][k].n_z<<" ";
      }
      out<<"\n";
    }
  }

  // curv
  for(int k=0; k<ICPZ; k++){
    for(int j=0; j<ICPY; j++){
      for(int i=0; i<ICPX; i++){
        out<<lset_cell[i][j][k].curv<<" ";
      }
      out<<"\n";
    }
  }

  // psi
  for(int k=0; k<ICPZ; k++){
    for(int j=0; j<ICPY; j++){
      for(int i=0; i<ICPX; i++){
        out<<lset_cell[i][j][k].psi<<" ";
      }
      out<<"\n";
    }
  }
  
    // interface
  for(int k=0; k<ICPZ; k++){
    for(int j=0; j<ICPY; j++){
      for(int i=0; i<ICPX; i++){
        out<<int(lset_cell[i][j][k].tag_interface)<<" ";
      }
      out<<"\n";
    }
  }
  out.close();
}
