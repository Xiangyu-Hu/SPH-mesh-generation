#include "glbfunc.h"
#include "level_set.h"
#include "sph_mesh_generation.h"

/***************************************************/
/*                                                 */
/*     Functions defined in class "Levelset"       */
/*                                                 */
/***************************************************/

//-----------------------------------------------------
// constructor
//-----------------------------------------------------
Levelset::Levelset()
{
  num_level = 0;
  glbl_total_num_particle = 0;
  glbl_total_num_pkg = 0;
  glbl_num_interface_cell = 0;
  Lmin = 0;
  Lmax = 0;
  n_out = 0;
  
  total_mass = 0.;
  total_mass_surface = 0.;
  total_mass_segment = 0.;
  total_volume = 0.;
  total_area = 0.;
  total_length = 0.;
  maximum_phi = 0.;
  maximum_psi = 0.;
  maximum_dl = 0.;
  middel_dl = 0.;
  minimum_dl = 0.;
  maximum_curv = 0.;
  minimum_curv = 1.e20;
  
  #ifdef _MPI_
  error_tolerance = 0.5;
    #if defined (_CVP_LSET_INIT_)
    num_color_tmp   = 1;
    #endif
  #endif
  lset_pkg.clear();
  interface_pkg.clear();
  lset_bc_pkg.clear();
}
//-----------------------------------------------------
// initialize
//-----------------------------------------------------
void Levelset::Initialize(SOLVER *sph, communicator &world)
{
  num_level = 1;
  glbl_total_num_particle = sph->glbl_total_num_particle;
  glbl_total_num_pkg = 0;
  
  ini_num_cell = sph->ini_num_cell;
  domain = sph->domain;
  box_l = sph->box_l;
  box_r = sph->box_r;
  Lmin = sph->Lmin;
  Lmax = sph->Lmax;

  #ifdef _READ_SDF_
  size_t _size = sizeof (sph->sdf_file);
  snprintf(sdf_file, _size, "%s", sph->sdf_file);
  cells_to_cut = sph->cells_to_cut;
  #endif
   
  lset_level_info = new p_Levelset_levelinfo[num_level];
  for (int i = 0; i < num_level; i++){
    lset_level_info[i] = new Levelset_levelinfo;
    lset_level_info[i]->Initialize(i, this, world);
  }
  for (int i = 0; i < num_level; i++){
    lset_level_info[i]->Set_cell_topology(this, world);
  }
  
  glbl_total_num_pkg = lset_pkg.size();
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Class Levelset is initialized\n";
  out.close();
#endif
}
//-----------------------------------------------------
// define level set function
//-----------------------------------------------------
Real Levelset::F_phi(my_real pos)
{
  Real phi = 0.;
  Real R;
  my_real center;
  
  center.i = DIM_X == 1 ? 50. : 0.;
  center.j = DIM_Y == 1 ? 50. : 0.;
  center.k = DIM_Z == 1 ? 50. : 0.;
  R = 43.;
  
  Real dist = R - get_distance_2p (pos, center);
  phi = dist;
  return phi;

  // my_real center;
  // center.i = 50.;
  // center.j = 50.;
  // center.k = 0.;

  // Real phi;
  // if(pos.i >= center.i)
  //   pos.i = 100. - pos.i;

  // Real R = 45.0;
  // Real W = 7.5;
  // Real H = 30.63;
  // Real low_Y  = center.j - sqrt(R*R-W*W);
  // Real left_X = center.i - W;
  // Real left_Y = center.j + H;
  
  // Real dist =get_distance (my_minus_data (pos, center));

  // if(pos.i <= center.i)
  // {
  //   if(pos.j >= left_Y)
  //   {
  //     if(dist >= R)
  //       phi = dist - R;
  //     else
  //     {
  //       if(pos.i <= left_X)
  //         phi = -AMIN1(R - dist,  sqrt((pos.i-left_X)*(pos.i-left_X) + (pos.j-left_Y)*(pos.j-left_Y))); 
  //       else
  //         phi = -AMIN1(R - dist,  fabs(pos.j - left_Y)); 
  //     }
  //   }else
  //   {
  //     if(pos.i <= left_X)
  //     {
  //       if(dist >= R)
  //         phi = dist - R;
  //       else
  //       {
  //         phi = -AMIN1(R - dist,  fabs(pos.i - left_X)); 
  //       }
  //     }else
  //     {
  //       if(pos.j >= low_Y)
  //         phi = AMIN1(fabs(pos.j - left_Y),  fabs(pos.i - left_X)); 
  //       else
  //         phi = sqrt((pos.i-left_X)*(pos.i-left_X) + (pos.j-low_Y)*(pos.j-low_Y));
  //     }

  //   }
  // }

  // phi = -phi;
  // return phi;

  // Real phi = 0.;
  // Real R = 25.;
  // Real r = 20.;
  // Real x = pos.i - 50.;
  // Real y = pos.j - 50.;
  // Real z = pos.k - 50.;
  
  // phi = z*z + pow((sqrt(x*x + y*y) - R), 2.) - r*r;
  
  // return -1.*sgn1(phi)*sqrt(fabs(phi));

}
//-----------------------------------------------------
// Calculate_density_field
//-----------------------------------------------------
void Levelset::Calculate_target_density_field(SOLVER *sph, communicator &world)
{
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Level_set is processing\n";
    n_out = 0;

    Copy_phi_value (world);
    
    Calculate_level_set(world);
    
    Find_interface_packages(world);
    
    // Output_level_set_vti(n_out++, world);

    // // Smoothing_curvarure_field (world);

    // Calculate_psi(world);
    
    // Output_level_set_vti(n_out++, world);

    // Calculate_global_effective_curv(world);
    
    Calculate_total_volume_mass(world);
    
    Output_level_set_vti(n_out++, world);
  
    cout<<"<<<<< Target density field is calculated\n";
    cout<<"**********************************************************\n";
  }

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Calculate_density_field finished\n";
  out.close();
#endif
}
//-----------------------------------------------------
// Copy_phi_value
//-----------------------------------------------------
void Levelset::Copy_phi_value(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      lset_pkg[i]->Copy_phi_value();
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Copy_phi_value finished\n";
  out.close();
#endif
}
//-----------------------------------------------------
// Calculate level set
//-----------------------------------------------------
void Levelset::Calculate_level_set(communicator &world)
{
  // TODO
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int i_level = lset_pkg[i]->level;
      lset_pkg[i]->Boundary(lset_level_info[i_level]);
    }
  }, ap);

  parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int i_level = lset_pkg[i]->level;
      lset_pkg[i]->Get_normal(lset_level_info[i_level]);
      lset_pkg[i]->Get_interface_tag(lset_level_info[i_level]);
      lset_pkg[i]->Get_curvature(lset_level_info[i_level]);
      lset_pkg[i]->Get_global_scale(lset_level_info[i_level]);
    }
  }, ap);

  parallel_for( blocked_range<int>(0, int(lset_bc_pkg.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int i_level = lset_bc_pkg[i]->level;
      lset_bc_pkg[i]->Get_BC_normal(lset_level_info[i_level]);
    }
  }, ap);

  Parallel_get_global_scale(lset_pkg.size(), lset_pkg, this);

  parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int i_level = lset_pkg[i]->level;
      lset_pkg[i]->Get_extended_interface_tag(lset_level_info[i_level]);
    }
  }, ap);

  // parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
  //   for(int i=r.begin(); i!=r.end(); ++i){
  //     int i_level = lset_pkg[i]->level;
  //     lset_pkg[i]->Clean_curvature();
  //   }
  // }, ap);
  
  if (num_singularity > 0){
    my_real Pos_singularity[num_singularity];

    // parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
    //   for(int i=r.begin(); i!=r.end(); ++i){
      for(int i=0; i!=glbl_total_num_pkg; ++i){
        lset_pkg[i]->Get_singularity_tag(this, Pos_singularity);
      }
    // }, ap);
  }

  if (num_segment > 0){
    parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
      for(int num=r.begin(); num!=r.end(); ++num)
        lset_pkg[num]->Get_segment_tag(this);
    }, ap);
  }

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Calculate_level_set finished\n";
  out.close();
#endif
}
//-----------------------------------------------------
// Get_scale
//-----------------------------------------------------
Real Levelset::Get_scale(Real curv, Real phi)
{
  Real scale = 0.;
  // if (curv > maximum_curv_artificial) curv = maximum_curv_artificial;
  
  if (phi > 0. || fabs(phi) < lset_level_info[0]->dl){

    /*****************************************************/
    // scale = middel_dl + pow(fabs(phi)/maximum_phi, 2.0)*(maximum_dl - middel_dl);

    // scale = AMAX1(middel_dl, AMIN1(scale, maximum_dl));

    // scale = minimum_dl + pow((fabs(curv) - maximum_curv)/(fabs(minimum_curv - maximum_curv) + 1.e-20), 2.0)*(scale - minimum_dl);

    /*****************************************************/
    // scale = minimum_dl + fabs(phi)/maximum_phi*(maximum_dl - minimum_dl);

    // scale = minimum_dl + pow(fabs(phi)/maximum_phi, 2.0)*(maximum_dl - minimum_dl);
    /*****************************************************/

    // scale = 1./(1./maximum_dl + 1./23.*(1./minimum_dl - 1./maximum_dl)*fabs(46.-phi-23.));

    /*****************************************************/
    Real factor = 2.5;
    scale = minimum_dl + tanh(factor*fabs(phi)/maximum_phi)/tanh(factor)*(maximum_dl - minimum_dl);
    /*****************************************************/
    // Real factor = 1.0;

    // scale = minimum_dl + fabs(fabs(curv) - maximum_curv)/(fabs(minimum_curv - maximum_curv) + 1.e-20)*(middel_dl - minimum_dl);
    // // scale = middel_dl + tanh(factor*fabs(phi)/maximum_phi)/tanh(factor)*(maximum_dl - middel_dl);

    // // scale = AMAX1(middel_dl, AMIN1(scale, maximum_dl));

    // scale = scale + pow(fabs(phi)/maximum_phi, 2.0)*(maximum_dl - scale);

  }else {
    scale = minimum_dl;
  }

  return scale;
}
//-----------------------------------------------------
// Get_positive_pkg
//-----------------------------------------------------
int Levelset::Get_positive_pkg(std::vector<p_Levelset_package> &lset_positive_pkg)
{
  for(int i = 0; i < glbl_total_num_pkg; i++){
    p_Levelset_package pkg_i = lset_pkg[i];
    // my_real coord = pkg_i->get_pkg_position(pkg_i->index.i, pkg_i->index.j, pkg_i->index.k, lset_level_info[0]->dpkg, box_l);
    // Real phi = Get_phi_at_position (coord);

    // if (phi>=0) lset_positive_pkg.push_back(pkg_i);

    if (pkg_i->total_mass > 0.) lset_positive_pkg.push_back(pkg_i);
  }

  return lset_positive_pkg.size();
}
//-----------------------------------------------------
// Find_interface_packages
//-----------------------------------------------------
Real Levelset::Get_phi_at_position(my_real coord)
{
  p_Levelset_levelinfo current_level = lset_level_info[0];
  
  my_real coord_shift = my_minus_data (coord, box_l);
  
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
  Real coeffecient              [1+DIM_X][1+DIM_Y][1+DIM_Z];
  
  for(int r=0; r<=DIM_X; r++)
    for(int s=0; s<=DIM_Y; s++)
      for(int t=0; t<=DIM_Z; t++){
        Target_smooth_infor_scale[r][s][t] = current_pkg->p_cell[r+im][s+jm][t+km]->phi;
        coeffecient[r][s][t] = (DIM_X==1 ? fabs(ii-(1.0-r)):1.0)*(DIM_Y==1 ? fabs(jj-(1.0-s)):1.0)*(DIM_Z==1 ? fabs(kk-(1.0-t)):1.0);
      }

  Real phi = 0. ;

  for(int r=0; r<=DIM_X; r++)
    for(int s=0; s<=DIM_Y; s++)
      for(int t=0; t<=DIM_Z; t++){
        phi    += Target_smooth_infor_scale[r][s][t]*coeffecient[r][s][t];
      }

  return phi;
}
//-----------------------------------------------------
// Find_interface_packages
//-----------------------------------------------------
void Levelset::Find_interface_packages(communicator &world)
{
  for (int i = 0; i < num_level; i++){
    lset_level_info[i]->Find_interface_packages(this, world);
  }
  
  int num_interface_cell = 0;
  glbl_num_interface_cell = 0;
  
  Parallel_get_num_of_interface_cell(interface_pkg.size(), interface_pkg, this, num_interface_cell);
  
  // JZ20181220 :: the levelset solver is not for MPI yet
  glbl_num_interface_cell = num_interface_cell;

  // reduce(world, num_interface_cell, glbl_num_interface_cell, std::plus<int>(), 0);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Find_interface_packages finished\n";
  out.close();
#endif
}
//-----------------------------------------------------
// Get_extended_cell_tags
//-----------------------------------------------------
void Levelset::Get_extended_cell_tags(communicator &world)
{
  for (int i = 0; i < num_level; i++){
    lset_level_info[i]->Get_extended_cell_tags(this, world);
  }

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Get_extended_cell_tags finished\n";
  out.close();
#endif
}
//-----------------------------------------------------
// Calculate_psi
//-----------------------------------------------------
void Levelset::Calculate_psi(communicator &world)
{
  int iterate = 0;
  int bnd = 50000;
  if (num_singularity == 0) goto here;
  
  // maximum_curv = maximum_curv_artificial;

  if (world.rank() == 0)
    cout<<"<<<<< The Psi field is started calculating:"<<endl;

  for(int n=0; n<bnd; n++){

    int tag_convergence = 0;

    Extend_psi (world);
    
    Reinitialize_psi (world, tag_convergence);
    
    iterate++;

    int n_out_ = 50;
    if (iterate%n_out_ ==0 ){
      if (world.rank() == 0)
        cout<<"<<< Current iteration number : "<<iterate<<endl;
      n_out += n_out_;
//       Output_level_set_dat(iterate, world);
      Output_level_set_vti(n_out, world);
    }
    
    if(tag_convergence == 0) break;
    
  }
  
  Get_psi_max (world);
  
  Redistribute_curv (world);
  
  Recalculate_max_curv (world);
  
  here:
    
  if (world.rank() == 0){
    cout<<"<<<<< The psi field is calculated with iteration : "<<iterate<<" ("<<bnd<<")"<<endl;
    cout<<"**********************************************************\n";
  }

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Calculate_psi finished\n";
  out.close();
#endif
}
//-----------------------------------------------------
// Extend_level_set
//-----------------------------------------------------
void Levelset::Extend_psi(communicator &world)
{
  static affinity_partitioner ap;

  //Extend state values among different nodes
  for(int n=0; n<20; n++){

    int NumofOverlapPkgs = interface_pkg.size();
  
    parallel_for( blocked_range<int>(0, NumofOverlapPkgs),[&](const blocked_range<int>& r){
      for(int num=r.begin(); num!=r.end(); ++num)
        interface_pkg[num]->Extend_psi(this);
    }, ap);

    parallel_for( blocked_range<int>(0, NumofOverlapPkgs),[&](const blocked_range<int>& r){
      for(int num=r.begin(); num!=r.end(); ++num)
        interface_pkg[num]->Update_extend_psi(this);
    }, ap);
  }
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Extend_level_set finished\n";
  out.close();
#endif
}
//-----------------------------------------------------
// Reinitialize_psi
//-----------------------------------------------------
void Levelset::Reinitialize_psi(communicator &world, int &tag_convergence)
{
  static affinity_partitioner ap;
  
  int NumofOverlapPkgs = interface_pkg.size();
  
  parallel_for( blocked_range<int>(0, NumofOverlapPkgs),[&](const blocked_range<int>& r){
    for(int num=r.begin(); num!=r.end(); ++num)
      interface_pkg[num]->Reinitialize_psi_backup(this);
  }, ap);

  parallel_for( blocked_range<int>(0, NumofOverlapPkgs),[&](const blocked_range<int>& r){
    for(int num=r.begin(); num!=r.end(); ++num)
      interface_pkg[num]->Reinitialize_psi_reset(this);
  }, ap);
  
  parallel_for( blocked_range<int>(0, NumofOverlapPkgs),[&](const blocked_range<int>& r){
    for(int num=r.begin(); num!=r.end(); ++num)
      interface_pkg[num]->Reinitialize_psi_increment(this);
  }, ap);
  
  parallel_for( blocked_range<int>(0, NumofOverlapPkgs),[&](const blocked_range<int>& r){
    for(int num=r.begin(); num!=r.end(); ++num)
      interface_pkg[num]->Update_extend_psi_reinitialize(this);
  }, ap);

  // TODO: add boundary update

  tbb::mutex gMutex;
  parallel_for( blocked_range<int>(0, NumofOverlapPkgs),[&](const blocked_range<int>& r){
    for(int num=r.begin(); num!=r.end(); ++num){
      if (interface_pkg[num]->Check_covergence_psi(this)){
        tbb::mutex::scoped_lock lock(gMutex);
        tag_convergence = 1;
      }
    }
  }, ap);
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Reinitialize_psi finished\n";
  out.close();
#endif
}
//-----------------------------------------------------
// Get_psi_max
//-----------------------------------------------------
void Levelset::Get_psi_max(communicator &world)
{
  static affinity_partitioner ap;
    
  int NumofOverlapPkgs = interface_pkg.size();
  maximum_psi = 0.;
  
  parallel_for( blocked_range<int>(0, NumofOverlapPkgs),[&](const blocked_range<int>& r){
    for(int num=r.begin(); num!=r.end(); ++num)
      interface_pkg[num]->Get_psi_max(this);
  }, ap);
  
  Parallel_get_psi_max(NumofOverlapPkgs, interface_pkg, this);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Get_psi_max finished\n";
  out.close();
#endif
}
//-----------------------------------------------------
// Redistribute_curv
//-----------------------------------------------------
void Levelset::Redistribute_curv(communicator &world)
{
  static affinity_partitioner ap;
    
  int NumofOverlapPkgs = interface_pkg.size();
  
  parallel_for( blocked_range<int>(0, NumofOverlapPkgs),[&](const blocked_range<int>& r){
    for(int num=r.begin(); num!=r.end(); ++num)
      interface_pkg[num]->Redistribute_curv(this);
  }, ap);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Redistribute_curv finished\n";
  out.close();
#endif
}
//-----------------------------------------------------
// Recalculate_max_curv
//-----------------------------------------------------
void Levelset::Recalculate_max_curv(communicator &world)
{
  static affinity_partitioner ap;
    
  parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int i_level = lset_pkg[i]->level;
      lset_pkg[i]->Get_global_scale(lset_level_info[i_level]);
    }
  }, ap);
  Parallel_get_global_scale(lset_pkg.size(), lset_pkg, this);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Recalculate_max_curv finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Smoothing_curvarure_field
//--------------------------------------------------
void Levelset::Smoothing_curvarure_field(communicator &world)
{
  static affinity_partitioner ap;
  //Extend state values among different nodes
  int iterate = 0;
  int n_iter = 50;
  int n_out_ = 50;
  int tag_convergence = 0;
  int NumofOverlapPkgs = interface_pkg.size();

  if (world.rank() == 0){
    cout<<"<<<<< Smoothing curvature for "<<n_iter<<" iterations .....:"<<endl;
    cout<<"<<<<<The max., min. curvature are        : "<<right<<setw(10)<<maximum_curv           <<", "<<right<<setw(10)<<minimum_curv<<endl;
  }
  for(int n=0; n<= n_iter; n++){
    tag_convergence = 0;

    Extend_curvature (world);

    Smooth_curvature (tag_convergence, world);
    
    iterate ++;
   
    if(iterate%n_out_ == 0){
      if (world.rank() == 0)
        cout<<"<<< Current iteration number : "<<iterate<<endl;
      n_out += n_out_;
//        Output_level_set_dat(iterate, world);
      Output_level_set_vti(n_out, world);
    }
  }
  
  parallel_for( blocked_range<int>(0, int(NumofOverlapPkgs)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int i_level = lset_pkg[i]->level;
      interface_pkg[i]->Get_global_scale(lset_level_info[i_level]);
    }
  }, ap);
  
  Parallel_get_global_scale(lset_pkg.size(), lset_pkg, this);
  
  if (world.rank() == 0){
    cout<<"<<<<< Curvature field smoothed."<<endl;
    cout<<"<<<<<The max., min. curvature are        : "<<right<<setw(10)<<maximum_curv           <<", "<<right<<setw(10)<<minimum_curv<<endl;
    cout<<"**********************************************************\n";
  }
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Smoothing_curvarure_field finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Extend_curvature
//--------------------------------------------------
void Levelset::Extend_curvature(communicator &world)
{
  static affinity_partitioner ap;

  //Extend state values among different nodes
  for(int n=0; n<20; n++){

    int NumofOverlapPkgs = interface_pkg.size();
  
    parallel_for( blocked_range<int>(0, NumofOverlapPkgs),[&](const blocked_range<int>& r){
      for(int num=r.begin(); num!=r.end(); ++num)
        interface_pkg[num]->Extend_curvature_for_smooth(this);
    }, ap);

    parallel_for( blocked_range<int>(0, NumofOverlapPkgs),[&](const blocked_range<int>& r){
      for(int num=r.begin(); num!=r.end(); ++num)
        interface_pkg[num]->Update_curvature_for_smooth(this);
    }, ap);
  }
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Extend_curvature finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Smooth_curvature
//--------------------------------------------------
void Levelset::Smooth_curvature(int &tag_convergence, communicator &world)
{
  static affinity_partitioner ap;
  int NumofOverlapPkgs = interface_pkg.size();

//      //set boundary
//     int num_bc_pkg = lset_bc_pkg.size();
//     parallel_for( blocked_range<int>(0, num_bc_pkg),[&](const blocked_range<int>& r){
//       for(int i=r.begin(); i!=r.end(); ++i){
//         int i_level = lset_bc_pkg[i]->level;
//         lset_bc_pkg[i]->Boundary_curv(lset_level_info[i_level]);
//       }
//     }, ap);

  // backup
  parallel_for( blocked_range<int>(0, int(NumofOverlapPkgs)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int i_level = lset_pkg[i]->level;
      interface_pkg[i]->Extend_curv_backup(lset_level_info[i_level]);
    }
  }, ap);    

  //get increments
  parallel_for( blocked_range<int>(0, int(NumofOverlapPkgs)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int i_level = lset_pkg[i]->level;
      interface_pkg[i]->Extend_curv(lset_level_info[i_level]);
    }
  }, ap);

  //Update values
  parallel_for( blocked_range<int>(0, int(NumofOverlapPkgs)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int i_level = lset_pkg[i]->level;
      interface_pkg[i]->Update_extend_curv_for_smoothing(lset_level_info[i_level]);
    }
  }, ap);

  // //judge the convergence
  // tbb::mutex gMutex;
  // parallel_for( blocked_range<int>(0, int(NumofOverlapPkgs)),[&](const blocked_range<int>& r){
  //   for(int i=r.begin(); i!=r.end(); ++i){
  //     if (interface_pkg[i]->Check_covergence_curv(this)){
  //       tbb::mutex::scoped_lock lock(gMutex);
  //       tag_convergence = 1;
  //     }
  //   }
  // }, ap);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Smooth_curvature finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Calculate_global_effective_curv
//--------------------------------------------------
void Levelset::Calculate_global_effective_curv(communicator &world)
{
  static affinity_partitioner ap;
  //Extend state values among different nodes
  int tag_convergence = 0;
  int iterate = 0;

  if (world.rank() == 0)
    cout<<"<<<<< The global effective curvature is started calculating:"<<endl;

  for(int n=0; n<= 100000; n++){
  // for(int n=0; n<= 2000; n++){
    tag_convergence = 0;
     //set boundary
    int num_bc_pkg = lset_bc_pkg.size();
    parallel_for( blocked_range<int>(0, num_bc_pkg),[&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        int i_level = lset_bc_pkg[i]->level;
        lset_bc_pkg[i]->Boundary_curv(lset_level_info[i_level]);
      }
    }, ap);
    // backup
    parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        int i_level = lset_pkg[i]->level;
        lset_pkg[i]->Extend_curv_backup(lset_level_info[i_level]);
      }
    }, ap);
    //get increments
    parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        int i_level = lset_pkg[i]->level;
        lset_pkg[i]->Extend_curv(lset_level_info[i_level]);
      }
    }, ap);
    //Update values
    parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        lset_pkg[i]->Update_extend_curv(this);
      }
    }, ap);
    //judge the convergence
    tbb::mutex gMutex;
    parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        if (lset_pkg[i]->Check_covergence_curv(this)){
          tbb::mutex::scoped_lock lock(gMutex);
          tag_convergence = 1;
        }
      }
    }, ap);
    
    int NumofOverlapPkgs = interface_pkg.size();
    parallel_for( blocked_range<int>(0, NumofOverlapPkgs),[&](const blocked_range<int>& r){
      for(int num=r.begin(); num!=r.end(); ++num)
        interface_pkg[num]->Copy_curvature(this);
    }, ap);
    
    iterate ++;
    
    int n_out_ = 200;
    if(iterate%n_out_ == 0){
      if (world.rank() == 0)
        cout<<"<<< Current iteration number : "<<iterate<<endl;
      n_out += n_out_;
//        Output_level_set_dat(iterate, world);
      Output_level_set_vti(n_out, world);
    }

    if(tag_convergence == 0){
      if (world.rank() == 0){
        cout<<"<<<<< The global effective curvature is calculated with iteration : "<<iterate<<endl;
        cout<<"**********************************************************\n";
      }
      break;
    }
  }
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Calculate_global_effective_curv finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Calculate_total_volume_mass
//--------------------------------------------------
void Levelset::Calculate_total_volume_mass(communicator &world)
{
  total_volume       = 0.;
  total_area         = 0.;
  total_length       = 0.;
  total_mass_segment = 0.;
  total_mass_surface = 0.;
  total_mass         = 0.;
  
  minimum_dl = fac_minimum_dl*lset_level_info[0]->dl;
  middel_dl  = fac_middel_dl*lset_level_info[0]->dl;
  maximum_dl = fac_maximum_dl*lset_level_info[0]->dl;

  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(glbl_total_num_pkg)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      lset_pkg[i]->Calculate_volume_mass(this);
    }
  }, ap);
  
  Parallel_get_total_volume_mass(lset_pkg.size(), lset_pkg, this);
  
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<<The max., mid. and min. scales are  : "<<right<<setw(10)<<maximum_dl             <<", "<<right<<setw(10)<<middel_dl   <<", "<<right<<setw(10)<<minimum_dl<<endl;
    cout<<"<<<<<The max., min. curvature are        : "<<right<<setw(10)<<maximum_curv           <<", "<<right<<setw(10)<<minimum_curv<<endl;
    cout<<"<<<<<The max.artificial curvature is     : "<<right<<setw(10)<<maximum_curv_artificial<<endl;
    cout<<"<<<<<The max. psi is                     : "<<right<<setw(10)<<maximum_psi            <<endl;
    cout<<"<<<<<The total mass and volume are       : "<<right<<setw(10)<<total_mass             <<", "<<right<<setw(10)<<total_volume<<endl;
    cout<<"<<<<<The total surf. mass and area are   : "<<right<<setw(10)<<total_mass_surface     <<", "<<right<<setw(10)<<total_area  <<endl;
    if (num_segment > 0)
    cout<<"<<<<<The total seg. mass and length are  : "<<right<<setw(10)<<total_mass_segment     <<", "<<right<<setw(10)<<total_length<<endl;
    cout<<"<<<<<The maximum dist                    : "<<right<<setw(10)<<maximum_phi            <<endl;
    cout<<"**********************************************************\n";
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Calculate_total_volume_mass finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Output levelset infomation
//--------------------------------------------------
void Levelset::Output_level_set_dat(int n, communicator &world)
{
    FILE    *fp;
    char    filename[256];
    sprintf(filename,"%s%d%s%d%s","./outdata/levelset_out.",n,".",world.rank(),".plt");

    ofstream out(filename, ios::trunc);
    
  //   for (int i = 0; i < glbl_total_num_pkg; i++){
  //     lset_pkg[i]->Output_level_set(this, world, filename, n);
  //   }
  //   
  //   for (int i = 0; i < lset_bc_pkg.size(); i++){
  //     lset_bc_pkg[i]->Output_level_set(this, world, filename, n);
  //   }
    
    lset_level_info[0]->Output_level_set(this, world, filename, n);

    out.close();

#ifdef _TEST_
  char    filename1[256];
  sprintf(filename1,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out1(filename1, ios::app);
  out1<<"<<<<<Output_level_set_dat finished\n";
  out1.close();
#endif
}
//--------------------------------------------------
// Output levelset infomation
//--------------------------------------------------
void Levelset::Output_level_set_vti(int n, communicator &world)
{
    mesh_writer.Output_lset(n, lset_level_info[0], world);

#ifdef _TEST_
  char    filename1[256];
  sprintf(filename1,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out1(filename1, ios::app);
  out1<<"<<<<<Output_level_set_dat finished\n";
  out1.close();
#endif
}
//--------------------------------------------------
// Save_lset_restart_file
//--------------------------------------------------
void Levelset::Save_lset_restart_file(char* filename, SOLVER *sph, communicator &world)
{
  ofstream out(filename, ios::app | ios::binary);
  
  out.flush();
  
  // out<<"//<<<<restart parameters: convergence_error | num_singularity | coordinate of sigularities | num_segments | start and end coordinates of segments | fac_maximum_dl | fac_middel_dl | fac_minimum_dl | maximum_curv_artificial | total_mass | total_mass_surface | total_mass_segment | total_volume | total_area | total_length | maximum_phi | maximum_psi | maximum_dl | middel_dl | minimum_dl | maximum_curv | minimum_curv | glbl_num_interface_cell"<<"\n";
  
  // out<<convergence_error<<" "<<num_singularity<<"\n";
       
  // for (int i= 0; i < num_singularity; i++ ){
  //   out<<singularity[i].i<<" "<<singularity[i].j<<" "<<singularity[i].k<<"\n";
  // }

  // out<<num_segment<<"\n";
       
  // for (int i= 0; i < num_segment; i++ ){
  //   out<<segment[i].first.i<<" "<<segment[i].first.j<<" "<<segment[i].first.k<<"\n";
  //   out<<segment[i].second.i<<" "<<segment[i].second.j<<" "<<segment[i].second.k<<"\n";
  // }
  
  // out<<fac_maximum_dl<<" "<<fac_middel_dl<<" "
  //      <<fac_minimum_dl<<" "<<maximum_curv_artificial<<"\n"
  //      <<total_mass<<" "<<total_mass_surface<<" "<<total_mass_segment<<"\n"
  //      <<total_volume<<" "<<total_area<<" "<<total_length<<"\n"
  //      <<maximum_phi<<" "<<maximum_psi<<"\n"
  //      <<maximum_dl<<" "<<middel_dl<<" "<<minimum_dl<<"\n"
  //      <<maximum_curv<<" "<<minimum_curv<<"\n"
  //      <<glbl_num_interface_cell<<"\n";

  out.write(reinterpret_cast<char*>(&convergence_error), sizeof(Real));
  out.write(reinterpret_cast<char*>(&num_singularity)  , sizeof(int));
       
  for (int i= 0; i < num_singularity; i++ ){
    out.write(reinterpret_cast<char*>(&singularity[i].i), sizeof(Real));
    out.write(reinterpret_cast<char*>(&singularity[i].j), sizeof(Real));
    out.write(reinterpret_cast<char*>(&singularity[i].k), sizeof(Real));
  }

  out.write(reinterpret_cast<char*>(&num_segment), sizeof(int));
       
  for (int i= 0; i < num_segment; i++ ){
    out.write(reinterpret_cast<char*>(&segment[i].first.i) , sizeof(Real));
    out.write(reinterpret_cast<char*>(&segment[i].first.j) , sizeof(Real));
    out.write(reinterpret_cast<char*>(&segment[i].first.k) , sizeof(Real));
    out.write(reinterpret_cast<char*>(&segment[i].second.i), sizeof(Real));
    out.write(reinterpret_cast<char*>(&segment[i].second.j), sizeof(Real));
    out.write(reinterpret_cast<char*>(&segment[i].second.k), sizeof(Real));
  }
  
  out.write(reinterpret_cast<char*>(&fac_maximum_dl)         , sizeof(Real));
  out.write(reinterpret_cast<char*>(&fac_middel_dl)          , sizeof(Real));
  out.write(reinterpret_cast<char*>(&fac_minimum_dl)         , sizeof(Real));
  out.write(reinterpret_cast<char*>(&maximum_curv_artificial), sizeof(Real));
  out.write(reinterpret_cast<char*>(&total_mass)             , sizeof(Real));
  out.write(reinterpret_cast<char*>(&total_mass_surface)     , sizeof(Real));
  out.write(reinterpret_cast<char*>(&total_mass_segment)     , sizeof(Real));
  out.write(reinterpret_cast<char*>(&total_volume)           , sizeof(Real));
  out.write(reinterpret_cast<char*>(&total_area)             , sizeof(Real));
  out.write(reinterpret_cast<char*>(&total_length)           , sizeof(Real));
  out.write(reinterpret_cast<char*>(&maximum_phi)            , sizeof(Real));
  out.write(reinterpret_cast<char*>(&maximum_psi)            , sizeof(Real));
  out.write(reinterpret_cast<char*>(&maximum_dl)             , sizeof(Real));
  out.write(reinterpret_cast<char*>(&middel_dl)              , sizeof(Real));
  out.write(reinterpret_cast<char*>(&minimum_dl)             , sizeof(Real));
  out.write(reinterpret_cast<char*>(&maximum_curv)           , sizeof(Real));
  out.write(reinterpret_cast<char*>(&minimum_curv)           , sizeof(Real));
  out.write(reinterpret_cast<char*>(&glbl_num_interface_cell), sizeof(int ));

  out.flush();
  
  lset_level_info[0]->Write_level_set_rstfile(this, world, filename);
}
//--------------------------------------------------
// Save_lset_restart_file
//--------------------------------------------------
void Levelset::Load_lset_restart_file(ifstream &load, SOLVER *sph, communicator &world)
{

  load.read(reinterpret_cast<char*>(&convergence_error), sizeof(Real));
  load.read(reinterpret_cast<char*>(&num_singularity)  , sizeof(int));
  
  singularity = new my_real[num_singularity];

  for (int i= 0; i < num_singularity; i++ ){
    load.read(reinterpret_cast<char*>(&singularity[i].i), sizeof(Real));
    load.read(reinterpret_cast<char*>(&singularity[i].j), sizeof(Real));
    load.read(reinterpret_cast<char*>(&singularity[i].k), sizeof(Real));
  }

  load.read(reinterpret_cast<char*>(&num_segment), sizeof(int));
  
  segment = new std::pair <my_real, my_real>[num_segment];

  for (int i= 0; i < num_segment; i++ ){
    load.read(reinterpret_cast<char*>(&segment[i].first.i) , sizeof(Real));
    load.read(reinterpret_cast<char*>(&segment[i].first.j) , sizeof(Real));
    load.read(reinterpret_cast<char*>(&segment[i].first.k) , sizeof(Real));
    load.read(reinterpret_cast<char*>(&segment[i].second.i), sizeof(Real));
    load.read(reinterpret_cast<char*>(&segment[i].second.j), sizeof(Real));
    load.read(reinterpret_cast<char*>(&segment[i].second.k), sizeof(Real));
  }
  
  load.read(reinterpret_cast<char*>(&fac_maximum_dl)         , sizeof(Real));
  load.read(reinterpret_cast<char*>(&fac_middel_dl)          , sizeof(Real));
  load.read(reinterpret_cast<char*>(&fac_minimum_dl)         , sizeof(Real));
  load.read(reinterpret_cast<char*>(&maximum_curv_artificial), sizeof(Real));
  load.read(reinterpret_cast<char*>(&total_mass)             , sizeof(Real));
  load.read(reinterpret_cast<char*>(&total_mass_surface)     , sizeof(Real));
  load.read(reinterpret_cast<char*>(&total_mass_segment)     , sizeof(Real));
  load.read(reinterpret_cast<char*>(&total_volume)           , sizeof(Real));
  load.read(reinterpret_cast<char*>(&total_area)             , sizeof(Real));
  load.read(reinterpret_cast<char*>(&total_length)           , sizeof(Real));
  load.read(reinterpret_cast<char*>(&maximum_phi)            , sizeof(Real));
  load.read(reinterpret_cast<char*>(&maximum_psi)            , sizeof(Real));
  load.read(reinterpret_cast<char*>(&maximum_dl)             , sizeof(Real));
  load.read(reinterpret_cast<char*>(&middel_dl)              , sizeof(Real));
  load.read(reinterpret_cast<char*>(&minimum_dl)             , sizeof(Real));
  load.read(reinterpret_cast<char*>(&maximum_curv)           , sizeof(Real));
  load.read(reinterpret_cast<char*>(&minimum_curv)           , sizeof(Real));
  load.read(reinterpret_cast<char*>(&glbl_num_interface_cell), sizeof(int ));

  // char buffer[500];
  // load.getline (buffer,500);
  // load.getline (buffer,500);

  // load>>convergence_error>>num_singularity;
  
  // singularity = new my_real[num_singularity];
       
  // for (int i= 0; i < num_singularity; i++ ){
  //   load>>singularity[i].i>>singularity[i].j>>singularity[i].k;
  // }

  // load>>num_segment;
  
  // segment = new std::pair <my_real, my_real>[num_segment];
       
  // for (int i= 0; i < num_segment; i++ ){
  //   load>>segment[i].first.i>>segment[i].first.j>>segment[i].first.k;
  //   load>>segment[i].second.i>>segment[i].second.j>>segment[i].second.k;
  // }
  
  // load>>fac_maximum_dl>>fac_middel_dl
  //      >>fac_minimum_dl>>maximum_curv_artificial
  //      >>total_mass>>total_mass_surface>>total_mass_segment
  //      >>total_volume>>total_area>>total_length
  //      >>maximum_phi>>maximum_psi
  //      >>maximum_dl>>middel_dl>>minimum_dl
  //      >>maximum_curv>>minimum_curv>>glbl_num_interface_cell;
       
  Initialize(sph, world);
       
  lset_level_info[0]->Load_level_set_rstfile(load, this, world);
  
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<<The max., mid. and min. scales are  : "<<right<<setw(10)<<maximum_dl        <<", "<<right<<setw(10)<<middel_dl   <<", "<<right<<setw(10)<<minimum_dl<<endl;
    cout<<"<<<<<The max., min. curvature are        : "<<right<<setw(10)<<maximum_curv      <<", "<<right<<setw(10)<<minimum_curv<<endl;
    cout<<"<<<<<The max. psi is                     : "<<right<<setw(10)<<maximum_psi       <<endl;
    cout<<"<<<<<The total mass and volume are       : "<<right<<setw(10)<<total_mass        <<", "<<right<<setw(10)<<total_volume<<endl;
    cout<<"<<<<<The total surf. mass and area are   : "<<right<<setw(10)<<total_mass_surface<<", "<<right<<setw(10)<<total_area  <<endl;
    if (num_segment > 0)
    cout<<"<<<<<The total seg. mass and length are  : "<<right<<setw(10)<<total_mass_segment<<", "<<right<<setw(10)<<total_length<<endl;
    cout<<"<<<<<The maximum dist                    : "<<right<<setw(10)<<maximum_phi <<endl;
    cout<<"**********************************************************\n";
  }
}
#if defined (_MPI_) && (defined (_CVP_LSET_) || defined (_CVP_LSET_INIT_))
//--------------------------------------------------
// CVP_for_initial_partitioning
//--------------------------------------------------
void Levelset::CVP_for_initial_partitioning (std::vector<my_real> & vp_coords, std::vector<Real> &vp_scale, communicator &world)
{
  lset_voronoi.Initialize (this, 1, world);

  lset_voronoi.Partitioning_for_initial_VP_distribution(this, world);

  lset_voronoi.Get_vp_position_and_hmin (vp_coords, vp_scale, world);

#ifdef _TEST_
  char    filename1[256];
  sprintf(filename1,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out1(filename1, ios::app);
  out1<<"<<<<<CVP_for_initial_partitioning finished\n";
  out1.close();
#endif
}
//--------------------------------------------------
// Save_VP_coord_hmin
//--------------------------------------------------
void Levelset::Save_VP_coord_hmin (std::vector<my_real> & vp_coords, std::vector<Real> &vp_scale, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s","./rstfile/VP_coords.csv");

  ofstream out_(filename, ios::trunc);

  out_<<"\"color\",\"coordX\",\"coordY\",\"coordZ\",\"hmin\""<<endl;

  int nVP = vp_coords.size();
  for (int i = 0; i < nVP; ++i)
  {
    out_<<i<<","<<vp_coords[i].i<<","<<vp_coords[i].j<<","<<vp_coords[i].k<<","<<vp_scale[i]<<endl;
  }
  out_.close();

#ifdef _TEST_
  char    filename1[256];
  sprintf(filename1,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out1(filename1, ios::app);
  out1<<"<<<<<CVP_for_initial_partitioning finished\n";
  out1.close();
#endif
}
#endif
