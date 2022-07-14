#include "glbfunc.h"
#include "sph_mesh_generation.h"

//--------------------------------------------------
// SPH mesh generation class
// initialization
//--------------------------------------------------
SPH_mesh_generation::SPH_mesh_generation(Initialization &Ini, communicator &world)
{
  count_full = 0;
  fill = true;
  current_stage = FILL;
  iterate_num    = 0;
  i_iter         = 0;
  max_v          = 1.e-7;
  glbl_timestep  = 0.0;
  glbl_max_v     = 0.0;
  timestep_shift = 0.0;
  run_time       = 0.0;
  safe_guard     = 1.5;

  angle_thres = PI/4.0;
  
  max_num_iteration = 0;
  
  local_id_record = 0;
  
  total_num_particle_mesh = 0.;
  glbl_total_num_particle_mesh = 0.;
  
  total_num_particle_surface = 0;
  glbl_total_num_particle_surface = 0.;
  
  total_num_particle_singularity = 0.;
  glbl_total_num_particle_singularity = 0.;
  
  total_num_particle_segment = 0.;
  glbl_total_num_particle_segment = 0.;
  
  find_surface_particle = false;
  find_singularity_particle = false;
  find_segment_particle = false;

  WVT_coeff  = 0.;
  APD_coeff  = 0.;

  damp_coeff_start   = 0.;
  damp_coeff_end     = damp_coeff_start;
  damp_ramping_start = Ini.start_time;
  damp_ramping_end   = Ini.end_time;

  nu_coeff_start   = 0.1;
  nu_coeff_end     = nu_coeff_start;
  nu_ramping_start = Ini.start_time;
  nu_ramping_end   = Ini.end_time;

  v_reini_start    = 1;
  v_reini_end      = 1;
  v_reini_change   = int(Ini.end_time);

  tet_delete_thres = 10.;

  t_for_post = Ini.end_time;
  n_post     = 0;

  #ifdef _MPI_
  icount = 0;
  sigma_coeff_start   = 4.;
  sigma_coeff_end     = 1.;
  sigma_ramping_start = Ini.start_time;
  sigma_ramping_end   = Ini.end_time;
  #endif
  
  my_set_const (body_force, 0.);
  
  #ifdef _READ_SDF_
  my_set_const (cells_to_cut, 0);
  #endif

  particle_singularity.clear();
  particle_segment .clear();
  particle_surface.clear();
  particle_bound.clear();
  dummy_particle.clear();
  ghost_particle.clear();
  real_particle.clear();
  num_particles_per_singularity.clear();

  time_for_total                   = new Real[world.size()+1];
  time_for_partition               = new Real[world.size()+1];
  time_for_graph                   = new Real[world.size()+1];
  time_for_simulation              = new Real[world.size()+1];
  time_for_force_calculation       = new Real[world.size()+1];
  time_for_density                 = new Real[world.size()+1];
  time_for_rest                    = new Real[world.size()+1];
  time_for_update                  = new Real[world.size()+1];
  time_for_reset_neighbor          = new Real[world.size()+1];
  time_for_mapping                 = new Real[world.size()+1];
  time_for_bulid_local_map         = new Real[world.size()+1];
  time_for_map_particle_to_tree    = new Real[world.size()+1];
  time_for_refresh_neighbor        = new Real[world.size()+1];
  time_for_clear_cell_list         = new Real[world.size()+1];
  time_for_release_particle        = new Real[world.size()+1];
  time_for_update_every_level_info = new Real[world.size()+1];
  time_for_communication0          = new Real[world.size()+1];
  time_for_communication1          = new Real[world.size()+1];
  time_for_communication2          = new Real[world.size()+1];
  time_for_communication_partition = new Real[world.size()+1];
  time_for_communication_flag1     = new Real[world.size()+1];
  time_for_communication_flag2     = new Real[world.size()+1];
  time_for_exchangebuffer1         = new Real[world.size()+1];
  time_for_exchangebuffer2         = new Real[world.size()+1];

  for (int i = 0; i < world.size()+1; i++){
    time_for_total[i]                   = 0.;
    time_for_partition[i]               = 0.;
    time_for_graph[i]                   = 0.;
    time_for_simulation[i]              = 0.;
    time_for_force_calculation[i]       = 0.;
    time_for_density[i]                 = 0.;
    time_for_rest[i]                    = 0.;
    time_for_update[i]                  = 0.;
    time_for_reset_neighbor[i]          = 0.;
    time_for_mapping[i]                 = 0.;
    time_for_bulid_local_map[i]         = 0.;
    time_for_map_particle_to_tree[i]    = 0.;
    time_for_refresh_neighbor[i]        = 0.;
    time_for_clear_cell_list[i]         = 0.;
    time_for_release_particle[i]        = 0.;
    time_for_update_every_level_info[i] = 0.;
    time_for_communication0[i]          = 0.;
    time_for_communication1[i]          = 0.;
    time_for_communication2[i]          = 0.;
    time_for_communication_partition[i] = 0.;
    time_for_communication_flag1[i]     = 0.;
    time_for_communication_flag2[i]     = 0.;
    time_for_exchangebuffer1[i]         = 0.;
    time_for_exchangebuffer2[i]         = 0.;
  }
#ifdef _MPI_
  need_for_rebuild_graph           = 0;
  need_for_partition_count         = 0;
  need_for_rebuild_local_tree      = 0;
#endif
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<class sph initialized\n";
  out.close();
#endif
}
//--------------------------------------------------
// initial condition for cases
//--------------------------------------------------
void SPH_mesh_generation::Load_case(communicator &world)
{
  /****************************************************/
  // user define area
  /****************************************************/
  // JZ20190104::Read initial signed distance function from 
  // open-source code SDFGen
  #ifdef _READ_SDF_

  if (DIM == 2 && world.rank() == 0){
    cout<<"<<<<< ERROR!!! Unable to read 2D SDF files"<<endl;
    world.abort(-1);
  }

  // sprintf(sdf_file,"%s","./SDF_file/bunny_watertight.sdf");
  // sprintf(sdf_file,"%s","./SDF_file/bunny_watertight_dp0p75.sdf");
  // sprintf(sdf_file,"%s","./SDF_file/wing02_dp0p2.sdf");
  // sprintf(sdf_file,"%s","./SDF_file/wing02_dp0p5_outer.sdf");
  sprintf(sdf_file,"%s","./SDF_file/wing02_dp1_outer.sdf");
  // sprintf(sdf_file,"%s","./SDF_file/bunny_watertight_dp1.sdf");
  // sprintf(sdf_file,"%s","./SDF_file/Gear_Spur_16T.sdf");
  // sprintf(sdf_file,"%s","./SDF_file/2Gears_high.sdf");
  // sprintf(sdf_file,"%s","./SDF_file/tyra.sdf");
  
  ifstream load(sdf_file);
  
  if (! load.is_open())  { 
    cout << "<<<<< Error opening .sdf file in rank "<<world.rank(); 
    world.abort (-1); 
  }

  my_int grid_size;
  Real dx = 0.;

  my_real shift;
  load>>grid_size.i>>grid_size.j>>grid_size.k;
  load>>shift.i    >>shift.j    >>shift.k;
  load>>dx;
  box_l = shift;

  cells_to_cut.i = grid_size.i % ICPX;
  cells_to_cut.j = grid_size.j % ICPY;
  cells_to_cut.k = grid_size.k % ICPZ;

  grid_size = my_minus_data (grid_size, cells_to_cut);

  domain.i = dx*Real(grid_size.i);
  domain.j = dx*Real(grid_size.j);
  domain.k = dx*Real(grid_size.k);

  my_set_const (box_l, 0.);
  box_r = my_add_data (box_l, domain);

  my_set_data (particle_box  , domain);
  my_set_data (particle_box_l, box_l );
  my_set_data (particle_box_r, box_r );

  resolution.i = grid_size.i / ICPX;
  resolution.j = grid_size.j / ICPY;
  resolution.k = grid_size.k / ICPZ;

  load.close();

  #else
  my_set_const (domain, 100.0);
  my_set_const (box_l, 0.);
  my_set_const (box_r, 100.);
  my_set_const (particle_box, 100.0);
  my_set_const (particle_box_l, 0.);
  my_set_const (particle_box_r, 100.);

  resolution.i = DIM_X == 1 ? 128 : 1;
  resolution.j = DIM_Y == 1 ? 128 : 1;
  resolution.k = DIM_Z == 1 ? 128 : 1;
  #endif
  
  level_set.convergence_error = 1.e-4;
  
  level_set.num_singularity = 0;
  level_set.singularity = new my_real[level_set.num_singularity];
  level_set.num_segment = 0;
  level_set.segment = new std::pair<my_real, my_real> [level_set.num_segment];
  
  // my_real center;
  // center.i = 50.;
  // center.j = 50.;
  // center.k = 0.;

  // Real R = 45.0;
  // Real W = 7.5;
  // Real H = 30.63;
  // Real low_Y  = center.j - sqrt(R*R-W*W);
  // Real left_X = center.i - W;
  // Real left_Y = center.j + H;
  // Real right_X = center.i + W;

  // level_set.num_singularity = 4;
  // level_set.singularity = new my_real[level_set.num_singularity];
  // level_set.singularity[0].i = left_X;
  // level_set.singularity[0].j = left_Y;
  // level_set.singularity[0].k = 0.;
  
  // level_set.singularity[1].i = right_X;
  // level_set.singularity[1].j = left_Y;
  // level_set.singularity[1].k = 0.;
  
  // level_set.singularity[2].i = left_X;
  // level_set.singularity[2].j = low_Y;
  // level_set.singularity[2].k = 0.;
  
  // level_set.singularity[3].i = right_X;
  // level_set.singularity[3].j = low_Y;
  // level_set.singularity[3].k = 0.;
  
  level_set.fac_maximum_dl = 20.;
  level_set.fac_middel_dl =  5.;
  level_set.fac_minimum_dl = 2.;
  level_set.maximum_curv_artificial = 5.;

  #if defined (_MPI_) &&  defined (_CVP_LSET_INIT_)
  level_set.error_tolerance = 0.1;
  level_set.num_color_tmp = 6;
  #endif

  #ifdef _READ_SDF_
    level_set.num_artifact_region    = 0;
    level_set.radi_artifact_region   = new Real   [level_set.num_artifact_region];
    level_set.center_artifact_region = new my_real[level_set.num_artifact_region];

    // level_set.num_artifact_region    = 3;
    // level_set.radi_artifact_region   = new Real   [level_set.num_artifact_region];
    // level_set.center_artifact_region = new my_real[level_set.num_artifact_region];

    // level_set.radi_artifact_region  [0]   =  20.;
    // level_set.center_artifact_region[0].i = 115.; 
    // level_set.center_artifact_region[0].j =  53.; 
    // level_set.center_artifact_region[0].k =  13.;

    // level_set.radi_artifact_region  [1]   =  15.;
    // level_set.center_artifact_region[1].i =  93.;
    // level_set.center_artifact_region[1].j =  36.; 
    // level_set.center_artifact_region[1].k =  13.;

    // level_set.radi_artifact_region  [2]   =  16.;
    // level_set.center_artifact_region[2].i =  63.;
    // level_set.center_artifact_region[2].j =  51.;
    // level_set.center_artifact_region[2].k =  10.;
  #endif
  
  /****************************************************/
  // level set operations
  /****************************************************/
  num_level      = 1;
  Lmax           = num_level-1;
  Lmin           = num_level-1;
  
  ini_scale = domain.i / resolution.i;
  max_scale = ini_scale;
  min_scale = ini_scale;

  ini_num_cell.i = DIM_X==1 ? int(domain.i/ini_scale): 1;
  ini_num_cell.j = DIM_Y==1 ? int(domain.j/ini_scale): 1;
  ini_num_cell.k = DIM_Z==1 ? int(domain.k/ini_scale): 1;
}
//--------------------------------------------------
// Load_rst_case
//--------------------------------------------------
void SPH_mesh_generation::Load_rst_case(communicator &world)
{
  /****************************************************/
  // read restart file
  /****************************************************/
  char fn[256];
  char buffer[500];
  sprintf(fn,"%s","./rstfile/levelset.rst");
  
  ifstream load(fn, ios::binary);
  
  if (! load.is_open())  { 
    cout << "<<<<< Error opening file in rank "<<world.rank(); 
    world.abort (-1); 
  }
  
  if (world.rank() == 0)
    cout<<"<<<<< Loading restart file........."<<endl;
  
  // load.getline (buffer,500);
  int tmp_int = 0;
  Real tmp_real = 0.;
  
  // load>>tmp_int;
  load.read(reinterpret_cast<char*>(&tmp_int), sizeof(int));

  if (tmp_int != int (DIM) && world.rank() == 0){
    cout<<"<<<<<ERROR!!! dimention is wrong!!! DIM = "<<tmp_int<<endl;
  }

  // load>>domain.i>>domain.j>>domain.k
  //       >>box_l.i>>box_l.j>>box_l.k
  //       >>box_r.i>>box_r.j>>box_r.k
  //       >>particle_box.i>>particle_box.j>>particle_box.k
  //       >>particle_box_l.i>>particle_box_l.j>>particle_box_l.k
  //       >>particle_box_r.i>>particle_box_r.j>>particle_box_r.k
  //       >>resolution.i>>resolution.j>>resolution.k
  //       >>Lmin>>Lmax
  //       >>ini_num_cell.i>>ini_num_cell.j>>ini_num_cell.k
  //       >>num_level;

  load.read(reinterpret_cast<char*>(&domain.i)        , sizeof(Real));
  load.read(reinterpret_cast<char*>(&domain.j)        , sizeof(Real));
  load.read(reinterpret_cast<char*>(&domain.k)        , sizeof(Real));
  load.read(reinterpret_cast<char*>(&box_l.i)         , sizeof(Real));
  load.read(reinterpret_cast<char*>(&box_l.j)         , sizeof(Real));
  load.read(reinterpret_cast<char*>(&box_l.k)         , sizeof(Real));
  load.read(reinterpret_cast<char*>(&box_r.i)         , sizeof(Real));
  load.read(reinterpret_cast<char*>(&box_r.j)         , sizeof(Real));
  load.read(reinterpret_cast<char*>(&box_r.k)         , sizeof(Real));
  load.read(reinterpret_cast<char*>(&particle_box.i)  , sizeof(Real));
  load.read(reinterpret_cast<char*>(&particle_box.j)  , sizeof(Real));
  load.read(reinterpret_cast<char*>(&particle_box.k)  , sizeof(Real));
  load.read(reinterpret_cast<char*>(&particle_box_l.i), sizeof(Real));
  load.read(reinterpret_cast<char*>(&particle_box_l.j), sizeof(Real));
  load.read(reinterpret_cast<char*>(&particle_box_l.k), sizeof(Real));
  load.read(reinterpret_cast<char*>(&particle_box_r.i), sizeof(Real));
  load.read(reinterpret_cast<char*>(&particle_box_r.j), sizeof(Real));
  load.read(reinterpret_cast<char*>(&particle_box_r.k), sizeof(Real));
  load.read(reinterpret_cast<char*>(&resolution.i)    , sizeof(int ));
  load.read(reinterpret_cast<char*>(&resolution.j)    , sizeof(int ));
  load.read(reinterpret_cast<char*>(&resolution.k)    , sizeof(int ));
  load.read(reinterpret_cast<char*>(&Lmin)            , sizeof(int ));
  load.read(reinterpret_cast<char*>(&Lmax)            , sizeof(int ));
  load.read(reinterpret_cast<char*>(&ini_num_cell.i)  , sizeof(int ));
  load.read(reinterpret_cast<char*>(&ini_num_cell.j)  , sizeof(int ));
  load.read(reinterpret_cast<char*>(&ini_num_cell.k)  , sizeof(int ));
  load.read(reinterpret_cast<char*>(&num_level)       , sizeof(int ));

  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Simulation parameters:\n"
          <<"<<<dmain          : "<<right<<setw(10)<<domain.i        <<" "<<right<<setw(10)<<domain.j        <<" "<<right<<setw(10)<<domain.k        <<"\n"
          <<"<<<box_l          : "<<right<<setw(10)<<box_l.i         <<" "<<right<<setw(10)<<box_l.j         <<" "<<right<<setw(10)<<box_l.k         <<"\n"
          <<"<<<box_r          : "<<right<<setw(10)<<box_r.i         <<" "<<right<<setw(10)<<box_r.j         <<" "<<right<<setw(10)<<box_r.k         <<"\n"
          <<"<<<particle_box   : "<<right<<setw(10)<<particle_box.i  <<" "<<right<<setw(10)<<particle_box.j  <<" "<<right<<setw(10)<<particle_box.k  <<"\n"
          <<"<<<particle_box_l : "<<right<<setw(10)<<particle_box_l.i<<" "<<right<<setw(10)<<particle_box_l.j<<" "<<right<<setw(10)<<particle_box_l.k<<"\n"
          <<"<<<particle_box_r : "<<right<<setw(10)<<particle_box_r.i<<" "<<right<<setw(10)<<particle_box_r.j<<" "<<right<<setw(10)<<particle_box_r.k<<"\n"
          <<"<<<resolution     : "<<right<<setw(10)<<resolution.i    <<" "<<right<<setw(10)<<resolution.j    <<" "<<right<<setw(10)<<resolution.k    <<"\n"
          <<"<<<Lmin | Lmax    : "<<right<<setw(10)<<Lmin            <<" "<<right<<setw(10)<<Lmax            <<"\n"
          <<"<<<ini_num_cell   : "<<right<<setw(10)<<ini_num_cell.i  <<" "<<right<<setw(10)<<ini_num_cell.j  <<" "<<right<<setw(10)<<ini_num_cell.k  <<"\n"
          <<"<<<num_level      : "<<right<<setw(10)<<num_level       <<"\n";
    cout<<"**********************************************************\n";
  }

  level_set.Load_lset_restart_file (load, this, world);
  
  load.close();

  if (world.rank() == 0)
    cout<<"<<<<< Restart file loaded"<<endl;
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Load_rst_case finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Load_data_for_post
//--------------------------------------------------
void SPH_mesh_generation::Load_data_for_post(communicator &world)
{
  /****************************************************/
  // read postpros file
  /****************************************************/
  char list[256];
  char buffer[500];
  sprintf(list,"%s","./rstfile_particles/list.csv");

  ifstream load_list(list);
  
  if (! load_list.is_open())  { 
    cout << "<<<<< Error opening list.csv in rank "<<world.rank(); 
    world.abort (-1); 
  }

  int num_steps = 0;
  while (load_list.getline (buffer,500))
    num_steps++;
  num_steps--;

  load_list.clear();
  load_list.seekg(0, ios::beg);
  
  int num_of_rank = 0;
  load_list>>num_of_rank;
  cout<<"num_of_rank: "<<num_of_rank<<endl;

  if (num_of_rank <= 0){
    cout<<"<<<<< ERROR number of ranks!!!"<<endl;
  }

  std::vector<int> step;       step.clear();
  std::vector<int> npoints; npoints.clear();

  for (int i = 0; i < num_steps; ++i)
  {
    int tmp_step = 0;
    int tmp_np  = 0;
    load_list>>tmp_step>>tmp_np;
       step.push_back(tmp_step);
    npoints.push_back(tmp_np);
  }

  if (step.size() != num_steps){
    cout<<"ERROR!!! number of post step is wrong!!!!"<<endl;
    world.abort(-1);
  }

  for (int i = 0; i < num_steps; ++i){

    cout<<"<<<<< Post_processing for step:   "<<step[i]<<endl;

    run_time = step[i];

    particle = new p_Particle[npoints[i]];
    total_num_particle = npoints[i];
    glbl_total_num_particle = total_num_particle;
    std::vector<p_Particle> particle_total; particle_total.clear();

    int pid = 0;
    for(int ip = 0; ip < total_num_particle; ip++){
      particle[pid] = particlepool.malloc();
      particle_total.push_back(particle[pid]);
      particle[pid]->Initialize (this, pid, REAL_PARTICLE, world.rank());
      pid++;
    }

    Asign_local_particle_index (world);

    pid = 0;
    for (int irank = 0; irank < num_of_rank; ++irank)
    {
      char fn[256];
      sprintf(fn,"%s%d%s%d%s","./rstfile_particles/step_",step[i],"_rank_",irank,".pst");

      ifstream load(fn, ios::binary);

      if (! load.is_open())  {
        cout << "<<<<< Error opening post file: "<<fn<<endl;
        world.abort (-1);
      }

      cout<<">>>>> Reading partices from rank :   "<<irank<<endl;

      int npp = 0;
      load.read(reinterpret_cast<char*>(&npp), sizeof(int));

      for (int ipp = 0; ipp < npp; ++ipp)
      {
        my_real coord_;
        load.read(reinterpret_cast<char*>(&coord_.i), sizeof(Real));
        load.read(reinterpret_cast<char*>(&coord_.j), sizeof(Real));
        load.read(reinterpret_cast<char*>(&coord_.k), sizeof(Real));

        my_set_data(particle[pid]->coord, coord_);
        particle[pid]->color =  irank;
        pid++;
      }
    }

    if (pid != total_num_particle){
      cout<<"<<<<< ERROR in number of particles read from .pst files"<<endl;
      world.abort(-1);
    }

    Get_particle_info (world);

    Set_particle_scale (world);

    Find_and_map_particle (world);

    if (DIM == 2){

      mesh.Remove_all_edges (world);
      mesh.Remove_all_verts (world);
      
      mesh.Initialize (particle_total.size(), world);

    }else if (DIM == 3){

      #ifdef _GEN_DUMMY_PART_
        Reconstruce_dummy_particles (particle_total, world);
      #endif

      // construct faces
      mesh.Reconstruct_tets (particle_total, this, world);

      Reset_particle_array (world);
      
      Check_total_number_of_particle (world);
    }

    // Output mesh
    mesh_writer.Output_mesh (n_post, particle_total, this, world);

    // Output report
    if (DIM == 3)  {
      mesh.Write_tet_quality_plt (n_post, run_time, world);

      #ifdef _GEN_DUMMY_PART_
        Free_dummy_particles (particle_total, world);
      #endif

    }
    n_post ++;
    particle_total.clear();

    for(int ip = 0; ip < total_num_particle; ip++){
      particle[ip]->Clearup();
      particlepool.free(particle[ip]);
    }
    delete [] particle;
  }

  cout<<"<<<<< Post_processing finished"<<endl;

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Load_data_for_post finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// initial condition for cases
//--------------------------------------------------
void SPH_mesh_generation::Initialize_case(communicator &world)
{
#ifdef _MPI_
  glbl_total_num_color = world.size();
#else
  glbl_total_num_color = 1;
#endif
  
#if !defined(_READ_TARGET_FIELD_) && !defined(_READ_TARGET_FIELD_UNCHANGE_) && !defined(_READ_TARGET_FIELD_FOR_POST_)
  /****************************************************/
  // Non-restart
  /****************************************************/

  if (world.rank() == 0){

    Load_case(world);

    level_set.Initialize (this, world);
    
    level_set.Calculate_target_density_field(this, world);

    Get_num_particles (world);

    #if defined (_MPI_) &&  defined (_CVP_LSET_INIT_)
      std::vector<my_real> vp_coords; vp_coords.clear();
      std::vector<Real>    vp_scale ; vp_scale.clear();

      level_set.CVP_for_initial_partitioning (vp_coords, vp_scale, world);

      level_set.Save_VP_coord_hmin (vp_coords, vp_scale, world);
    #endif

    Save_restart_file(world);
  }

  world.barrier();
  world.abort(-1);

#else
  /****************************************************/
  // Restart
  /****************************************************/
  Load_rst_case(world);

  world.barrier();

  #if defined (_READ_TARGET_FIELD_FOR_POST_)
  // do check
  if (world.rank()==0 && world.size() > 1){
    cout<<"<<<<< ERROR: please set ntasks=1 for post only!!!"<<endl;
    world.abort(-1);    
  }

  Load_data_for_post(world);

  world.barrier();
  world.abort(-1);
  #endif
  
  #if defined(_READ_TARGET_FIELD_)
  /****************************************************/
  // change the target density field by 
  // overwriting the parameters accordingly
  /****************************************************/
    // level_set.Output_level_set_dat(0, world);
    level_set.fac_maximum_dl = 16.0;
    level_set.fac_middel_dl =  5.0;
    level_set.fac_minimum_dl = 2.0;
    level_set.maximum_curv_artificial = 1.0;
    level_set.maximum_curv_artificial = AMIN1 (level_set.maximum_curv_artificial, level_set.maximum_curv);
    
    level_set.Calculate_total_volume_mass(world);
  #endif

  level_set.Get_extended_cell_tags(world);

  if (world.rank() == 0)
    level_set.Output_level_set_vti(level_set.n_out++, world);
  // level_set.Output_level_set_dat(1, world);
  // world.barrier();
  // world.abort(-1);
  
#endif
  /****************************************************/
  // particle mesh generation solver initialization
  /****************************************************/
  max_num_iteration = 1e6;
  fill_coeff = 2.0;
  max_count_full = 500000000;
  max_count_relax = 500;
  
  WVT_coeff  = 0.0;
  APD_coeff  = 0.0;

  t_for_post = 0.0*t_for_post;
  angle_thres = PI*0.4;
  
  nu_coeff_start   = 0.1;
  nu_coeff_end     = 0.1;
  nu_ramping_start = nu_ramping_end*0.5;
  nu_ramping_end   = nu_ramping_end*0.52;
  
  damp_coeff_start   = 0.;
  damp_coeff_end     = 0.;
  damp_ramping_start = damp_ramping_end*0.3;
  damp_ramping_end   = damp_ramping_end*0.35;
  
  v_reini_start      = 1;
  v_reini_end        = 1;
  v_reini_change     = int(v_reini_change*0.5);

  tet_delete_thres   = 10.0;
  
  #ifdef _MPI_
  // control of ST force
  sigma_coeff_start   = 4.;
  sigma_coeff_end     = 1.;
  sigma_ramping_start = sigma_ramping_end*0.18;
  sigma_ramping_end   = sigma_ramping_end*0.23;

    // JZ20190124::different initialization method for creating VP positions
    #if defined (_CVP_LSET_)

      level_set.error_tolerance = 0.05;

    #elif defined (_CLOSE_PACKING_)

      // JZ20181221::using CVP for the initial partitioning
      my_real CVP_box_l;
      my_real CVP_box;
      my_set_const (CVP_box_l, 20.);
      my_set_const (CVP_box  , 60.);
      // CVP_box_l.i = 75.;
      // CVP_box_l.j = 30.;
      // CVP_box_l.k = 25.;
      // CVP_box.i = 50.;
      // CVP_box.j = 60.;
      // CVP_box.k = 70.;
      
      int dimX = 2;
      int dimY = 3;
      int dimZ = 1;
      
      if (dimX*dimY*dimZ != world.size()){
        if (world.rank() == 0) cout<<"<<<<< ERROR number of dimX/Y/Z been set"<<endl;
        world.barrier();
        world.abort(-1);
      }
    #endif
  #endif

  Get_num_particles (world);

  if (world.size() > 1)
    Recalculate_local_particle_number (world);

  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< APD_coeff             : "<<right<<setw(10)<<APD_coeff         <<endl;
    cout<<"<<<<< WVT_coeff             : "<<right<<setw(10)<<WVT_coeff         <<endl;
    
    cout<<"<<<<< damp_coeff_start      : "<<right<<setw(10)<<damp_coeff_start  <<endl;
    cout<<"<<<<< damp_coeff_end        : "<<right<<setw(10)<<damp_coeff_end    <<endl;
    cout<<"<<<<< damp_ramping_start    : "<<right<<setw(10)<<damp_ramping_start<<endl;
    cout<<"<<<<< damp_ramping_end      : "<<right<<setw(10)<<damp_ramping_end  <<endl;
    
    cout<<"<<<<< nu_coeff_start        : "<<right<<setw(10)<<nu_coeff_start    <<endl;
    cout<<"<<<<< nu_coeff_end          : "<<right<<setw(10)<<nu_coeff_end      <<endl;
    cout<<"<<<<< nu_ramping_start      : "<<right<<setw(10)<<nu_ramping_start  <<endl;
    cout<<"<<<<< nu_ramping_end        : "<<right<<setw(10)<<nu_ramping_end    <<endl;
    
    cout<<"<<<<< v_reini_start         : "<<right<<setw(10)<<v_reini_start     <<endl;
    cout<<"<<<<< v_reini_end           : "<<right<<setw(10)<<v_reini_end       <<endl;
    cout<<"<<<<< v_reini_change        : "<<right<<setw(10)<<v_reini_change    <<endl;
    #ifdef _MPI_
    cout<<"<<<<< sigma_coeff_start     : "<<right<<setw(10)<<sigma_coeff_start  <<endl;
    cout<<"<<<<< sigma_coeff_end       : "<<right<<setw(10)<<sigma_coeff_end    <<endl;
    cout<<"<<<<< sigma_ramping_start   : "<<right<<setw(10)<<sigma_ramping_start<<endl;
    cout<<"<<<<< sigma_ramping_end     : "<<right<<setw(10)<<sigma_ramping_end  <<endl;    
    #endif
    cout<<"<<<<< max_num_iteration     : "<<right<<setw(10)<<max_num_iteration  <<endl;
    cout<<"<<<<< count_full            : "<<right<<setw(10)<<count_full         <<endl;
    cout<<"<<<<< max_count_full        : "<<right<<setw(10)<<max_count_full     <<endl;
    cout<<"<<<<< max_count_relax       : "<<right<<setw(10)<<max_count_relax    <<endl;
    cout<<"<<<<< fill_coeff            : "<<right<<setw(10)<<fill_coeff         <<endl;
    cout<<"**********************************************************\n";
  }  
  max_scale = level_set.maximum_dl*CELL_RATIO;
  min_scale = level_set.minimum_dl*CELL_RATIO;
  ini_scale = max_scale;

  // Define_level_infor (world);

  Define_level_infor_tight (world);
  
  Allocate_particles (world);
  
  #ifdef _MPI_
    sph_voronoi.Initialize(this, 1, world);

    #if defined (_CVP_LSET_)

      std::vector<my_real> vp_coords; vp_coords.clear();
      std::vector<Real>    vp_scale ; vp_scale.clear();

      if (world.rank() == 0)
        level_set.CVP_for_initial_partitioning (vp_coords, vp_scale, world);

      sph_voronoi.Set_VP_coords_and_hmin(vp_coords, vp_scale, world);

      Random_particle_distribution_within_VPi (world);

    #elif defined (_CLOSE_PACKING_)

  //     sph_voronoi.Random_particle_distribution_inbox(CVP_box_l, CVP_box, world);
      Real R = sph_voronoi.Hex_close_packing_inbox(CVP_box_l, CVP_box, dimX, dimY, dimZ, world);

      Random_particle_distribution_within_R(R, world);

    #elif defined (_READ_VP_)
      
      std::vector<my_real> vp_coords; vp_coords.clear();
      std::vector<Real>    vp_scale ; vp_scale.clear();

      if (world.rank() == 0){

        char   vp_file[256];
        char   buffer[500];
        sprintf(vp_file,"%s","./rstfile/VP_coords.csv");

        ifstream load(vp_file);
        
        if (! load.is_open())  { 
          cout << "<<<<< Error opening VP_coords.csv in rank "<<world.rank(); 
          world.abort (-1); 
        }

        int num_color_from_file = 0;
        while (load.getline (buffer,500))
          num_color_from_file++;
        num_color_from_file--;

        if (num_color_from_file != world.size()){
          cout<<"<<<<< ERROR num_color_from_file"<<endl;
          world.abort(-1);
        }

        load.clear();
        load.seekg(0, ios::beg);

        load.getline (buffer,500);
        num_color_from_file = 0;

        for (int i = 0; i < world.size(); ++i)
        {
          int tmp = 0;
          my_real coord_;
          Real    h_min_;
          char    dummy;
          load>>tmp>>dummy>>coord_.i>>dummy>>coord_.j>>dummy>>coord_.k>>dummy>>h_min_;
          // cout<<"<<< "<<tmp<<" "<<dummy<<" "<<coord_.i<<" "<<dummy<<" "<<coord_.j<<" "<<dummy<<" "<<coord_.k<<" "<<dummy<<" "<<h_min_<<endl;
          vp_coords.push_back(coord_);
          vp_scale.push_back(h_min_);
        }
      }

      sph_voronoi.Set_VP_coords_and_hmin(vp_coords, vp_scale, world);

      Random_particle_distribution_within_VPi (world);

      // Random_particle_distribution (world);

    #endif
  #else
    Random_particle_distribution(world);
  #endif

  world.barrier();

  // Output_paraview_file(0, 1, world);
  
  Initialize_system(world);

  // Output_paraview_file(1, 1, world);
  
  Asign_local_particle_index (world);
  
  // sph_voronoi.Output_vtk_file(0, 0, world);
  // world.barrier();
  // world.abort(-1);

  Pre_run (world);

  Output_simulation_infor (0, world);
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<class sph initialized\n";
  out.close();
#endif
}
//--------------------------------------------------
// Get_num_particles
//--------------------------------------------------
void SPH_mesh_generation::Get_num_particles (communicator &world)
{
  // total_num_particle_singularity = level_set.num_singularity;
  // total_num_particle_segment     = int(level_set.total_mass_segment*1.) - 0.5*total_num_particle_singularity;
  // total_num_particle_surface     = int(level_set.total_mass_surface*1.) - 0.5*total_num_particle_segment;
  // // Approx. number of mesh points inside the volume
  // total_num_particle_mesh        = int(level_set.total_mass*1.)         - 0.5*total_num_particle_surface;

  // 20181217::Approx. number of mesh points inside the volume
  total_num_particle_singularity = level_set.num_singularity;
  total_num_particle_segment     = AMAX1(int(level_set.total_mass_segment*1. - 0.5*level_set.num_singularity), 0);
  total_num_particle_surface     = AMAX1(int(level_set.total_mass_surface*1. - 0.5*level_set.total_mass_segment*1.), 0);
  total_num_particle_mesh        = AMAX1(int(level_set.total_mass*1.         - 0.5*level_set.total_mass_surface*1.), 0);

  int n_interface_cell = level_set.glbl_num_interface_cell;
  
  if (n_interface_cell < total_num_particle_surface){
    // TODO find a solution to this
    cout<<"<<<<< WARNING!!: Scales are not set correct for surface particles. n_interface_cell: "<<n_interface_cell<<" total_num_particle_surface: "<<glbl_total_num_particle_surface<<endl;
    cout<<"<<<<< total_num_particle_surface is set to "<<n_interface_cell<<endl;
    total_num_particle_surface = n_interface_cell;
  }
  
  total_num_particle = total_num_particle_mesh + total_num_particle_surface + total_num_particle_singularity + total_num_particle_segment;

  glbl_total_num_particle_mesh        = total_num_particle_mesh;
  glbl_total_num_particle_surface     = total_num_particle_surface;
  glbl_total_num_particle_segment     = total_num_particle_segment;
  glbl_total_num_particle_singularity = total_num_particle_singularity;
  
  glbl_total_num_particle = glbl_total_num_particle_mesh + glbl_total_num_particle_surface + glbl_total_num_particle_segment + glbl_total_num_particle_singularity;

  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< The total mesh particle is           : "<<right<<setw(10)<<glbl_total_num_particle_mesh       <<endl;
    cout<<"<<<<< The total surface particle is        : "<<right<<setw(10)<<glbl_total_num_particle_surface    <<endl;
    cout<<"<<<<< The total segment particle is        : "<<right<<setw(10)<<glbl_total_num_particle_segment    <<endl;
    cout<<"<<<<< The total characteristic particle is : "<<right<<setw(10)<<glbl_total_num_particle_singularity<<endl;
    cout<<"<<<<< The total number particle is         : "<<right<<setw(10)<<glbl_total_num_particle            <<endl;
    cout<<"**********************************************************\n";
  }
}
//--------------------------------------------------
// Recalculate_local_particle_number
//--------------------------------------------------
void SPH_mesh_generation::Recalculate_local_particle_number (communicator &world)
{
  int tmp = int(floor(glbl_total_num_particle / world.size()));
  int rest = glbl_total_num_particle - tmp*world.size();

  for (int i = 0; i < rest; ++i)
  {
    if (world.rank() == i) tmp ++;
  }

  int glbl_tmp = 0;

  all_reduce( world, tmp, glbl_tmp, std::plus<int>());

  if (glbl_tmp != glbl_total_num_particle){
    if (world.rank() == 0) cout<<"ERROR in glbl_total_num_particle!!!"<<endl;
  }else{
    total_num_particle = tmp;
    
    total_num_particle_mesh = 0;
    total_num_particle_surface = 0;
    total_num_particle_segment = 0;
    total_num_particle_singularity = 0;

    all_gather( world, total_num_particle, nparticle_in_each_cpu);
  }

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Recalculate_local_particle_number finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Pre_run
//--------------------------------------------------
void SPH_mesh_generation::Set_timestep (communicator &world, int stage)
{
  local_timestep = 1.e20;
  max_v    = -1.e20;

  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      Real temp;
      Particle *current = particle[i];
      current->Set_timestep(stage, this);
    }
  }, ap);
  ParallelGetTimestep(total_num_particle,particle,local_timestep, max_v);
  all_reduce(world, local_timestep, glbl_timestep, mpi::minimum<Real>());
  all_reduce(world, max_v, glbl_max_v, mpi::maximum<Real>());
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Set_timestep finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Pre_run
//--------------------------------------------------
void SPH_mesh_generation::Pre_run (communicator &world)
{  

  Get_particle_info (world);

  Set_particle_scale (world);

  Reset_particle_array (world);
  
  Check_total_number_of_particle (world);
  
  Find_and_map_particle (world);
  
  Reset_particle_force (world, 1);
  
  Prepare_for_the_current_step (world);
  
  Accumulate_particle_force (world, 1);
    
  Set_timestep(world, 1);
  
  Reset_particle_force (world, 1);
   
  Accumulate_particle_force (world, 2);
  
  Set_timestep(world, 2);
 
  Reset_particle_force (world, 2);

  world.barrier();
//       run_time ++;
//   Output_paraview_file(0, 1, world);
//       world.barrier();
//   world.abort(-1);
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Pre_run finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Allocate memory for SPH particles
//--------------------------------------------------
void SPH_mesh_generation::Allocate_particles (communicator &world)
{
  particle = new p_Particle[total_num_particle];
  
  int pid_start = 0;
  for (int irank = 0; irank < world.rank(); irank++){
    pid_start += nparticle_in_each_cpu[irank];
  }

  int pid = 0;
  // JZ20181220 :: set all particles as Real (inside) particles
  for(int i = 0; i < total_num_particle; i++){
    particle[pid] = particlepool.malloc();
    real_particle.push_back(particle[pid]);
    particle[pid]->Initialize (this, pid+pid_start, REAL_PARTICLE, world.rank());
    pid++;
  }

  // JZ20181220 this is only an initial implementation
  // for(int i = 0; i < total_num_particle_mesh; i++){
  //   particle[pid] = particlepool.malloc();
  //   real_particle.push_back(particle[pid]);
  //   particle[pid]->Initialize (this, pid, REAL_PARTICLE);
  //   pid++;
  // }
  // for(int i = 0; i < total_num_particle_surface; i++){
  //   particle[pid] = particlepool.malloc();
  //   real_particle.push_back(particle[pid]);
  //   particle[pid]->Initialize (this, pid, REAL_PARTICLE);
  //   pid++;
  // }
  // for(int i = 0; i < total_num_particle_segment; i++){
  //   particle[pid] = particlepool.malloc();
  //   real_particle.push_back(particle[pid]);
  //   particle[pid]->Initialize (this, pid, REAL_PARTICLE);
  //   pid++;
  // }
  // for(int i = 0; i < total_num_particle_singularity; i++){
  //   particle[pid] = particlepool.malloc();
  //   real_particle.push_back(particle[pid]);
  //   particle[pid]->Initialize (this, pid, REAL_PARTICLE);
  //   pid++;
  // }
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Allocate_particles finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Asign_local_particle_index
//--------------------------------------------------
void SPH_mesh_generation::Asign_local_particle_index (communicator &world)
{
  local_id_record = 0;
  
  int pid = 0;
  for(int i = 0; i < total_num_particle; i++){
    particle[i]->Set_local_index (i);
    local_id_record ++;
  }
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Asign_local_particle_index finished\n";
  out.close();
#endif
}
#ifdef _MPI_
  #if defined (_CVP_LSET_) || defined (_READ_VP_)
//--------------------------------------------------
// Random_particle_distribution_within_VPi
//--------------------------------------------------
void SPH_mesh_generation::Random_particle_distribution_within_VPi (communicator &world)
{
  my_real center = sph_voronoi.Get_coord_of_VP_at_iRank(world.rank());
  Real         R = sph_voronoi.Get_h_min_of_VP_at_iRank(world.rank());

  if (R < 1.e-6){
    cout<<"<<< ["<<world.rank()<<"]-> ERROR in h_min:"<<R<<endl;
    world.abort(-1);
  }
  cout<<"<<< ["<<world.rank()<<"]->R:"<<R<<endl;

  my_real _box_l;
  my_real _box_r;
  my_real _box;
  
  _box_l.i = DIM_X ? center.i - R : 0.;
  _box_l.j = DIM_Y ? center.j - R : 0.;
  _box_l.k = DIM_Z ? center.k - R : 0.;

  _box_r.i = DIM_X ? center.i + R : 0.;
  _box_r.j = DIM_Y ? center.j + R : 0.;
  _box_r.k = DIM_Z ? center.k + R : 0.;

  my_real lset_box_l = level_set.box_l;
  my_real lset_box_r = level_set.box_r;
  Real padding = level_set.lset_level_info[0]->scale;

  my_set_const (_box, 2*R);
  
  srand((unsigned)time(NULL)); 
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      my_real position;

      bool positive_phase = false;
      bool in_R = false;

      while (!positive_phase || !in_R){
        position.i = 0.;
        position.j = 0.;
        position.k = 0.;
        positive_phase = false;
        in_R = false;

    #if (DIM_X)
        position.i = rand()/ double(RAND_MAX+1.0);
        particle[i]->coord.i = position.i*_box.i + rand()/ double(RAND_MAX+1.0)*1.e-3 + _box_l.i;
        particle[i]->coord.i = AMAX1(_box_l.i+1.e-10,AMIN1(particle[i]->coord.i, _box_r.i-1.e-10));
        particle[i]->coord.i = AMAX1(lset_box_l.i+padding,AMIN1(lset_box_r.i-padding,particle[i]->coord.i));
    #endif
    #if (DIM_Y)
        position.j = rand()/ double(RAND_MAX+1.0);
        particle[i]->coord.j = position.j*_box.j + rand()/ double(RAND_MAX+1.0)*1.e-3 + _box_l.j;
        particle[i]->coord.j = AMAX1(_box_l.j+1.e-10,AMIN1(particle[i]->coord.j, _box_r.j-1.e-10)); 
        particle[i]->coord.j = AMAX1(lset_box_l.j+padding,AMIN1(lset_box_r.j-padding,particle[i]->coord.j));
    #endif
    #if (DIM_Z)
        position.k = rand()/ double(RAND_MAX+1.0);
        particle[i]->coord.k = position.k*_box.k + rand()/ double(RAND_MAX+1.0)*1.e-3 + _box_l.k;
        particle[i]->coord.k = AMAX1(_box_l.k+1.e-10,AMIN1(particle[i]->coord.k, _box_r.k-1.e-10));
        particle[i]->coord.k = AMAX1(lset_box_l.k+padding,AMIN1(lset_box_r.k-padding,particle[i]->coord.k));
    #endif
        particle[i]->Calculate_particle_infor(this);
        particle[i]->Set_particle_info(this);

        // if (particle[i]->phi > level_set.maximum_dl) positive_phase = true;
        if (particle[i]->phi > particle[i]->h) positive_phase = true;
        
        Real dist = get_distance_2p (particle[i]->coord, center);

        if (dist <= R) in_R = true;
      }
    }
  }, ap);
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Random_particle_distribution_within_R finished\n";
  out.close();
#endif
}
  #endif
//--------------------------------------------------
// Random_particle_distribution_within_R
//--------------------------------------------------
void SPH_mesh_generation::Random_particle_distribution_within_R(Real R, communicator &world)
{
  my_real center = sph_voronoi.Get_coord_of_VP_at_iRank(world.rank());
  my_real _box_l;
  my_real _box_r;
  my_real _box;
  
  _box_l.i = DIM_X ? center.i - R : 0.;
  _box_l.j = DIM_Y ? center.j - R : 0.;
  _box_l.k = DIM_Z ? center.k - R : 0.;

  _box_r.i = DIM_X ? center.i + R : 0.;
  _box_r.j = DIM_Y ? center.j + R : 0.;
  _box_r.k = DIM_Z ? center.k + R : 0.;

  my_set_const (_box, 2*R);

  my_real lset_box_l = level_set.box_l;
  my_real lset_box_r = level_set.box_r;
  Real padding = level_set.lset_level_info[0]->scale;
  
  srand((unsigned)time(NULL)); 
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      my_real position;

      bool positive_phase = false;
      bool in_R = false;

      while (!positive_phase || !in_R){
        position.i = 0.;
        position.j = 0.;
        position.k = 0.;
        positive_phase = false;
        in_R = false;

    #if (DIM_X)
        position.i = rand()/ double(RAND_MAX+1.0);
        particle[i]->coord.i = position.i*_box.i + rand()/ double(RAND_MAX+1.0)*1.e-3 + _box_l.i;
        particle[i]->coord.i = AMAX1(_box_l.i+1.e-10,AMIN1(particle[i]->coord.i, _box_r.i-1.e-10));
        particle[i]->coord.i = AMAX1(lset_box_l.i+padding,AMIN1(lset_box_r.i-padding,particle[i]->coord.i));
    #endif
    #if (DIM_Y)
        position.j = rand()/ double(RAND_MAX+1.0);
        particle[i]->coord.j = position.j*_box.j + rand()/ double(RAND_MAX+1.0)*1.e-3 + _box_l.j;
        particle[i]->coord.j = AMAX1(_box_l.j+1.e-10,AMIN1(particle[i]->coord.j, _box_r.j-1.e-10)); 
        particle[i]->coord.j = AMAX1(lset_box_l.j+padding,AMIN1(lset_box_r.j-padding,particle[i]->coord.j));
    #endif
    #if (DIM_Z)
        position.k = rand()/ double(RAND_MAX+1.0);
        particle[i]->coord.k = position.k*_box.k + rand()/ double(RAND_MAX+1.0)*1.e-3 + _box_l.k;
        particle[i]->coord.k = AMAX1(_box_l.k+1.e-10,AMIN1(particle[i]->coord.k, _box_r.k-1.e-10));
        particle[i]->coord.k = AMAX1(lset_box_l.k+padding,AMIN1(lset_box_r.k-padding,particle[i]->coord.k));
    #endif
        particle[i]->Calculate_particle_infor(this);
        particle[i]->Set_particle_info(this);

        // if (particle[i]->phi > level_set.maximum_dl) positive_phase = true;
        if (particle[i]->phi > particle[i]->h) positive_phase = true;
        
        Real dist = get_distance_2p (particle[i]->coord, center);

        if (dist <= R) in_R = true;
      }
    }
  }, ap);
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Random_particle_distribution_within_R finished\n";
  out.close();
#endif
}
#endif
//--------------------------------------------------
// Random_particle_distribution
//--------------------------------------------------
void SPH_mesh_generation::Random_particle_distribution(communicator &world)
{

  my_real lset_box_l = level_set.box_l;
  my_real lset_box_r = level_set.box_r;
  Real padding = level_set.lset_level_info[0]->scale;

  srand((unsigned)time(NULL)); 
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      my_real position;
      // if (particle[i]->type == SINGULARITY_PARTICLE) continue;

      bool positive_phase = false;
      bool in_the_current_Voronoi = false;

      while (!positive_phase || !in_the_current_Voronoi){
        position.i = 0.;
        position.j = 0.;
        position.k = 0.;
        positive_phase = false;
        in_the_current_Voronoi = false;

    #if (DIM_X)
        position.i = rand()/ double(RAND_MAX+1.0);
        particle[i]->coord.i = position.i*particle_box.i + rand()/ double(RAND_MAX+1.0)*1.e-3 + particle_box_l.i;
        particle[i]->coord.i = AMAX1(particle_box_l.i+1.e-10,AMIN1(particle[i]->coord.i, particle_box_r.i-1.e-10));
        particle[i]->coord.i = AMAX1(lset_box_l.i+padding,AMIN1(lset_box_r.i-padding,particle[i]->coord.i));
    #endif
    #if (DIM_Y)
        position.j = rand()/ double(RAND_MAX+1.0);
        particle[i]->coord.j = position.j*particle_box.j + rand()/ double(RAND_MAX+1.0)*1.e-3 + particle_box_l.j;
        particle[i]->coord.j = AMAX1(particle_box_l.j+1.e-10,AMIN1(particle[i]->coord.j, particle_box_r.j-1.e-10)); 
        particle[i]->coord.j = AMAX1(lset_box_l.j+padding,AMIN1(lset_box_r.j-padding,particle[i]->coord.j));
    #endif
    #if (DIM_Z)
        position.k = rand()/ double(RAND_MAX+1.0);
        particle[i]->coord.k = position.k*particle_box.k + rand()/ double(RAND_MAX+1.0)*1.e-3 + particle_box_l.k;
        particle[i]->coord.k = AMAX1(particle_box_l.k+1.e-10,AMIN1(particle[i]->coord.k, particle_box_r.k-1.e-10));
        particle[i]->coord.k = AMAX1(lset_box_l.k+padding,AMIN1(lset_box_r.k-padding,particle[i]->coord.k));
    #endif
        particle[i]->Calculate_particle_infor(this);
        particle[i]->Set_particle_info(this);

        // if (particle[i]->phi > level_set.maximum_dl) positive_phase = true;
        if (particle[i]->phi > particle[i]->h) positive_phase = true;

        if      (world.size() == 1) in_the_current_Voronoi = true;
    #ifdef _MPI_
        else if (sph_voronoi.Get_CPU_id(particle[i]->coord) == world.rank()) {
          in_the_current_Voronoi = true;
        }
    #endif
      }
    }
  }, ap);
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Random_particle_distribution finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// SPH_run_mesh_generation
//--------------------------------------------------
void SPH_mesh_generation::SPH_run_mesh_generation(communicator &world)
{
    time_simulation.restart();
    time_update.restart();
    
  Update_coord_full_step (world);
  
    time_for_update[world.rank()] += time_update.elapsed();
    time_rest.restart();
    
  Get_particle_info (world);
  
    time_for_rest[world.rank()] += time_rest.elapsed();
    time_mapping.restart();
    
  Find_and_map_particle (world);
  
    time_for_mapping[world.rank()] += time_mapping.elapsed();
    time_rest.restart();
    
  Set_particle_scale (world);
  
    time_for_rest[world.rank()] += time_rest.elapsed();
    time_for_simulation[world.rank()] += time_simulation.elapsed();
      
  Prepare_for_the_current_step (world);
  
    time_simulation.restart();
    time_rest.restart();
    
  Reset_particle_force (world, 1);
  
    time_for_rest[world.rank()] += time_rest.elapsed();  
    time_density.restart();
    
  Get_kernel_summation (world);
  
    time_for_density[world.rank()] += time_density.elapsed();
    time_for_simulation[world.rank()] += time_simulation.elapsed();

  Refresh_particle_infor(world, 1, 1);

    time_simulation.restart();
    time_force_calculation.restart();
    
  Accumulate_particle_force (world, 1);
  
    time_for_force_calculation[world.rank()] += time_force_calculation.elapsed();
    time_rest.restart();  
    
  Set_timestep(world, 1);
  
    time_for_rest[world.rank()] += time_rest.elapsed();  
    time_update.restart();
    
  Update_velocity_half_step (world);
  
    time_for_update[world.rank()] += time_update.elapsed();
    time_for_simulation[world.rank()] += time_simulation.elapsed();

  Refresh_particle_infor(world, 2, 2);

    time_simulation.restart(); 
    time_rest.restart();  
    
  Reset_particle_force (world, 2);
  
    time_for_rest[world.rank()] += time_rest.elapsed();  
    time_force_calculation.restart();

  Accumulate_particle_force (world, 2);

    time_for_force_calculation[world.rank()] += time_force_calculation.elapsed();
    time_rest.restart();

  Set_timestep(world, 2);

    time_for_rest[world.rank()] += time_rest.elapsed();  
    time_update.restart();
  
  Update_velocity_half_step (world);
  
    time_for_update[world.rank()] += time_update.elapsed();
    time_for_simulation[world.rank()] += time_simulation.elapsed();

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<SPH_run_mesh_generation finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Get_particle_info
//--------------------------------------------------
void SPH_mesh_generation::Get_particle_info(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle[i]->Calculate_particle_infor(this);
    }
  }, ap);
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Get_particle_info finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Set_particle_scale
//--------------------------------------------------
void SPH_mesh_generation::Set_particle_scale(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle[i]->Set_particle_info(this);
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Set_particle_scale finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Find_and_map_particle
//--------------------------------------------------
void SPH_mesh_generation::Find_and_map_particle (communicator &world)
{
  if (current_stage == FILL){
    int num_interface_cell = level_set.glbl_num_interface_cell;
    find_surface_particle = false;
    find_singularity_particle = false;
    find_segment_particle = false;

    static affinity_partitioner ap;
    parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        int type_ = particle[i]->type;
        
        particle[i]->Find_and_map_particle (this);

        if (particle[i]->type == SEGMENT_PARTICLE || particle[i]->type == SURFACE_PARTICLE){
          // small hack here to avoid issues in initial mesh generation
          particle[i]->Get_dispersed_position(0.00001*particle[i]->h);
        }

        if (type_ == REAL_PARTICLE && particle[i]->type == SURFACE_PARTICLE){
          find_surface_particle = true;
        }
        if (type_ != SEGMENT_PARTICLE && particle[i]->type == SEGMENT_PARTICLE){
          find_segment_particle = true;
        }
        if (type_ != SINGULARITY_PARTICLE && particle[i]->type == SINGULARITY_PARTICLE){
          find_singularity_particle = true;
        }

        particle[i]->Calculate_particle_infor (this);

        if (particle[i]->phi < 0. && particle[i]->type == REAL_PARTICLE){

          Remap_this_particle_to_positive_phase (particle[i]);
        
        }
      }
    }, ap);
  }
  // JZ20181129::seems that the relax full stage is not needed anymore, just keep it for future usage
  // else if (current_stage == FULL || (current_stage == RELAX && count_full < max_count_relax - 50 /*int((max_count_relax-100)*0.75)*/)){
  //   int num_interface_cell = level_set.glbl_num_interface_cell;
  //   find_singularity_particle = false;
  
  //   int current_local_num_char_particle = particle_characteristic.size();
  //   current_glbl_total_num_particle_singularity = 0;
  //   reduce(world, current_local_num_char_particle, current_glbl_total_num_particle_singularity, std::plus<int>(), 0);
    
  //   if (current_glbl_total_num_particle_singularity < glbl_total_num_particle_singularity){
  //     int current_num_surf_particle = particle_surface.size();
  //     static affinity_partitioner ap;
  //     parallel_for( blocked_range<int>(0, current_num_surf_particle ),[&](const blocked_range<int>& r){
  //       for(int i=r.begin(); i!=r.end(); ++i){
  //         int type_ = particle_surface[i]->type;
  //         particle_surface[i]->Find_characteristic_particle (this);

  //         if (type_ == SURFACE_PARTICLE && particle_surface[i]->type == SINGULARITY_PARTICLE){
  //           find_singularity_particle = true;
  //           particle_characteristic.push_back(particle_surface[i]);

  //           particle_surface[i]->Calculate_particle_infor(this);
  //           particle_surface[i]->Set_particle_info(this);
  //         }
  //       }
  //     }, ap);
    
  //     if (find_singularity_particle){
  //       current_local_num_char_particle = particle_characteristic.size();
  //       current_glbl_total_num_particle_singularity = 0;
  //       reduce(world, current_local_num_char_particle, current_glbl_total_num_particle_singularity, std::plus<int>(), 0);
    
  //       concurrent_vector<std::pair<int, p_Particle>> index_cha; index_cha.clear();
  //       parallel_for( blocked_range<int>(0, current_local_num_char_particle ),[&](const blocked_range<int>& r){
  //         for(int i=r.begin(); i!=r.end(); ++i){
  //           int tag_interface = NORMAL_CELL;
  //           int tag_characteristic = NORMAL_CELL;
  
  //           particle_characteristic[i]->Get_current_particle_level_set_cell_tag (this, tag_interface, tag_characteristic);
  
  //           my_int  pos_pkg;
  //           my_int  pos_cell;
        
  //           my_real coord_shift = my_minus_data (particle_characteristic[i]->coord, box_l);
        
  //           level_set.lset_level_info[0]->Get_id_pkg_cell (coord_shift, pos_pkg, pos_cell);
        
  //           int unique_index = level_set.lset_level_info[0]->Get_unique_cell_id (pos_pkg, pos_cell);
        
  //           std::pair<int, p_Particle> pair_tmp;
  //           pair_tmp.first = unique_index;
  //           pair_tmp.second = particle_characteristic[i];
        
  //           index_cha.push_back(pair_tmp);
            
  //           if (tag_characteristic != SINGULARITY_CELL){
  //             cout<<"ERROR!!!: not a singularity cell\n";
  //             world.abort(-1);
  //           }
  //         }
  //       }, ap);
    
  //       sort_by_key (index_cha, current_local_num_char_particle);
  
  //       for (int i = 0; i < current_local_num_char_particle-1; i++){
  //         if (index_cha[i].first == index_cha[i+1].first){
  //           index_cha[i].second->type = SURFACE_PARTICLE;
  //           particle_surface.push_back(index_cha[i].second);

  //           index_cha[i].second->Calculate_particle_infor(this);
  //           index_cha[i].second->Set_particle_info(this);
            
  //           index_cha[i+1].second->Map_characteristic_particle(this);
  //           index_cha[i+1].second->Calculate_particle_infor(this);
  //           index_cha[i+1].second->Set_particle_info(this);
  //         }
  //       }
  //     }
  //   }
  // }

  int do_reset = 0;
  int glbl_do_reset = 0;
  if ((find_surface_particle || find_segment_particle || find_singularity_particle)){    
    do_reset = 1;
  }
  
  all_reduce(world, do_reset, glbl_do_reset, mpi::maximum<int>());
  
  if (glbl_do_reset){
    
    Reset_particle_array (world);
    
    Check_total_number_of_particle (world);
    
    find_surface_particle = false;
    find_segment_particle = false;
    find_singularity_particle = false;
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Find_and_map_particle finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Check and pop surface particle 
// when simulation becomes stable
//--------------------------------------------------
void SPH_mesh_generation::Check_and_pop_surface_particle (communicator &world)
{
  if (int(particle_surface.size()) < int (glbl_total_num_particle_surface)){
    // allow to fill
    count_full = 0;
    current_stage = FILL;
  }
  if (current_stage == FILL){
    if (int(particle_surface.size()) > int (glbl_total_num_particle_surface)){
      count_full ++;
    }
    if ((count_full >= max_count_full && int(particle_surface.size()) > int (glbl_total_num_particle_surface)) || int(particle_surface.size()) > int (glbl_total_num_particle_surface*fill_coeff)){
    // stop fill
      current_stage = RELAX;
      count_full = max_count_relax;
          
      int current_local_num_surf_particle = particle_surface.size();
      current_glbl_total_num_particle_surface_old = 0;
      reduce(world, current_local_num_surf_particle, current_glbl_total_num_particle_surface_old, std::plus<int>(), 0);
    }
  }

  // JZ20181129::seems that the relax full stage is not needed anymore, just keep it for future usage
  // if (current_stage == RELAX){
  //   //start to pop
  //   int num_particle_to_pop = particle_surface.size() - glbl_total_num_particle_surface;
  //   count_full --;
    
  //   if ((particle_surface.size()-num_particle_to_pop >= glbl_total_num_particle_surface)){
  //     // ramdomly pop particle out
  //     srand((unsigned)time(NULL));
  //     std::vector<int> id_pop; id_pop.clear();
  //     for (int i = 0; i < num_particle_to_pop; i++){
  //       int idd = AMIN1(int(particle_surface.size()-1), AMAX1(0, int(rand()/ double(RAND_MAX+1.0)*particle_surface.size())));
  //       id_pop.push_back(idd);
  //     }
  //     // check if two numbers are the same
  //     bool rrun = true;
  //     while(rrun){
  //       rrun = false;
  //       std::sort( id_pop.begin(), id_pop.end());
  //       for (int i = 1; i < num_particle_to_pop-1; i++){
  //         if (id_pop[i] == id_pop[i-1] || id_pop[i] == id_pop[i+1] || particle_surface[id_pop[i]]->h > level_set.middel_dl){
  //           int idd = AMIN1(int(particle_surface.size()-1), AMAX1(0, int(rand()/ double(RAND_MAX+1.0)*particle_surface.size())));
  //           id_pop[i] = idd;
  //           rrun = true;
  //         }
  //       }
  //     }
  //     // now pop particles
  //     static affinity_partitioner ap;
  //     parallel_for( blocked_range<int>(0, num_particle_to_pop),[&](const blocked_range<int>& r){
  //       for(int i=r.begin(); i!=r.end(); ++i){
  //         p_Particle current_particle = particle_surface[id_pop[i]];
          
  //         Remap_this_particle_to_positive_phase (current_particle );
  
  //         real_particle.push_back(current_particle);
  //       }
  //     }, ap);
      
  //     particle_surface.clear();
  //     parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
  //       for(int i=r.begin(); i!=r.end(); ++i){
  //         if (particle[i]->type == SURFACE_PARTICLE){
  //           particle_surface.push_back(particle[i]);
  //         }
  //       }
  //     }, ap);

  //     Check_total_number_of_particle (world);
      
  //     int current_local_num_surf_particle = particle_surface.size();
  //     current_glbl_total_num_particle_surface = 0;
  //     reduce(world, current_local_num_surf_particle, current_glbl_total_num_particle_surface, std::plus<int>(), 0);
      
  //     if (current_glbl_total_num_particle_surface != glbl_total_num_particle_surface){
  //       cout<<"<<<<<ERROR: current_glbl_total_num_particle_surface: "<<current_glbl_total_num_particle_surface<<"  glbl_total_num_particle_surface: "<<glbl_total_num_particle_surface<<endl;
  //     }
  //     if (count_full % 10 == 0)
  //       Output_simulation_infor(1, world);
  //   }
  //   if (count_full == 0){
  //     current_stage = FULL;
  //     max_count_relax = AMAX1(100, max_count_relax/2);
  //     fill_coeff = fill_coeff - (fill_coeff-1.)/2;
  //     Output_simulation_infor(1, world);
  //   }
  // }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Check_and_pop_surface_particle finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Check_total_number_of_particle
//--------------------------------------------------
void SPH_mesh_generation::Check_total_number_of_particle (communicator &world)
{
  int current_total_num_particle = real_particle.size() + particle_surface.size() + particle_singularity.size() + particle_segment.size();
  int current_glbl_total_num_particle= 0;

  all_reduce( world, current_total_num_particle, current_glbl_total_num_particle, std::plus<int>());

  if (current_glbl_total_num_particle != glbl_total_num_particle){

    if (world.rank() == 0){
      cout<<"<<<<< ERROR!! wrong total number of particles: previous "<<glbl_total_num_particle<<" | now "<<current_glbl_total_num_particle<<endl;
      cout<<"<<<<<The total mesh particle is        : "<<real_particle.size()       <<endl;
      cout<<"<<<<<The total surface particle is     : "<<particle_surface.size()    <<endl;
      cout<<"<<<<<The total segment particle is     : "<<particle_segment.size()    <<endl;
      cout<<"<<<<<The total singularity particle is : "<<particle_singularity.size()<<endl;      
    }

  }else{    
    total_num_particle_mesh        = real_particle.size();
    total_num_particle_surface     = particle_surface.size();
    total_num_particle_segment     = particle_singularity.size();
    total_num_particle_singularity = particle_segment.size();

    current_glbl_total_num_particle_mesh        = 0;
    current_glbl_total_num_particle_surface     = 0;
    current_glbl_total_num_particle_segment     = 0;
    current_glbl_total_num_particle_singularity = 0;

    reduce(world, total_num_particle_mesh       , current_glbl_total_num_particle_mesh       , std::plus<int>(), 0);
    reduce(world, total_num_particle_surface    , current_glbl_total_num_particle_surface    , std::plus<int>(), 0);
    reduce(world, total_num_particle_segment    , current_glbl_total_num_particle_segment    , std::plus<int>(), 0);
    reduce(world, total_num_particle_singularity, current_glbl_total_num_particle_singularity, std::plus<int>(), 0);
  }
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Check_total_number_of_particle finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Remap_particle
//--------------------------------------------------
void SPH_mesh_generation::Remap_particle (communicator &world)
{
  srand((unsigned)time(NULL));
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle[i]->Calculate_particle_infor (this);
      if (particle[i]->type == SURFACE_PARTICLE || particle[i]->type == SEGMENT_PARTICLE){
        
        int tag_interface = NORMAL_CELL;
        int tag_characteristic = NORMAL_CELL;
        int idx_characteristic = -1;
        
        particle[i]->Get_level_set_char_cell_tag_and_idx (this, tag_interface, tag_characteristic, idx_characteristic);
        // particle[i]->Get_extended_char_cell_tag_and_idx (this, tag_interface, tag_characteristic, idx_characteristic);
        
        if (tag_interface != CUT_CELL && tag_characteristic == NORMAL_CELL){
          particle[i]->type = REAL_PARTICLE;
        }
      }
      
      if (particle[i]->phi < 0. && particle[i]->type == REAL_PARTICLE){
        
        Remap_this_particle_to_positive_phase (particle[i]);

      }
    }
  }, ap);
  
  Reset_particle_array (world);

  Check_total_number_of_particle (world);
  
  if (int(particle_surface.size()) < int (glbl_total_num_particle_surface)){
    // allow to fill
    count_full = 0;
    current_stage = FILL;
  }
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Remap_negtive_particle finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Reset_particle_array
//--------------------------------------------------
void SPH_mesh_generation::Reset_particle_array (communicator &world)
{
  real_particle.clear();
  particle_surface.clear();
  particle_segment.clear();
  particle_singularity.clear();
  
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      if (particle[i]->type == REAL_PARTICLE){
        real_particle.push_back(particle[i]);
      }else if (particle[i]->type == SURFACE_PARTICLE){
        particle_surface.push_back(particle[i]);
      }else if (particle[i]->type == SEGMENT_PARTICLE){
        particle_segment.push_back(particle[i]);
      }else{
        particle_singularity.push_back(particle[i]);
      }
    }
  }, ap);
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Reset_particle_array finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Remap_this_particle_to_positive_phase
//--------------------------------------------------
void SPH_mesh_generation::Remap_this_particle_to_positive_phase (p_Particle particle)
{
  my_real position;
 
  bool positive_phase = false;
  
  my_real box_;
  my_real box_l_;
  my_real box_r_;
  
  #ifdef _MPI_
  my_real local_domain = my_minus_data (local_box_r, local_box_l);
  my_real local_center = my_add_data (local_box_l, my_multiply_const (local_domain, 0.5));
  my_set_data (box_  , my_multiply_const (local_domain, 0.5));
  my_set_data (box_l_, my_minus_data     (local_center, my_multiply_const (box_, 0.5)));
  my_set_data (box_r_, my_add_data       (local_center, my_multiply_const (box_, 0.5)));
  #else
  my_set_data (box_  , particle_box);
  my_set_data (box_l_, particle_box_l);
  my_set_data (box_r_, particle_box_r);
  #endif
  
  while (!positive_phase ){
    #if (DIM_X)
    position.i = rand()/ double(RAND_MAX+1.0);
    particle->coord.i = position.i*box_.i + rand()/ double(RAND_MAX+1.0)*1.e-3 + box_l_.i;
    particle->coord.i = AMAX1(box_l_.i+1.e-10,AMIN1(particle->coord.i, box_r_.i-1.e-10));
    #endif
    #if (DIM_Y)
    position.j = rand()/ double(RAND_MAX+1.0);
    particle->coord.j = position.j*box_.j + rand()/ double(RAND_MAX+1.0)*1.e-3 + box_l_.j;
    particle->coord.j = AMAX1(box_l_.j+1.e-10,AMIN1(particle->coord.j, box_r_.j-1.e-10)); 
    #endif
    #if (DIM_Z)
    position.k = rand()/ double(RAND_MAX+1.0);
    particle->coord.k = position.k*box_.k + rand()/ double(RAND_MAX+1.0)*1.e-3 + box_l_.k;
    particle->coord.k = AMAX1(box_l_.k+1.e-10,AMIN1(particle->coord.k, box_r_.k-1.e-10));
    #endif
    particle->Calculate_particle_infor(this);
    if (particle->phi > particle->h) positive_phase = true;
  }
  particle->type = REAL_PARTICLE;
  particle->Set_particle_info(this);
}
//--------------------------------------------------
// Reset_particle_force
//--------------------------------------------------
void SPH_mesh_generation::Reset_particle_force (communicator &world, int stage)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
        particle[i]->Reset_particle_force (this, stage);
    }
  }, ap);  
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Reset_particle_force finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Set_local_ghost_index
//--------------------------------------------------
void SPH_mesh_generation::Set_local_ghost_index (communicator &world)
{
#if SYM_DIM != 0
  int num_sym_particle = int(particle_sym.size());
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, num_sym_particle),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i)
      particle_sym[i]->Set_local_index (i+local_id_record);
  }, ap);
  local_id_record += num_sym_particle;
#endif
#if PERI_DIM != 0
  int num_peri_particle = int(particle_peri.size());
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, num_peri_particle),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i)
      particle_peri[i]->Set_local_index (i+local_id_record);
  }, ap);
  local_id_record += num_peri_particle;
#endif

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Set_local_ghost_index finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Prepare_for_the_current_step
//--------------------------------------------------
void SPH_mesh_generation::Prepare_for_the_current_step (communicator &world)
{
/**********************************************************/
  // JZ20181230::adding for local mesh generation
  local_id_record = total_num_particle; 
#ifdef _MPI_
  if (need_for_partition_count % 20 == 0 || need_for_partition_count == 2){
      time_rest.restart();
    Check_need_for_partition(world);
      time_for_rest[world.rank()] += time_rest.elapsed();
      Output_global_info(1, world);
      Output_runtime_info(0, world);
  }
  // JZ20181220 :: CVP is not used for dynamic partitioning
  #ifdef _CVP_LB_
  if (need_for_partition == 1){
    #ifdef _V_MASS_CENTER_
    Update_vp_position_mass_center(world);
    #endif
    Partitioning(world); need_for_partition = 0;
  }
  #endif
    time_simulation.restart();
  if (need_for_rebuild_local_tree % 10 == 0){
      time_mapping.restart();
      time_bulid_local_map.restart();
    Build_local_tree(world, 2);
      time_for_bulid_local_map[world.rank()] += time_bulid_local_map.elapsed();
      time_for_mapping[world.rank()] += time_mapping.elapsed();
  }
    time_for_simulation[world.rank()] += time_simulation.elapsed();
#endif
  /**********************************************************/
    time_simulation.restart();
    time_density.restart();
  // Refresh particle scale
  Refresh_particle_scale(world);
    time_for_density[world.rank()] += time_density.elapsed();
#ifdef _MPI_
    time_release_particle.restart();
  // Release buffer particle memory
  Release_memory_of_buffer_paticle(world);
    time_for_release_particle[world.rank()] += time_release_particle.elapsed();
#endif
    time_mapping.restart();
    time_clear_cell_list.restart();
  // reset the dynamic cell list on the tree
  Reset_cell_list_info(world);
    time_for_clear_cell_list[world.rank()] += time_clear_cell_list.elapsed();

  // map the physical domain to the tree data structure
  // and update all the tree infomation
    time_map_particle_to_tree.restart();
  Map_the_particle_to_tree(world);
  
    time_for_map_particle_to_tree[world.rank()] += time_map_particle_to_tree.elapsed();
    time_for_mapping[world.rank()] += time_mapping.elapsed();

#if PERI_DIM != 0
    time_release_particle.restart();
  // Release buffer particle memory
  Release_memory_of_pbc_paticle(world);
    time_for_release_particle[world.rank()] += time_release_particle.elapsed();
  // construct periodical particles in buffer area
    time_rest.restart();
  Construct_periodical_BC_particles(world, 1);
    time_for_rest[world.rank()] += time_rest.elapsed();

#endif
#if SYM_DIM != 0
    time_release_particle.restart();
  // Release buffer particle memory
  Release_memory_of_sbc_paticle(world);
    time_for_release_particle[world.rank()] += time_release_particle.elapsed();

    time_rest.restart();
  // construct symetric particles in buffer area
  Construct_symmetric_BC_particles(world);

    time_for_rest[world.rank()] += time_rest.elapsed();
#endif
    time_mapping.restart();
    time_update_every_level_info.restart();
  // update level info to map the cell list from low level to high level
  Update_every_level_info(1, world);

#ifdef _MPI_
  //update exchange information tag system
  Update_exchange_infor(world);
#endif
  time_for_update_every_level_info[world.rank()] += time_update_every_level_info.elapsed();
  time_for_mapping[world.rank()] += time_mapping.elapsed();
  time_for_simulation[world.rank()] += time_simulation.elapsed();
/**********************************************************/
#ifdef _MPI_
  if (need_for_rebuild_graph % 10 == 0){
      time_graph.restart();
    // Graph operation
    #ifdef _NARROW_BAND_GRAPH_
      Construct_graph_modified(world);
    #else
      Construct_graph(world);
    #endif
      time_for_graph[world.rank()] += time_graph.elapsed();
  }

    time_release_particle.restart();
  // release ghost particle memory
  Release_memory_of_ghost_paticle(world);
    time_for_release_particle[world.rank()] += time_release_particle.elapsed();
    time_communication.restart();
  // Exchange buffer particle information
  Data_communication(0, world);
    time_for_communication0[world.rank()] += time_communication.elapsed();
    time_simulation.restart();
    time_mapping.restart();
    time_clear_cell_list.restart();
  // reset the dynamic cell list on the tree

  Reset_cell_list_info(world);
    time_for_clear_cell_list[world.rank()] += time_clear_cell_list.elapsed();

  // map the physical domain to the tree data structure
  // and update all the tree infomation
    time_map_particle_to_tree.restart();
  Map_the_particle_to_tree(world);

  // map buffer particles to the tree data structure and update all the tree infomation
  Map_the_buffer_particle_to_tree(world);

#if PERI_DIM != 0
  Construct_periodical_BC_particles(world, 2);
  // map pbc particle to the tree
  Map_the_peri_particle_to_tree(world);
#endif

#if SYM_DIM != 0
  // map pbc particle to the tree
  Map_the_sym_particle_to_tree(world);
#endif
    time_for_map_particle_to_tree[world.rank()] += time_map_particle_to_tree.elapsed();
    time_update_every_level_info.restart();
  // update level info to map the cell list from low level to high level
  Update_every_level_info(1, world);

  Get_particle_info_for_ghosts(world);

    time_for_update_every_level_info[world.rank()] += time_update_every_level_info.elapsed();
    time_for_mapping[world.rank()] += time_mapping.elapsed();
    time_for_simulation[world.rank()] += time_simulation.elapsed();
#endif
/**********************************************************/
    time_rest.restart();
  // set local index for ghost particles
  Set_local_ghost_index (world);
    time_for_rest[world.rank()] += time_rest.elapsed();
    time_simulation.restart();
    time_reset_neighbor.restart();

  // reset the neighbor infomation from the tree data
  Reset_particle_neighbor_info(world);

    time_for_reset_neighbor[world.rank()] += time_reset_neighbor.elapsed();
    time_refresh_neighbor.restart();

  // get particle density(rho_n+1) and calculate hydrostatic pressure (P_n+1)
  Refresh_neighbor_info(2, world);

    time_for_refresh_neighbor[world.rank()]  += time_refresh_neighbor.elapsed();
    time_for_simulation[world.rank()] += time_simulation.elapsed();

  // Refresh_particle_infor(world, 0, 1);

#ifdef _MPI_
  need_for_partition_count ++;
  need_for_rebuild_graph ++;
  need_for_rebuild_local_tree ++;
#endif

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Prepare_for_the_current_step finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Get_kernel_summation
//--------------------------------------------------
void SPH_mesh_generation::Get_kernel_summation (communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle[i]->Get_kernel_summation(this);
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Get_kernel_summation finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Accumulate_particle_force
//--------------------------------------------------
void SPH_mesh_generation::Accumulate_particle_force (communicator &world, int stage)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      if (particle[i]->type == REAL_PARTICLE){
        particle[i]->Accumulate_particle_force_real(this, stage);
      }else if (particle[i]->type == SURFACE_PARTICLE){
        particle[i]->Accumulate_particle_force_surf(this, stage);
        particle[i]->Map_surface_particle_force(this);
      }else if (particle[i]->type == SEGMENT_PARTICLE){
        particle[i]->Accumulate_particle_force_seg(this, stage);
        particle[i]->Map_surface_particle_force(this);
      }
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Accumulate_particle_force finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Map_surface_particle_force
//--------------------------------------------------
void SPH_mesh_generation::Map_surface_particle_force (communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      if (particle[i]->type == SURFACE_PARTICLE || particle[i]->type == SEGMENT_PARTICLE){
        particle[i]->Map_surface_particle_force(this);
      }
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Map_surface_particle_force finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Update_velocity_half_step
//--------------------------------------------------
void SPH_mesh_generation::Update_velocity_half_step (communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle[i]->Update_velocity_half_step(this);
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_velocity_half_step finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Update_coord_full_step
//--------------------------------------------------
void SPH_mesh_generation::Update_coord_full_step (communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      if (particle[i]->type != SINGULARITY_PARTICLE){
        particle[i]->Update_coord_full_step(this);
      }
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_coord_full_step finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Free_dummy_particles
//--------------------------------------------------
void SPH_mesh_generation::Free_dummy_particles(std::vector<p_Particle> &particle_total, communicator &world)
{
  int num_dummy = dummy_particle.size();

  for (int i = 0; i < num_dummy; ++i)
  {
    particlepool.free(dummy_particle[i]);
  }

  dummy_particle.clear();

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Free_dummy_particles finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Reconstruce_dummy_particles
//--------------------------------------------------
void SPH_mesh_generation::Reconstruce_dummy_particles(std::vector<p_Particle> &particle_total, communicator &world)
{

  // if (run_time < v_reini_change) return;

  // static affinity_partitioner ap;
  // int num_surface = particle_surface.size();
  // parallel_for( blocked_range<int>(0, num_surface),[&](const blocked_range<int>& r){
  //   for(int i=r.begin(); i!=r.end(); ++i){
  //     // if (particle_surface[i]->type != SURFACE_PARTICLE){
  //     //   cout<<"<<<<< ERROR!!! current particle is not surface particle\n";
  //     //   world.abort(-1);
  //     // }else{

  //     // }
  //   }
  // }, ap);

  int ncount = 0;
  for (int i = 0; i < total_num_particle; ++i)
  {
    if (particle[i]->type == SURFACE_PARTICLE || particle[i]->type == SEGMENT_PARTICLE || particle[i]->type == SINGULARITY_PARTICLE){
      p_Particle tmp;
      tmp = particlepool.malloc();
      bool pass = false;

      tmp->Initialize (this, ncount+total_num_particle, DUMMY_PARTICLE, world.rank());

      pass = particle[i]->Reconstruce_dummy_particles(this, tmp);

      if (pass){
        dummy_particle.push_back(tmp);
        tmp->Set_local_index ( ncount+local_id_record);
        ncount++;        
      }else{
        particlepool.free(tmp);
      }
    }
  }

  local_id_record += ncount;

  // for (int i = 0; i < ncount; ++i)
  // {
  //   particle_total.push_back (dummy_particle[i]);
  // }

  // cout<<"dummy_particle.size(): "<<dummy_particle.size()<<endl;
  // cout<<"particle_total.size(): "<<particle_total.size()<<endl;

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Reconstruce_dummy_particles finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Post processing
//--------------------------------------------------
void SPH_mesh_generation::Post_processing(int _nout, communicator &world)
{
  Save_particle_post_file (_nout, world);

  if (run_time < t_for_post) return;

  Find_and_map_particle (world);

  Get_particle_info (world);
    
  Set_particle_scale (world);
  
  Prepare_for_the_current_step (world);

  // prepare the total particle array
  std::vector<p_Particle> particle_total;
  for(int i = 0; i < total_num_particle; i++){
    particle_total.push_back(particle[i]);
  }

  #ifdef _MPI_
  int num_ghost_buffer_particle = particle_buffer.size();
  for(int i = 0; i < num_ghost_buffer_particle; i++){
    particle_total.push_back(particle_buffer[i]);
  }
  #endif
  
  #if SYM_DIM != 0
  int num_ghost_bc_particle = particle_sym.size();
  
  for(int i = 0; i < num_ghost_bc_particle; i++){
    particle_total.push_back(particle_sym[i]);
  }
  #endif
  #if PERI_DIM != 0
  int num_ghost_bc_particle = particle_peri.size();
  for(int i = 0; i < num_ghost_bc_particle; i++){
    particle_total.push_back(particle_peri[i]);
  }
  #endif
  
  mesh.Remove_all_edges (world);
  mesh.Remove_all_verts (world);
    
  mesh.Initialize (particle_total.size(), world);

  if (DIM == 2){
    // get local voronoi tesselation and add edges

    static affinity_partitioner ap;
    parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
          particle[i]->Local_voronoi_diagram(this);
      }
    }, ap);
    #ifdef _MPI_
    parallel_for( blocked_range<int>(0, int(num_ghost_buffer_particle)),[&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
          particle_buffer[i]->Local_voronoi_diagram(this);
      }
    }, ap);
    #endif
    // construct faces
    mesh.Reconstruct_tris (particle_total, this, world);

  }else if (DIM == 3){

    #ifdef _GEN_DUMMY_PART_
      Reconstruce_dummy_particles (particle_total, world);
    #endif

    // construct faces
    mesh.Reconstruct_tets (particle_total, this, world);

    Reset_particle_array (world);
    
    Check_total_number_of_particle (world);
  }

  // cout<<"dummy_particle.size(): "<<dummy_particle.size()<<endl;
  // cout<<"particle_total.size(): "<<particle_total.size()<<endl;

  // Output mesh
  mesh_writer.Output_mesh (_nout, particle_total, this, world);

  // Output report
  if (DIM == 2)       mesh.Write_tri_quality_plt (n_post, run_time, world);
  else if (DIM == 3)  {
    mesh.Write_tet_quality_plt (n_post, run_time, world);

    #ifdef _GEN_DUMMY_PART_
      Free_dummy_particles (particle_total, world);
    #endif

  }

  particle_total.clear();

  n_post ++;

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Post_processing finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Output particle infomation
//--------------------------------------------------
void SPH_mesh_generation::Output_paraview_file(int n, int flag, communicator &world)
{
  // prepare the total particle array
  // TODO::make it more efficient
  std::vector<p_Particle> particle_total;
  for(int i = 0; i < total_num_particle; i++){
    particle_total.push_back(particle[i]);
  }

//   #ifdef _MPI_
//   int num_ghosts = particle_buffer.size();
//   for(int i = 0; i < num_ghosts; i++){
//     particle_total.push_back(particle_buffer[i]);
//   }
//   #endif

  particle_writer.Output_particle (n, particle_total, this, world);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Output_paraview_file "<<n<<" finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Output particle infomation
//--------------------------------------------------
void SPH_mesh_generation::Output_plt_file(int n, int flag, communicator &world)
{
  // flag = 1: normal output
  // flag = 2: check partitioning
  FILE    *fp;
  char    filename[256];
  if (flag == 1)
    sprintf(filename,"%s%d%s%d%s","./outdata/simulation_sp.",n,".",world.rank(),".plt");
  else if (flag == 2){
    sprintf(filename,"%s%d%s%d%s%d%s","./outdata/CVP_particle_iter.", i_iter,".", n,".",world.rank(),".plt");
  }

  ofstream out(filename, ios::trunc);
//  if (world.rank() == 0)
  out<<"VARIABLES = \"x\",\"y\",\"z\",\"local_id\",\"color\",\"level\",\"scale\",\"type\",\"NN\",\"n_x\",\"n_y\",\"a_x\",\"a_y\",\"rho\",\"h\",\"phi\",\"curv\"\n";
  for (int i = 0; i < total_num_particle; i++){
    Real     v = 0.0;
    Particle *current_particle = particle[i];
    my_real p_coord; my_set_const(p_coord, 0.);
    my_shift_coordinate(p_coord, current_particle->coord, domain, box_l, box_r);
    out<<p_coord.i<<" "
       <<p_coord.j<<"  "
       <<p_coord.k<<"  "
       <<current_particle->local_id<<" "
       <<current_particle->color<<" "
       <<current_particle->level<<" "
       <<current_particle->scale<<" "
       <<current_particle->type<<" "
       <<current_particle->neighbor.size()<<" "
       <<current_particle->norm.i<<"  "
       <<current_particle->norm.j<<"  "
       <<current_particle->a.i<<"  "
       <<current_particle->a.j<<"  "
       <<current_particle->rho<<"  "
       <<current_particle->h<<"  "
       <<current_particle->phi<<"  "       
       <<current_particle->curv<<"\n";
  }
  out.close();
}
#if SYM_DIM != 0
//--------------------------------------------------
// Output particle infomation
//--------------------------------------------------
void SPH_mesh_generation::Output_SBC_plt_file(int n, communicator &world)
{
  FILE    *fp;
  char    filename[256];
    sprintf(filename,"%s%d%s%d%s","./outdata/SBC_particle.",n,".",world.rank(),".plt");

  int total_num_sym_particle = int(particle_sym.size());

  ofstream out(filename, ios::trunc);

  out<<"VARIABLES = \"x\",\"y\",\"z\",\"local_id\",\"color\",\"level\",\"scale\",\"type\",\"NN\",\"n_x\",\"n_y\",\"a_x\",\"a_y\",\"rho\",\"h\",\"phi\",\"curv\"\n";
  for (int i = 0; i < total_num_sym_particle; i++){
    Real     v = 0.0;
    Particle *current_particle = particle_sym[i];
    my_real p_coord; my_set_const(p_coord, 0.);
    my_shift_coordinate(p_coord, current_particle->coord, domain, box_l, box_r);
    out<<p_coord.i<<" "
       <<p_coord.j<<"  "
       <<p_coord.k<<"  "
       <<current_particle->local_id<<" "
       <<current_particle->color<<" "
       <<current_particle->level<<" "
       <<current_particle->scale<<" "
       <<current_particle->type<<" "
       <<current_particle->neighbor.size()<<" "
       <<current_particle->norm.i<<"  "
       <<current_particle->norm.j<<"  "
       <<current_particle->a.i<<"  "
       <<current_particle->a.j<<"  "
       <<current_particle->rho<<"  "
       <<current_particle->h<<"  "
       <<current_particle->phi<<"  "       
       <<current_particle->curv<<"\n";
  }
  out.close();
}
#endif
//--------------------------------------------------
// Output simulation infor
//--------------------------------------------------
void SPH_mesh_generation::Output_simulation_infor(int n, communicator &world)
{
  Get_energy(world);
  
  //TODO
//   Get_max_min_h(world);
  
  if (world.rank() == 0){
    char    filename[256];
    sprintf(filename,"%s","./outdata/simulation_infor.csv");
    ofstream out(filename, ios::app);
    if (n == 0){
      out<<"\"run_time\",\"max_v\",\"E_total\",\"E_k\",\"np_real\",\"np_surf\",\"np_seg\",\"np_sing\",\"calc_real\",\"calc_surf\",\"calc_seg\",\"calc_sing\",\"fill_coeff\",\"current_stage\",\"count_full\",\"max_count_relax\"\n";
    }

    out<<run_time<<","
       <<glbl_max_v<<","
       <<glbl_Etotal<<","
       <<glbl_Ek<<","
       <<current_glbl_total_num_particle_mesh<<","
       <<current_glbl_total_num_particle_surface<<","
       <<current_glbl_total_num_particle_segment<<","
       <<current_glbl_total_num_particle_singularity<<","
       <<glbl_total_num_particle_mesh<<","
       <<glbl_total_num_particle_surface<<","
       <<glbl_total_num_particle_segment<<","
       <<glbl_total_num_particle_singularity<<","
       <<fill_coeff<<","
       <<current_stage<<","
       <<count_full<<","
       <<max_count_relax<<"\n";
  }
}
//--------------------------------------------------
// Save_particle_post_file
//--------------------------------------------------
void SPH_mesh_generation::Save_particle_post_file(int n, communicator &world)
{
  if (world.rank() == 0){
    char    fn[256];
    sprintf(fn,"%s","./rstfile_particles/list.csv");
    if (n == 0){
      ofstream out_list(fn, ios::trunc);
      // out_list<<"#nranks"<<"\n";
      out_list<<world.size()<<"\n";
      // out_list<<"#step and number of particles"<<"\n";
      out_list<<int(run_time)<<" "<<glbl_total_num_particle<<"\n";
      out_list.close();
    }else{
      ofstream out_list(fn, ios::app);
      out_list<<int(run_time)<<" "<<glbl_total_num_particle<<"\n";
      out_list.close();
    }
  }

  char pp[256];
  sprintf(pp,"%s%d%s%d%s","./rstfile_particles/step_",int(run_time),"_rank_",world.rank(),".pst");

  ofstream out_pp(pp, ios::trunc | ios::binary);

  out_pp.write(reinterpret_cast<char*>(&total_num_particle), sizeof(int));
  for (int i = 0; i < total_num_particle; ++i)
  {
    out_pp.write(reinterpret_cast<char*>(&particle[i]->coord.i), sizeof(Real));
    out_pp.write(reinterpret_cast<char*>(&particle[i]->coord.j), sizeof(Real));
    out_pp.write(reinterpret_cast<char*>(&particle[i]->coord.k), sizeof(Real));
  }

  out_pp.close();

  world.barrier();

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Save_particle_post_file finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Save_restart_file
//--------------------------------------------------
void SPH_mesh_generation::Save_restart_file(communicator &world)
{
  FILE    *fp;
  char    fn[256];
  sprintf(fn,"%s","./rstfile/levelset.rst");
  
  ofstream out_(fn, ios::trunc | ios::binary);
  
  // out_<<"//<<<<restart parameters: DIM | domain | box_l | box_r | particle_box | particle_box_l | particle_box_r | resolution | Lmin | Lmax | ini_num_cell | num_level \n";
  
  // out_<<int(DIM)<<"\n"
  //       <<domain.i<<" "<<domain.j<<" "<<domain.k<<"\n"
  //       <<box_l.i<<" "<<box_l.j<<" "<<box_l.k<<"\n"
  //       <<box_r.i<<" "<<box_r.j<<" "<<box_r.k<<"\n"
  //       <<particle_box.i<<" "<<particle_box.j<<" "<<particle_box.k<<"\n"
  //       <<particle_box_l.i<<" "<<particle_box_l.j<<" "<<particle_box_l.k<<"\n"
  //       <<particle_box_r.i<<" "<<particle_box_r.j<<" "<<particle_box_r.k<<"\n"
  //       <<resolution.i<<" "<<resolution.j<<" "<<resolution.k<<"\n"
  //       <<Lmin<<" "<<Lmax<<"\n"
  //       <<ini_num_cell.i<<" "<<ini_num_cell.j<<" "<<ini_num_cell.k<<"\n"
  //       <<num_level<<"\n";

  int var = int(DIM);
  out_.write(reinterpret_cast<char*>(&var)             , sizeof(int));
  out_.write(reinterpret_cast<char*>(&domain.i)        , sizeof(Real));
  out_.write(reinterpret_cast<char*>(&domain.j)        , sizeof(Real));
  out_.write(reinterpret_cast<char*>(&domain.k)        , sizeof(Real));
  out_.write(reinterpret_cast<char*>(&box_l.i)         , sizeof(Real));
  out_.write(reinterpret_cast<char*>(&box_l.j)         , sizeof(Real));
  out_.write(reinterpret_cast<char*>(&box_l.k)         , sizeof(Real));
  out_.write(reinterpret_cast<char*>(&box_r.i)         , sizeof(Real));
  out_.write(reinterpret_cast<char*>(&box_r.j)         , sizeof(Real));
  out_.write(reinterpret_cast<char*>(&box_r.k)         , sizeof(Real));
  out_.write(reinterpret_cast<char*>(&particle_box.i)  , sizeof(Real));
  out_.write(reinterpret_cast<char*>(&particle_box.j)  , sizeof(Real));
  out_.write(reinterpret_cast<char*>(&particle_box.k)  , sizeof(Real));
  out_.write(reinterpret_cast<char*>(&particle_box_l.i), sizeof(Real));
  out_.write(reinterpret_cast<char*>(&particle_box_l.j), sizeof(Real));
  out_.write(reinterpret_cast<char*>(&particle_box_l.k), sizeof(Real));
  out_.write(reinterpret_cast<char*>(&particle_box_r.i), sizeof(Real));
  out_.write(reinterpret_cast<char*>(&particle_box_r.j), sizeof(Real));
  out_.write(reinterpret_cast<char*>(&particle_box_r.k), sizeof(Real));
  out_.write(reinterpret_cast<char*>(&resolution.i)    , sizeof(int ));
  out_.write(reinterpret_cast<char*>(&resolution.j)    , sizeof(int ));
  out_.write(reinterpret_cast<char*>(&resolution.k)    , sizeof(int ));
  out_.write(reinterpret_cast<char*>(&Lmin)            , sizeof(int ));
  out_.write(reinterpret_cast<char*>(&Lmax)            , sizeof(int ));
  out_.write(reinterpret_cast<char*>(&ini_num_cell.i)  , sizeof(int ));
  out_.write(reinterpret_cast<char*>(&ini_num_cell.j)  , sizeof(int ));
  out_.write(reinterpret_cast<char*>(&ini_num_cell.k)  , sizeof(int ));
  out_.write(reinterpret_cast<char*>(&num_level)       , sizeof(int ));

  out_.flush();
  
  level_set.Save_lset_restart_file (fn, this, world);
  
  cout<<"<<<<< Restart file is writen...."<<endl;
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Restart file saved\n";
  out.close();
#endif
}
#ifdef _MPI_
//--------------------------------------------------
// Get_particle_info_for_ghosts
//--------------------------------------------------
void SPH_mesh_generation::Get_particle_info_for_ghosts(communicator &world)
{
  static affinity_partitioner ap;
  int num_ghosts = particle_buffer.size();

  parallel_for( blocked_range<int>(0, num_ghosts),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle_buffer[i]->Calculate_particle_infor(this);
      particle_buffer[i]->Set_local_index (i+local_id_record);
    }
  }, ap);
  
  local_id_record += num_ghosts;

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Get_particle_info_for_ghosts finished\n";
  out.close();
#endif
}
#endif
