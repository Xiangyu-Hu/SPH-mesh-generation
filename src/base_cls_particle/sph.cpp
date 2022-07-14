#include <cmath>
#include "glbfunc.h"
#include "sph.h"

/***************************************************/
/*                                                 */
/*           Functions defined in class "SPH"      */
/*                                                 */
/***************************************************/

//--------------------------------------------------
// initialize system for the 1st time
//--------------------------------------------------
void SPH::Initialize_system(communicator &world){
#ifdef _MPI_
  need_for_partition = 0;
  error_tolerance = 0.1;
  imbalance_tolerance  = 0.1;
  // if number of cpu is larger than 1, graph is constructed for data communication
  sph_graph.Initialize(world);
#endif

  // initialize level info
  level_info = new p_Level_info[num_level];
  for(int i =0; i<num_level; i++)
    level_info[i] = new Level_info;
#ifndef _MPI_
  for(int i =0; i<num_level; i++)
    level_info[i]->Initialize(i, this, world);
  Refresh_particle_scale(world);
#else
  Build_local_tree(world, 1);
  need_for_rebuild_local_tree ++;
  Refresh_particle_scale(world);
#endif
#ifdef _MPI_
  #if defined(_DIST_NUMBER_) || defined(_WEIGHTED_PARTITION_)
  // reset the dynamic cell list on the tree
  Reset_cell_list_info(world);

  // map the physical domain to the tree data structure 
  // and update all the tree infomation
  Map_the_particle_to_tree(world);

  // update level info to map the cell list from low level to high level
  Update_every_level_info(1, world);
    #ifdef _SCLL_
    Build_every_level_subcell_list(world);
    #endif
  #endif

  Get_patitioning_mass(0, world);

  #ifdef _WEIGHTED_PARTITION_
  glbl_total_neighbor_size = Get_global_neighbor_size(world, 0);
  glbl_total_dist_number = Get_global_dist_number(world);
  weight = 0.5;
  Get_weighted_patitioning_mass(0, world);
  #endif

  total_mass = Get_total_p_mass_local(world);

  glbl_total_mass = 0.;
  all_reduce(world, total_mass, glbl_total_mass, std::plus<Real>());
  total_mass /= glbl_total_mass;

  Normalize_particle_mass(world);

  glbl_total_mass = 1.;

  #ifdef _V_MASS_CENTER_
    Update_vp_position_mass_center(world);
  #endif

  // Partitioning(world);
  need_for_partition_count ++;
#endif
}
//--------------------------------------------------
// define level infor
//--------------------------------------------------
void SPH::Define_level_infor (communicator &world)
{
  Real scale_temp = domain.i;
#if DIM_Y == 1
  scale_temp = AMIN1(domain.j,scale_temp);
#endif
#if DIM_Z == 1
  scale_temp = AMIN1(domain.k,scale_temp);
#endif

  while (scale_temp > ini_scale){
    scale_temp /= SCALE_RATIO;
  }
  scale_temp    *= SCALE_RATIO;
  ini_scale      = scale_temp;
  max_scale       = ini_scale;
  
  num_level = 1;
  while ((scale_temp+1.e-6) > min_scale){
    scale_temp /= SCALE_RATIO;
    num_level ++;
  }
  min_scale       = scale_temp*SCALE_RATIO;
  num_level      -= 1;
  Lmax            = num_level-1;
  Lmin            = 0;

  my_real scale_axis;
  my_set_data (scale_axis, domain);
#if DIM_X == 1
  while (scale_axis.i >= max_scale){
    scale_axis.i /= SCALE_RATIO;
  }
  scale_axis.i    *= SCALE_RATIO;
#endif
#if DIM_Y == 1
  while (scale_axis.j >= max_scale){
    scale_axis.j /= SCALE_RATIO;
  }
  scale_axis.j    *= SCALE_RATIO;
#endif
#if DIM_Z == 1
  while (scale_axis.k >= max_scale){
    scale_axis.k /= SCALE_RATIO;
  }
  scale_axis.k    *= SCALE_RATIO;
#endif

  ini_num_cell.i = DIM_X==1 ? int(domain.i/scale_axis.i): 1;
  ini_num_cell.j = DIM_Y==1 ? int(domain.j/scale_axis.j): 1;
  ini_num_cell.k = DIM_Z==1 ? int(domain.k/scale_axis.k): 1;
  
  if (world.rank() == 0){
    cout<<"<<<<< min_scale : "<<min_scale<<"\n";
    cout<<"<<<<< max_scale : "<<max_scale<<"\n";
    cout<<"<<<<< Lmax : "<<Lmax<<"\n";
    cout<<"<<<<< Lmin : "<<Lmin<<"\n";
  }
}
//--------------------------------------------------
// Define_level_infor_tight
//--------------------------------------------------
void SPH::Define_level_infor_tight (communicator &world)
{
  Real scale_temp = domain.i;
#if DIM_Y == 1
  scale_temp = AMIN1(domain.j,scale_temp);
#endif
#if DIM_Z == 1
  scale_temp = AMIN1(domain.k,scale_temp);
#endif

  int temp;
  temp = floor(scale_temp/ini_scale);
  ini_scale = scale_temp/Real(temp);
  max_scale = ini_scale;
  scale_temp = ini_scale;
  
  num_level = 1;
  while ((scale_temp+1.e-6) > min_scale){
    scale_temp /= SCALE_RATIO;
    num_level ++;
  }
  min_scale       = scale_temp*SCALE_RATIO;
  num_level      -= 1;
  Lmax            = num_level-1;
  Lmin            = 0;

  ini_num_cell.i = DIM_X==1 ? int(domain.i/max_scale): 1;
  ini_num_cell.j = DIM_Y==1 ? int(domain.j/max_scale): 1;
  ini_num_cell.k = DIM_Z==1 ? int(domain.k/max_scale): 1;
  
  if (world.rank() == 0){
    cout<<"<<<<< min_scale : "<<min_scale<<"\n";
    cout<<"<<<<< max_scale : "<<max_scale<<"\n";
    cout<<"<<<<< Lmax      : "<<Lmax<<"\n";
    cout<<"<<<<< Lmin      : "<<Lmin<<"\n";
  }
}
//--------------------------------------------------
// Set current time
//--------------------------------------------------
void SPH::Set_current_time(Real Time){
  run_time = Time;
}
//--------------------------------------------------
// Set_timestep
//--------------------------------------------------
void SPH::Set_timestep(communicator &world)
{
  local_timestep = 1.e20;
  max_v    = -1.e20;

  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      Real temp;
      Particle *current = particle[i];
      current->Set_timestep();
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
// Reset_timestep
//--------------------------------------------------
void SPH::Reset_timestep(Real timestep_new)
{
  glbl_timestep = timestep_new;
}
//--------------------------------------------------
// Get_timestep
//--------------------------------------------------
Real SPH::Get_timestep()
{
  return glbl_timestep;
}
//--------------------------------------------------
// reset all the cell list
//--------------------------------------------------
void SPH::Reset_cell_list_info(communicator &world)
{
  for (int i = 0; i < num_level; i++)
    level_info[i]->Reset_cell_list_info(world);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<reset all the cell list finished\n";
  out.close();
#endif
}
//--------------------------------------------------
//  map the particle to the tree according 
//  to the physical position
//--------------------------------------------------
void SPH::Map_the_particle_to_tree(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int current_level = 0;
      int flag = 0;
      for( int j = num_level-1; j>=0; j--){
        if(fabs(particle[i]->scale - level_info[j]->scale) <= 1.e-10){
          flag = 1;
          current_level = j;
        }
        else if (particle[i]->scale > level_info[j]->scale &&
                 particle[i]->scale < level_info[j-1]->scale){
          flag = 1;
          current_level = j-1;
        }
        if (flag == 1) break;
      }
      level_info[current_level]->Update_cell_list(particle[i], world);
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Map_the_particle_to_tree finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// update the tree data structure
//--------------------------------------------------
void SPH::Update_every_level_info(int flag, communicator &world)
{
  for( int i=num_level-1; i>=0; i--)
    level_info[i]->Reset_tags(world);

    level_info[num_level-1]->Update_leaf_particle_status();

  for( int i=num_level-2; i>=0; i--){
    level_info[i+1]->Update_exist_status( );
    level_info[i]->Update_leaf_particle_status();
    level_info[i]->Update_every_level_info(flag, level_info, i+1, world);
  }
  level_info[0]->Update_exist_status( );

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_every_level_info finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// reset all the particle neighbor info
//--------------------------------------------------
void SPH::Reset_particle_neighbor_info(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle[i]->Reset_neighbor_info();
    }
  }, ap);
#if SYM_DIM != 0
  parallel_for( blocked_range<int>(0, int(particle_sym.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle_sym[i]->Reset_neighbor_info();
    }
  }, ap);
#endif
#if PERI_DIM != 0
  parallel_for( blocked_range<int>(0, int(particle_peri.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle_peri[i]->Reset_neighbor_info();
    }
  }, ap);
#endif
#ifdef _MPI_
  parallel_for( blocked_range<int>(0, int(particle_buffer.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle_buffer[i]->Reset_neighbor_info();
    }
  }, ap);
#endif
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Reset_particle_neighbor_info finished\n";
  out.close();
#endif
}
//--------------------------------------------------------
// refresh the neighbor infomation for all the particles
//-------------------------------------------------------
void SPH::Refresh_neighbor_info(int flag, communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle[i]->Refresh_neighbor_info(this, flag);
    }
  }, ap);
#if SYM_DIM != 0
  parallel_for( blocked_range<int>(0, int(particle_sym.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle_sym[i]->Refresh_neighbor_info(this, flag);
    }
  }, ap);
#endif
#if PERI_DIM != 0
  parallel_for( blocked_range<int>(0, int(particle_peri.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle_peri[i]->Refresh_neighbor_info(this, flag);
    }
  }, ap);
#endif
#ifdef _MPI_
  parallel_for( blocked_range<int>(0, int(particle_buffer.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle_buffer[i]->Refresh_neighbor_info(this, flag);
    }
  }, ap);
#endif
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Refresh_neighbor_info finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Add a coarse level at the top of the list
//--------------------------------------------------
void SPH::Add_coarse_level(communicator &world){

  // reset array for level_info
  p_Level_info    *p_level_tmp;
  p_level_tmp    = new p_Level_info[num_level+1];
  p_level_tmp[0] = new Level_info;
  for( int i=0; i < num_level; i++){
    p_level_tmp[i+1] = level_info[i];
  }
  num_level ++;
  Lmin --;
  max_scale *= SCALE_RATIO;
  p_level_tmp[0]->Initialize(Lmin, this, world);

  delete [] level_info;
  level_info = new p_Level_info[num_level];

  for( int i=0; i < num_level; i++){
    level_info[i] = p_level_tmp[i];
    level_info[i]->Lmin = Lmin;
  }
  delete [] p_level_tmp;

#if PERI_DIM != 0 || SYM_DIM != 0 || defined(_MPI_)
  for( int i=1; i < num_level; i++){
    // multi_array need to be resized
    level_info[i]->Reinitialize_cell_list(world, this);
  }
#endif
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Add_coarse_level finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Add a finer level at the bottem of the list
//--------------------------------------------------
void SPH::Add_finer_level(communicator &world){
  // reset array for level_info
  p_Level_info    *p_level_tmp;
  p_level_tmp    = new p_Level_info[num_level+1];
  for( int i=0; i < num_level; i++){
    p_level_tmp[i] = level_info[i];
  }
  p_level_tmp[num_level] = new Level_info;
  num_level ++;
  Lmax ++;
  min_scale /= SCALE_RATIO;
  p_level_tmp[num_level-1]->Initialize(Lmax, this, world);

  delete [] level_info;
  level_info = new p_Level_info[num_level];
  for( int i=0; i < num_level; i++){
    level_info[i] = p_level_tmp[i];
    level_info[i]->Lmax = Lmax;
  }
  delete [] p_level_tmp;
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Add_finer_level finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Refresh particle scale
//--------------------------------------------------
void SPH::Refresh_particle_scale(communicator &world)
{
  int flag = 0;
  //flag used to decide whether a new level is required
  //-1 : add a coarse level
  //1  : add a finer level
  //0  : not required

  for(int i=0; i!=total_num_particle; ++i){
    particle[i]->Set_particle_scale(this, flag);
  }

  int glbl_flag_1 = 0;
  int glbl_flag_2 = 0;
  all_reduce(world, flag, glbl_flag_1, mpi::minimum<int>());
  all_reduce(world, flag, glbl_flag_2, mpi::maximum<int>());
  
  if (glbl_flag_1 == -1){
    Add_coarse_level(world);
  }
  if (glbl_flag_2 == 1){
    Add_finer_level(world);
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Refresh particle scale finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Refresh particle information
//--------------------------------------------------
void SPH::Refresh_particle_infor(communicator &world, int flag1, int flag2)
{
   time_simulation.restart();
    time_rest.restart();
#if PERI_DIM != 0
  // Refresh periodic bc particle information
  Refresh_pbc_particle_info(world, flag1);
#endif
#if SYM_DIM != 0
  // Refresh symmetric bc particle information
  Refresh_sbc_particle_info(world, flag1);
#endif
    time_for_rest[world.rank()] += time_rest.elapsed();
    time_for_simulation[world.rank()] += time_simulation.elapsed();

#ifdef _MPI_
    time_communication.restart();

  // Exchange buffer particle information
  Data_communication(flag2, world);

  if (flag2 == 0) time_for_communication0[world.rank()] += time_communication.elapsed();
  else if (flag2 == 1) time_for_communication1[world.rank()] += time_communication.elapsed();
  else if (flag2 == 2) time_for_communication2[world.rank()] += time_communication.elapsed();

#if PERI_DIM != 0
  // Refresh periodic bc particle information
  Refresh_pbc_particle_info(world, flag1);
#endif

#endif
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Refresh_particle_infor finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Get total energy of the system
//--------------------------------------------------
void SPH::Get_energy(communicator &world)
{
  Ek = Ep = Ef = Etotal = 0.;
  for( int i=0; i<total_num_particle; i++){
    Real iEk, iEp, iEf;
    iEk = iEp = iEf = 0.;
    particle[i]->Get_energy(iEk, iEp, iEf);
    Ek += iEk;
    Ep += iEp;
    Ef += iEf;
  }
  Etotal = Ek + Ep + Ef;

  if (world.rank() == 0){
    reduce( world, Etotal, glbl_Etotal, std::plus<Real>(), 0);
    reduce( world, Ek, glbl_Ek, std::plus<Real>(), 0);
    cout<<"<<<Total energy: "<<glbl_Etotal<<" Kinetic energy: "<<glbl_Ek<<"\n";
  }
  else{
    reduce( world, Etotal, std::plus<Real>(), 0);
    reduce( world, Ek, std::plus<Real>(), 0);
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Get_energy finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Get local total force
//--------------------------------------------------
void SPH::Check_total_force(communicator &world)
{
  my_real local_total_force;
  my_set_const(local_total_force, 0.);
  my_set_const(glbl_total_force, 0.);

  ParallelAccumulation(total_num_particle, particle, local_total_force);

  reduce(world, local_total_force.i, glbl_total_force.i, std::plus<Real>(), 0);
  reduce(world, local_total_force.j, glbl_total_force.j, std::plus<Real>(), 0);
  reduce(world, local_total_force.k, glbl_total_force.k, std::plus<Real>(), 0);

  if (world.rank() == 0){
    glbl_force_error = get_distance (glbl_total_force);
    if (glbl_force_error > 1.e-10){
      my_cout (glbl_total_force, "<< Global total force ");
      world.abort(-1);
    }
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Check_total_force finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Output runtime information
//--------------------------------------------------
void SPH::Output_runtime_info(int flag, communicator &world){
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s","./outdata/runtime_info.plt");
  ofstream out(filename, ios::app);

  if (world.rank() == 0){
    Real tmp[13];
    tmp[0] = time_for_graph[world.rank()];
    tmp[1] = time_for_partition[world.rank()];
    tmp[2] = time_for_communication0[world.rank()];
    tmp[3] = time_for_communication1[world.rank()];
    tmp[4] = time_for_communication2[world.rank()];
    tmp[5] = time_for_communication_partition[world.rank()];
    tmp[6] = time_for_simulation[world.rank()];
    tmp[7] = time_for_mapping[world.rank()];
    tmp[8] = time_for_refresh_neighbor[world.rank()];
    tmp[9] = time_for_force_calculation[world.rank()];
    tmp[10] = time_for_density[world.rank()];
    tmp[11] = time_for_rest[world.rank()];
    tmp[12] = time_for_update[world.rank()];
    gather(world, tmp[0], time_for_graph, 0);
    gather(world, tmp[1], time_for_partition, 0);
    gather(world, tmp[2], time_for_communication0, 0);
    gather(world, tmp[3], time_for_communication1, 0);
    gather(world, tmp[4], time_for_communication2, 0);
    gather(world, tmp[5], time_for_communication_partition, 0);
    gather(world, tmp[6], time_for_simulation, 0);
    gather(world, tmp[7], time_for_mapping, 0);
    gather(world, tmp[8], time_for_refresh_neighbor, 0);
    gather(world, tmp[9], time_for_force_calculation, 0);
    gather(world, tmp[10], time_for_density, 0);
    gather(world, tmp[11], time_for_rest, 0);
    gather(world, tmp[12], time_for_update, 0);
#ifdef _MPI_
    int *n_buffer;
    n_buffer = new int[world.size()+1];
    int total_buffer_size = int(particle_buffer.size());
#if PERI_DIM != 0
    total_buffer_size += particle_peri.size();
#endif
#if SYM_DIM != 0
    total_buffer_size += particle_sym.size();
#endif
    gather(world, total_buffer_size, n_buffer, 0);

    n_buffer[world.size()] = 0;
    for (int i = 0; i < world.size(); i++){
      n_buffer[world.size()] += n_buffer[i];
    }
    n_buffer[world.size()] /= world.size();
#endif
    out<<"VARIABLES = \"icpu\",\"t_graph\",\"t_partition\",\"t_communication0\",\"t_communication1\",\"t_communication2\",\"t_communication_partition\",\"t_simulation\",\"t_mapping\",\"t_refresh_neighbor\",\"t_force\",\"t_density\",\"t_rest\",\"t_update\",\"num_buffer\"\n";
    for(int i = 0; i <= world.size(); i++){
      out<<i<<" "
         <<time_for_graph[i]<<" "
         <<time_for_partition[i]<<" "
         <<time_for_communication0[i]<<" "
         <<time_for_communication1[i]<<" "
         <<time_for_communication2[i]<<" "
         <<time_for_communication_partition[i]<<" "
         <<time_for_simulation[i]<<" "
         <<time_for_mapping[i]<<" "
         <<time_for_refresh_neighbor[i]<<" "
         <<time_for_force_calculation[i]<<" "
         <<time_for_density[i]<<" "
         <<time_for_rest[i]<<" "
         <<time_for_update[i]<<" "
#ifdef _MPI_
         <<n_buffer[i]<<"\n";
#else
         <<"0\n";
#endif
    }
    out<<"\n";
  }else{
    gather(world, time_for_graph[world.rank()], 0);
    gather(world, time_for_partition[world.rank()], 0);
    gather(world, time_for_communication0[world.rank()], 0);
    gather(world, time_for_communication1[world.rank()], 0);
    gather(world, time_for_communication2[world.rank()], 0);
    gather(world, time_for_communication_partition[world.rank()], 0);
    gather(world, time_for_simulation[world.rank()], 0);
    gather(world, time_for_mapping[world.rank()], 0);
    gather(world, time_for_refresh_neighbor[world.rank()], 0);
    gather(world, time_for_force_calculation[world.rank()], 0);
    gather(world, time_for_density[world.rank()], 0);
    gather(world, time_for_rest[world.rank()], 0);
    gather(world, time_for_update[world.rank()], 0);
#ifdef _MPI_
    int total_buffer_size = int(particle_buffer.size());
#if PERI_DIM != 0
    total_buffer_size += particle_peri.size();
#endif
#if SYM_DIM != 0
    total_buffer_size += particle_sym.size();
#endif
    gather(world, total_buffer_size, 0);
#endif
  }
  out.close();
}
#ifdef _SCLL_
//--------------------------------------------------
// traverse the tree to build subcell list
//--------------------------------------------------
void SPH::Build_every_level_subcell_list(communicator &world){

  for( int i = 0; i < num_level; i++)
    level_info[i]->Build_subcell_list(world);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Build_every_level_subcell_list\n";
  out.close();
#endif
}
#endif
#if SYM_DIM != 0
//--------------------------------------------------
// release pbc particle memory
//--------------------------------------------------
void SPH::Release_memory_of_sbc_paticle(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<concurrent_vector<p_Particle>::iterator>(particle_sym.begin(), particle_sym.end()),
    [&](const blocked_range<concurrent_vector<p_Particle>::iterator>& r){
    for(concurrent_vector<p_Particle>::iterator it=r.begin(); it!=r.end(); ++it){
      if (NULL != *it){
        delete (*it); // release memory of exchange_particles in heap
        *it = NULL;
      }
    }  
  }, ap);

  particle_sym.clear();
  particle_sym.shrink_to_fit();
  if (int(particle_sym.capacity()) != 0){
    cout<<"Ghost particle size is wrong!!!\n"; world.abort(-1);
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Release_memory_of_sbc_paticle finished\n";
  out<<"particle_sym.size "<<particle_sym.size()<<" capacity: "<<particlepool.capacity()<<" available nodes "<<particlepool.available_node()<<" \n";
  out.close();
#endif
}
//--------------------------------------------------
// construct symmetric particles in buffer area
//--------------------------------------------------
void SPH::Construct_symmetric_BC_particles(communicator &world)
{
  for( int i=0; i<num_level; i++)
    level_info[i]->Construct_symmetric_BC_particles(particle_sym);

  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(particle_sym.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i)
      particle_sym[i]->id += glbl_total_num_particle;
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Construct_symmetric_BC_particles finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Refresh symmetric bc particle information
//--------------------------------------------------
void SPH::Refresh_sbc_particle_info(communicator &world, int flag)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(particle_sym.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      p_Particle cp = particle_sym[i]->copy[0];
      particle_sym[i]->Refresh_sbc_particle_info(cp, flag, this);
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Refresh_sbc_particle_info finished\n";
  out.close();
#endif
}
//--------------------------------------------------
//  map the sym particle to the tree according 
//  to the physical position
//--------------------------------------------------
void SPH::Map_the_sym_particle_to_tree(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(particle_sym.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int current_level = 0;
      int flag = 0;
      for( int j = num_level-1; j>=0; j--){
        if(particle_sym[i]->scale == level_info[j]->scale){
          flag = 1;
          current_level = j;
        }
        else if (particle_sym[i]->scale > level_info[j]->scale &&
                 particle_sym[i]->scale < level_info[j-1]->scale){
          flag = 1;
          current_level = j-1;
        }
        if (flag == 1) break;
      }
      level_info[current_level]->Update_cell_list(particle_sym[i], world);
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Map_the_sym_particle_to_tree finished\n";
  out.close();
#endif
}
#endif
#if PERI_DIM != 0
//--------------------------------------------------
// release pbc particle memory
//--------------------------------------------------
void SPH::Release_memory_of_pbc_paticle(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<concurrent_vector<p_Particle>::iterator>(particle_peri.begin(), particle_peri.end()),
    [&](const blocked_range<concurrent_vector<p_Particle>::iterator>& r){
    for(concurrent_vector<p_Particle>::iterator it=r.begin(); it!=r.end(); ++it){
      if (NULL != *it){
        delete (*it); // release memory of exchange_particles in heap
        *it = NULL;
      }
    }  
  }, ap);

  particle_peri.clear();
  particle_peri.shrink_to_fit();
  if (int(particle_peri.capacity()) != 0){
    cout<<"Ghost particle size is wrong!!!\n"; world.abort(-1);
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Release_memory_of_pbc_paticle finished\n";
  out<<"particle_peri.size "<<particle_peri.size()<<" capacity: "<<particlepool.capacity()<<" available nodes "<<particlepool.available_node()<<" \n";
  out.close();
#endif
}
//--------------------------------------------------
// construct periodical particles in buffer area
//--------------------------------------------------
void SPH::Construct_periodical_BC_particles(communicator &world, int flag)
{
  for( int i=0; i<num_level; i++){
    if (flag == 1){
      level_info[i]->Construct_periodical_BC_particles_1(particle_peri);
    }
    else if (flag == 2){
#ifdef _MPI_
      level_info[i]->Construct_periodical_BC_particles_2(particle_peri, world);
#endif
    }
  }

  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(particle_peri.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i)
      particle_peri[i]->id = i + glbl_total_num_particle;
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Construct_periodical_BC_particles finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Refresh periodical bc particle information
//--------------------------------------------------
void SPH::Refresh_pbc_particle_info(communicator &world, int flag)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(particle_peri.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      p_Particle cp = particle_peri[i]->copy[0];
      particle_peri[i]->Refresh_pbc_particle_info(cp, flag);
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Refresh_pbc_particle_info finished\n";
  out.close();
#endif
}
//--------------------------------------------------
//  map the peri particle to the tree according 
//  to the physical position
//--------------------------------------------------
void SPH::Map_the_peri_particle_to_tree(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(particle_peri.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int current_level = 0;
      int flag = 0;
      for( int j = num_level-1; j>=0; j--){
        if(particle_peri[i]->scale == level_info[j]->scale){
          flag = 1;
          current_level = j;
        }
        else if (particle_peri[i]->scale > level_info[j]->scale &&
                 particle_peri[i]->scale < level_info[j-1]->scale){
          flag = 1;
          current_level = j-1;
        }
        if (flag == 1) break;
      }
      level_info[current_level]->Update_cell_list(particle_peri[i], world);
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Map_the_peri_particle_to_tree finished\n";
  out.close();
#endif
}
#endif
#ifdef _MPI_
#if defined(_DIST_NUMBER_) || defined(_MASS_) || defined(_WEIGHTED_PARTITION_)
//--------------------------------------------------------
// get the mass for CVP partitioning
//-------------------------------------------------------
Real SPH::Get_total_p_mass_local(communicator &world){
  Real tmp = 0.;
  for(int i=0; i!=total_num_particle; ++i){
    tmp += particle[i]->p_mass;
  }
  return tmp;
}
#endif
#ifdef _WEIGHTED_PARTITION_
//--------------------------------------------------------
// get the mass for CVP partitioning
//-------------------------------------------------------
Real SPH::Get_global_neighbor_size(communicator &world, int flag){
  Real tmp = 0.;
  if (flag == 0){
    for(int i=0; i!=total_num_particle; ++i){
      tmp += particle[i]->mass;
    }
  }else if (flag == 1){
    for(int i=0; i!=total_num_particle; ++i){
      tmp += particle[i]->neighbor.size();
    }
  }
  Real tmp_global = 0.;

  all_reduce(world, tmp, tmp_global, std::plus<Real>());

  return tmp_global;
}
//--------------------------------------------------------
// get the mass for CVP partitioning
//-------------------------------------------------------
Real SPH::Get_global_dist_number(communicator &world){

  Real tmp = Get_total_p_mass_local(world);

  Real tmp_global = 0.;

  all_reduce(world, tmp, tmp_global, std::plus<Real>());

  return tmp_global;
}
//--------------------------------------------------------
// get the mass for CVP partitioning
//-------------------------------------------------------
void SPH::Get_weighted_patitioning_mass(int flag, communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle[i]->Get_weighted_patitioning_mass(this, flag, particle[i]->neighbor.size());
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Get_weighted_patitioning_mass finished\n";
  out.close();
#endif
}
#endif
//--------------------------------------------------------
// Normalize_particle_mass
//-------------------------------------------------------
void SPH::Normalize_particle_mass(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle[i]->Normalize_particle_mass(this);
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Normalize_particle_mass finished\n";
  out.close();
#endif
}
//--------------------------------------------------------
// get the mass for CVP partitioning
//-------------------------------------------------------
void SPH::Get_patitioning_mass(int flag, communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle[i]->Get_patitioning_mass(this, flag, particle[i]->neighbor.size());
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Get_patitioning_mass finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// update exchange information tag system
//--------------------------------------------------
void SPH::Update_exchange_infor(communicator &world)
{
  for( int i=0; i<num_level; i++)
    level_info[i]->Init_exchange_status();

  for( int i=1; i<num_level; i++)
    level_info[i]->Update_exchange_status(level_info, i-1);

  for( int i=0; i<num_level; i++)
    level_info[i]->Merge_exchange_status();
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_exchange_infor finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Release buffer particle memory
//--------------------------------------------------
void SPH::Release_memory_of_buffer_paticle(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<concurrent_vector<p_Particle>::iterator>(particle_buffer.begin(), particle_buffer.end()),
    [&](const blocked_range<concurrent_vector<p_Particle>::iterator>& r){
    for(concurrent_vector<p_Particle>::iterator it=r.begin(); it!=r.end(); ++it){
      if (NULL != *it){
        delete (*it); // release memory of exchange_particles in heap
        *it = NULL;
      }
    }
  }, ap);
  particle_buffer.clear();
  particle_buffer.shrink_to_fit();
  if (int(particle_buffer.capacity()) != 0){
    cout<<"Buffer particle size is wrong!!!\n"; world.abort(-1);
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Release_memory_of_buffer_paticle finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Release ghost particle memory
//--------------------------------------------------
void SPH::Release_memory_of_ghost_paticle(communicator &world)
{
  particle_ghost.clear();
  concurrent_vector <p_Particle>(particle_ghost).swap(particle_ghost);
  particle_ghost.shrink_to_fit();
  if (int(particle_ghost.capacity()) != 0){
    cout<<"Ghost particle size is wrong!!!\n"; world.abort(-1);
  }

  particle_ghost_start.clear();
  concurrent_vector <int>(particle_ghost_start).swap(particle_ghost_start);
  particle_ghost_start.shrink_to_fit();
  if (int(particle_ghost_start.capacity()) != 0){
    cout<<"Particle_ghost_start size is wrong!!!\n"; world.abort(-1);
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Release_memory_of_ghost_paticle finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// construct graph
//--------------------------------------------------
void SPH::Construct_graph_modified(communicator &world)
{
  if (world.rank() == 0){
    for( int i=0; i< num_level; i++)
      level_info[i]->Reset_color_list();
  }
#ifdef _NARROW_BAND_GRAPH_
  for( int i=0; i< num_level; i++)
    level_info[i]->Reset_narrow_band_list(world);
#endif
  // construct the graph
//  sph_graph.Update_graph_SPH(this, world);
  sph_graph.Update_graph_SPH_modified(this, world);

  // edge coloring
  sph_graph.Handle_graph_edge_coloring(world);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Construct_graph finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// construct graph
//--------------------------------------------------
void SPH::Construct_graph(communicator &world)
{
  if (world.rank() == 0){
    for( int i=0; i< num_level; i++)
      level_info[i]->Reset_color_list();
  }

  // construct the graph
  sph_graph.Update_graph_SPH(this, world);

  // edge coloring
  sph_graph.Handle_graph_edge_coloring(world);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Construct_graph finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// construct graph
//--------------------------------------------------
void SPH::Data_communication(int flag, communicator &world)
{
  serialization_vector <std::pair<int,int>> edge_pool;
  static affinity_partitioner ap;

  num_buffer_particle = 0;

  if (flag == 0)
    particle_ghost_start.push_back(0);

  if(world.rank()==0){
    // begin communication substeps
    for(size_t iter = 0; iter < sph_graph.graph_iteration; iter++){
      // reset the container
      edge_pool.Vector.clear();
      edge_pool.mem_size = 0;
      edge_pool.tag = 0;

      //select edege for communication
      static affinity_partitioner ap;
      parallel_for( blocked_range<int>(0, sph_graph.total_edges),
               [&](const blocked_range<int>& r){
        for(int i=r.begin(); i!=r.end(); ++i){
          if(sph_graph.edge_color[i] == iter)
            edge_pool.Vector.push_back(sph_graph.edge_pool_colored[i]);
        }
      }, ap);

      edge_pool.mem_size = int(edge_pool.Vector.size());
      edge_pool.tag = 1;

      // broadcast the current communication relationship
      broadcast(world,edge_pool,0);

      // findout target processor infor
      target_processor = -1;
      parallel_for( blocked_range<int>(0, edge_pool.Vector.size()),
               [&](const blocked_range<int>& r){
        for(int i=r.begin(); i!=r.end(); ++i){
          if(edge_pool.Vector[i].first == world.rank())
            target_processor = edge_pool.Vector[i].second;
          else if(edge_pool.Vector[i].second == world.rank())
            target_processor = edge_pool.Vector[i].first;
        }
      }, ap);

      // check for consistency
      if(world.rank() == target_processor)
        cout<<world.rank()<<" "<<"can never communicate with itself !!!"<<endl;

      // begin communication
      Communication_between_processor_pair(flag, iter, world);
    }
  }else{

    // begin communication substeps
    for(int iter = 0; iter < sph_graph.graph_iteration; iter++){
      // reset the container
      edge_pool.Vector.clear();
      edge_pool.mem_size = 0;
      edge_pool.tag = 0;
      
      // broadcast the current communication relationship
      broadcast(world,edge_pool,0);

      // findout target processor infor
      target_processor = -1;
      parallel_for( blocked_range<int>(0, edge_pool.Vector.size()),
               [&](const blocked_range<int>& r){
        for(int i=r.begin(); i!=r.end(); ++i){
          if(edge_pool.Vector[i].first == world.rank())
            target_processor = edge_pool.Vector[i].second;
          else if(edge_pool.Vector[i].second == world.rank())
            target_processor = edge_pool.Vector[i].first;
        }
      }, ap);

      // check for consistency
      if(world.rank() == target_processor)
        cout<<world.rank()<<" "<<"can never communicate with itself !!!"<<endl;

      // begin communication
      Communication_between_processor_pair(flag, iter, world);
    }
  }
  edge_pool.Vector.clear();
  if (icount == 0){
    exchange_local_old = Real(particle_buffer.size())/Real(total_num_particle);
    icount = 1;
  }

  if (particle_ghost_start.size() != (sph_graph.graph_iteration + 1)){
    cout<<"Num of graph iteration is wrong!!"<<particle_ghost_start.size()<<" "<<sph_graph.graph_iteration + 1<<"\n";
    world.abort(-1);
  }

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Data communication finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// communication between processor pair
//--------------------------------------------------
void SPH::Communication_between_processor_pair(int flag, int iter, communicator &world)
{
  // flag = 0: first data communication
  // flag = 1: second data communication
  // flag = 2: third data communication

  // exchange the topology information and construct buffer particles
  // for input
  serialization_vector <p_Particle> exchange_buffer_particles_in;
  exchange_buffer_particles_in.tag = 0;
  exchange_buffer_particles_in.Vector.clear();
  // for output
  serialization_vector <p_Particle> exchange_buffer_particles_out;
  exchange_buffer_particles_out.Vector.clear();

  // prepare the data for communication
  time_exchangebuffer.restart();

  if (flag == 0){
    for( int i=0; i<num_level; i++){
      level_info[i]->Exchange_buffer_between_pairs(world, target_processor, particle_ghost, flag);
    }
    int num_ghost = int(particle_ghost.size());
    particle_ghost_start.push_back(num_ghost);
  }
  int start = particle_ghost_start[iter];
  int end   = particle_ghost_start[iter+1];
  for (int i = start; i < end; i++){
    p_Particle current_particle = particle_ghost[i];
    current_particle->tag = flag;
    exchange_buffer_particles_out.Vector.push_back(current_particle);
  }
  exchange_buffer_particles_out.mem_size = int(exchange_buffer_particles_out.Vector.size());
  exchange_buffer_particles_out.tag = 1;

  if (flag == 0)
    time_for_exchangebuffer1[world.rank()] += time_exchangebuffer.elapsed();
  else
    time_for_exchangebuffer2[world.rank()] += time_exchangebuffer.elapsed(); 

  if(target_processor != -1){

    //nonblocking communication
    mpi::request reqs[2];
    reqs[0] = world.isend(target_processor, 0, exchange_buffer_particles_out);
    reqs[1] = world.irecv(target_processor, 0, exchange_buffer_particles_in);
    mpi::wait_all(reqs, reqs + 2);

    time_communication_flag.restart();

    if (flag == 0){
      static affinity_partitioner ap;
      parallel_for( blocked_range<int>(0, int(exchange_buffer_particles_in.Vector.size())),
                                   [&](const blocked_range<int>& r){
        for(int i=r.begin(); i!=r.end(); ++i){
          exchange_buffer_particles_in.Vector[i]->Reset_color(target_processor);
#if PERI_DIM != 0
          exchange_buffer_particles_in.Vector[i]->Shift_coordinate(domain, box_l, box_r);
#endif
        }
      }, ap);
    }

    for(int i=0; i!=int(exchange_buffer_particles_in.Vector.size()); ++i){
      if (flag == 0){
        particle_buffer.push_back(exchange_buffer_particles_in.Vector[i]);
      }
      else if (flag == 2 || flag == 1 || flag == 4){
        p_Particle temp = exchange_buffer_particles_in.Vector[i];
        if (temp->id != particle_buffer[num_buffer_particle]->id){
          cout<<"The sequence of received buffer particle is WRONG on cpu "<<world.rank()<<" ";
          cout<<temp->id<<" "<<temp->color<<" "<<particle_buffer[num_buffer_particle]->id<<" "<<particle_buffer[num_buffer_particle]->color<<" "<<"\n";
          world.abort(-1);
        }
        particle_buffer[num_buffer_particle]->Set_buffer_particle_info(temp, flag);
      }
      num_buffer_particle++;
    }

    if (flag == 0)
      time_for_communication_flag1[world.rank()] += time_communication_flag.elapsed();
    else if (flag == 2 || flag == 1 || flag == 4)
      time_for_communication_flag2[world.rank()] += time_communication_flag.elapsed();

  }
  exchange_buffer_particles_out.Vector.clear();
  concurrent_vector <p_Particle>(exchange_buffer_particles_out.Vector).swap(exchange_buffer_particles_out.Vector);
  exchange_buffer_particles_out.Vector.shrink_to_fit();

  // if the buffer memory are used, just distroy pointers;
  // else distroy the memory as well
  if (flag == 0){
    exchange_buffer_particles_in.Vector.clear();
    concurrent_vector <p_Particle>(exchange_buffer_particles_in.Vector).swap(exchange_buffer_particles_in.Vector);
    exchange_buffer_particles_in.Vector.shrink_to_fit();
  }
  else if (flag == 2 || flag == 1 || flag == 4){
    static affinity_partitioner ap;
    parallel_for( blocked_range<concurrent_vector<p_Particle>::iterator>(exchange_buffer_particles_in.Vector.begin(), exchange_buffer_particles_in.Vector.end()),
        [&](const blocked_range<concurrent_vector<p_Particle>::iterator>& r){
      for(concurrent_vector<p_Particle>::iterator it=r.begin(); it!=r.end(); ++it){
        if (NULL != *it){
          delete (*it); // release memory of exchange_particles in heap
          *it = NULL;
        }
      }  
    }, ap);
    exchange_buffer_particles_in.Vector.clear();
    exchange_buffer_particles_in.Vector.shrink_to_fit();
    if (int(exchange_buffer_particles_in.Vector.capacity()) != 0){
      cout<<"Vector memory is not released!!!\n"; world.abort(-1);
    }
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Communication_between_processor_pair "<<world.rank()<<" & "<<target_processor<<" finished\n";
  out.close();
#endif
}
//--------------------------------------------------
//  collect partition mass info
//--------------------------------------------------
void SPH::Calculate_total_partitioning_mass(communicator &world)
{
  Get_patitioning_mass(1, world);
#ifdef _WEIGHTED_PARTITION_
  glbl_total_dist_number   = Get_global_dist_number(world);
  glbl_total_neighbor_size = Get_global_neighbor_size(world, 1);

  Real t_force_avg1 = 0.;
  Real t_force_avg2 = 0.;
  Real t_neighbor_avg = 0.;
  all_reduce(world, time_for_refresh_neighbor[world.rank()], t_neighbor_avg,  std::plus<Real>());
  all_reduce(world, time_for_force_calculation[world.rank()], t_force_avg1,  std::plus<Real>());
  all_reduce(world, time_for_density[world.rank()], t_force_avg2,  std::plus<Real>());
  t_neighbor_avg /= world.size();
  t_force_avg1 /= world.size();
  t_force_avg2 /= world.size();

  t_force            = t_force_avg1 + t_force_avg2 - t_force_old;
  t_neighbor         = t_neighbor_avg - t_neighbor_old;

  t_force_old        = t_force_avg1 + t_force_avg2;
  t_neighbor_old     = t_neighbor_avg;

  weight = AMAX1(AMIN1(t_force/(t_force + t_neighbor + 1.e-20), 1.), 0.);

  Get_weighted_patitioning_mass(1, world);
#endif
  total_mass = Get_total_p_mass_local(world);

  glbl_total_mass = 0.;
  all_reduce(world, total_mass, glbl_total_mass, std::plus<Real>());
  total_mass /= glbl_total_mass;
  Normalize_particle_mass(world);
  glbl_total_mass = 1.;
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Calculate_total_partitioning_mass "<<world.rank()<<" & "<<target_processor<<" finished\n";
  out.close();
#endif
}
//--------------------------------------------------
//  check whether partition is needed
//--------------------------------------------------
void SPH::Check_need_for_partition(communicator &world)
{
  exchange_local = Real(particle_buffer.size())/Real(total_num_particle);

  Calculate_total_partitioning_mass(world);

  if (need_for_partition != 1){
    int local_need_for_partition = 0;

    error_local = exchange_local-exchange_local_old;

    imbalance_local = (total_mass - total_mass_old)/total_mass_old;

    if (error_local >= error_tolerance || imbalance_local >= imbalance_tolerance){
      local_need_for_partition = 1;
    }
    all_reduce(world, error_local, error_global, mpi::maximum<Real>());
    all_reduce(world, imbalance_local, imbalance_global, mpi::maximum<Real>());
    all_reduce(world, local_need_for_partition, need_for_partition, mpi::maximum<int>());
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Check_need_for_partition finished\n";
  out.close();
#endif
}
//--------------------------------------------------
//  map the buffer particle to the tree according 
//  to the physical position
//--------------------------------------------------
void SPH::Map_the_buffer_particle_to_tree(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(particle_buffer.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int current_level = 0;
      int flag = 0;
      for( int j = num_level-1; j>=0; j--){
        if(particle_buffer[i]->scale == level_info[j]->scale){
          flag = 1;
          current_level = j;
        }
        else if (particle_buffer[i]->scale > level_info[j]->scale &&
                 particle_buffer[i]->scale < level_info[j-1]->scale){
          flag = 1;
          current_level = j-1;
        }
        if (flag == 1) break;
      }
#if PERI_DIM != 0
      level_info[current_level]->Update_cell_list_and_shift_coord(particle_buffer[i], world);
#else
      level_info[current_level]->Update_cell_list(particle_buffer[i], world);
#endif
    }
  }, ap);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Map_the_buffer_particle_to_tree finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// find the out bound only in the coarsest level
// in each dimension
//--------------------------------------------------
void SPH::Find_bound_in_coarsest_level(communicator &world)
{
  my_set_const (local_box_l, 0.);
  my_set_const (local_box_r, 0.);

  ParallelGetLocalBound (total_num_particle, particle, local_box_l, local_box_r);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Find_bound_in_coarsest_level finished\n";
  out.close();
#endif
}
#if PERI_DIM != 0
//--------------------------------------------------
// shift particle positions if all the particles are
// out of domain
//--------------------------------------------------
void SPH::Check_PBC_particle_positions(communicator &world)
{
  int flag = 0;
  my_int dim; dim.i = dim.j = dim.k = 0;
  if ((local_box_l.i > box_r.i + 1.e-20 && DIM_X == 1)){
    flag  = 1;
    dim.i = -1;
  }else if((local_box_r.i < box_l.i - 1.e-20 && DIM_X == 1)){
    flag  = 1;
    dim.i = 1;
  }
  if ((local_box_l.j > box_r.j + 1.e-20 && DIM_Y == 1)){
    flag  = 1;
    dim.j = -1;
  }else if ((local_box_r.j < box_l.j - 1.e-20 && DIM_Y == 1)){
    flag  = 1;
    dim.j = 1;
  }
  if ((local_box_l.k > box_r.k + 1.e-20 && DIM_Z == 1)){
    flag  = 1;
    dim.k = -1;
  }else if((local_box_r.k < box_l.k - 1.e-20 && DIM_Z == 1)){
    flag  = 1;
    dim.k = 1;
  }

  if (flag != 0){
    static affinity_partitioner ap;
    parallel_for( blocked_range<int>(0, int(total_num_particle)),[&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i)
        particle[i]->Shift_coordinate_constrained(domain, box_l, box_r, dim);
    }, ap);
//    sph_voronoi.Shift_coordinate_constrained(dim, this, world);
    Find_bound_in_coarsest_level (world);
  }

//  sph_voronoi.Sync_VP_positions(this, world);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Check_PBC_particle_positions finished\n";
  out.close();
#endif
}
#endif
//--------------------------------------------------
// Build the local tree top-to-down
//--------------------------------------------------
void SPH::Traverse_and_construct_local_tree(communicator &world, int flag)
{
  for( int i=0; i< num_level; i++){
    if (flag == 1)
      level_info[i]->Initialize( i+Lmin, this, world);
    else if (flag == 2) 
      level_info[i]->Reinitialize_cell_list(world, this);
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Traverse_and_construct_local_tree finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Build dynamic essential local tree
// Currently only roughly local from the coarsest
// level
//--------------------------------------------------
void SPH::Build_local_tree(communicator &world, int flag)
{
  // flag = 1: initialize
  // flag = 2: reinitialize

  Find_bound_in_coarsest_level (world);
#if PERI_DIM != 0
  Check_PBC_particle_positions (world);
#endif
  Traverse_and_construct_local_tree (world, flag);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Build_local_tree finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// migrate SPH particles according to the graph
//--------------------------------------------------
void SPH::Migrate_SPH_particle(communicator &world){
  serialization_vector <std::pair<int,int>> edge_pool;
  static affinity_partitioner ap;

  // store incoming particles
  concurrent_vector <p_Particle> incoming_particles;
  incoming_particles.clear();
  if(world.rank()==0){
    // begin communication substeps
    for(size_t iter = 0; iter < sph_voronoi.vo_graph.graph_iteration; iter++){
      // reset the container
      edge_pool.Vector.clear();
      edge_pool.mem_size = 0;
      edge_pool.tag = 0;

      //select edege for communication
      parallel_for( blocked_range<int>(0, sph_voronoi.vo_graph.total_edges),
               [&](const blocked_range<int>& r){
        for(int i=r.begin(); i!=r.end(); ++i){
          if(sph_voronoi.vo_graph.edge_color[i] == iter)
            edge_pool.Vector.push_back(sph_voronoi.vo_graph.edge_pool_colored[i]);
        }
      }, ap);

      edge_pool.mem_size = int(edge_pool.Vector.size());
      edge_pool.tag = 1;

      // broadcast the current communication relationship
      broadcast(world,edge_pool,0);

      // findout target processor infor
      target_processor = -1;
      parallel_for( blocked_range<int>(0, edge_pool.Vector.size()),
               [&](const blocked_range<int>& r){
        for(int i=r.begin(); i!=r.end(); ++i){
          if(edge_pool.Vector[i].first == world.rank())
            target_processor = edge_pool.Vector[i].second;
          else if(edge_pool.Vector[i].second == world.rank())
            target_processor = edge_pool.Vector[i].first;
        }
      }, ap);

      // check for consistency
      if(world.rank() == target_processor)
        cout<<world.rank()<<" "<<"can never communicate with itself !!!"<<endl;

      // begin migration
      Migration_between_processor_pair(iter, incoming_particles, world);
    }
  }else{
    // begin communication substeps
    for(int iter = 0; iter < sph_voronoi.vo_graph.graph_iteration; iter++){
      // reset the container
      edge_pool.Vector.clear();
      edge_pool.mem_size = 0;
      edge_pool.tag = 0;
      
      // broadcast the current communication relationship
      broadcast(world,edge_pool,0);

      // findout target processor infor
      target_processor = -1;
      parallel_for( blocked_range<int>(0, edge_pool.Vector.size()),
               [&](const blocked_range<int>& r){
        for(int i=r.begin(); i!=r.end(); ++i){
          if(edge_pool.Vector[i].first == world.rank())
            target_processor = edge_pool.Vector[i].second;
          else if(edge_pool.Vector[i].second == world.rank())
            target_processor = edge_pool.Vector[i].first;
        }
      }, ap);

      // check for consistency
      if(world.rank() == target_processor)
        cout<<world.rank()<<" "<<"can never communicate with itself !!!"<<endl;

      // begin communication
      Migration_between_processor_pair(iter, incoming_particles, world);
    }
  }

  migrate_local = Real(total_num_particle);
  edge_pool.Vector.clear();
  concurrent_vector <unsigned int> record;
  record.clear();
  for (int i = 0; i < total_num_particle; i++){
    if (particle[i]->color == world.rank()){
      record.push_back(i);
    }
  }
  for (int i = 0; i < record.size(); i++){
      p_Particle temp = particle[record[i]];
      particle[record[i]] = particle[i];
      particle[i] = temp;
  }
  for (int i = record.size(); i < total_num_particle; i++){
      particle[i]->Clearup();
      particlepool.free(particle[i]);
  }
  total_num_particle = record.size() + incoming_particles.size();
  migrate_local = fabs(migrate_local-Real(record.size()))/Real(total_num_particle);
  p_Particle *tmp;
  tmp = new p_Particle[total_num_particle];
  parallel_for( blocked_range<int>(0, record.size()),
           [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      tmp[i] = particle[i];
    }
  }, ap);
  parallel_for( blocked_range<int>(0, incoming_particles.size()),
           [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      tmp[i+record.size()] = particlepool.malloc();
      tmp[i+record.size()]->Set_particle_info(incoming_particles[i]);
    }
  }, ap);
  parallel_for( blocked_range<concurrent_vector<p_Particle>::iterator>(incoming_particles.begin(), incoming_particles.end()),
    [&](const blocked_range<concurrent_vector<p_Particle>::iterator>& r){
    for(concurrent_vector<p_Particle>::iterator it=r.begin(); it!=r.end(); ++it){
      if (NULL != *it){
        delete (*it); // release memory of exchange_particles in heap
        *it = NULL;
      }
    }  
  }, ap);
  incoming_particles.clear();
  incoming_particles.shrink_to_fit();
  if (int(incoming_particles.capacity()) != 0){
    cout<<"Vector memory is not released!!!\n"; world.abort(-1);
  }
  delete [] particle;
  particle = new p_Particle[total_num_particle];
  parallel_for( blocked_range<int>(0, total_num_particle),
           [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle[i] = tmp[i];
    }
  }, ap);
  delete [] tmp;

  int glbl_total_num_particle_temp;
  reduce(world, total_num_particle, glbl_total_num_particle_temp, std::plus<int>(), 0);

  if (world.rank() == 0 && glbl_total_num_particle_temp != glbl_total_num_particle){
    cout<<"Total number of particle is wrong!!"<<glbl_total_num_particle_temp<<" "<<glbl_total_num_particle<<"\n";
  }

  record.clear();
  concurrent_vector <unsigned int>(record).swap(record);
  record.shrink_to_fit();
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Migrate_SPH_particle finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// particle migration between processor pair
//--------------------------------------------------
void SPH::Migration_between_processor_pair(int iter, concurrent_vector <p_Particle> &incoming_particles, communicator &world)
{
  // exchange the topology information and construct migrating particles
  // for output
  serialization_vector <p_Particle> exchange_migrate_particles_in;
  exchange_migrate_particles_in.Vector.clear();
  // for input
  serialization_vector <p_Particle> exchange_migrate_particles_out;
  exchange_migrate_particles_out.Vector.clear();

  // prepare the data for communication

  if(target_processor != -1){

    static affinity_partitioner ap;
    parallel_for( blocked_range<int>(0, total_num_particle),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        if (particle[i]->color == target_processor){
          particle[i]->tag = 3;
          exchange_migrate_particles_out.Vector.push_back(particle[i]);
        }
      }
    }, ap);
  
    exchange_migrate_particles_out.mem_size = int(exchange_migrate_particles_out.Vector.size());
    exchange_migrate_particles_out.tag = 1;

    //nonblocking communication
    mpi::request reqs[2];
    reqs[0] = world.isend(target_processor, 0, exchange_migrate_particles_out);
    reqs[1] = world.irecv(target_processor, 0, exchange_migrate_particles_in);
    mpi::wait_all(reqs, reqs + 2);

    parallel_for( blocked_range<int>(0, int(exchange_migrate_particles_in.Vector.size())),
                                 [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        if (exchange_migrate_particles_in.Vector[i]->color != world.rank()){
          cout<<"Color is wrong for migrated particle!!\n"; world.abort(-1);
        }
      	incoming_particles.push_back(exchange_migrate_particles_in.Vector[i]);
      }
    }, ap);
  }

  exchange_migrate_particles_in.Vector.clear();
  concurrent_vector <p_Particle>(exchange_migrate_particles_in.Vector).swap(exchange_migrate_particles_in.Vector);
  exchange_migrate_particles_in.Vector.shrink_to_fit();
  if (int(exchange_migrate_particles_in.Vector.capacity()) != 0){
    cout<<"Vector memory is not released!!!\n"; world.abort(-1);
  }
  exchange_migrate_particles_out.Vector.clear();
  concurrent_vector <p_Particle>(exchange_migrate_particles_out.Vector).swap(exchange_migrate_particles_out.Vector);
  exchange_migrate_particles_out.Vector.shrink_to_fit();
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Migration_between_processor_pair "<<world.rank()<<" & "<<target_processor<<" finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// CVP
//--------------------------------------------------
void SPH::Partitioning(communicator &world){

    time_release_particle.restart();
  // Release buffer particle memory
  Release_memory_of_buffer_paticle(world);
#if PERI_DIM != 0
    // Release buffer particle memory
    Release_memory_of_pbc_paticle(world);
#endif
#if SYM_DIM != 0
    // Release buffer particle memory
    Release_memory_of_sbc_paticle(world);
#endif
    time_for_release_particle[world.rank()] += time_release_particle.elapsed();

#if PERI_DIM != 0
  // shift back SPH particle position
  sph_voronoi.Shift_back_SPH_particle_position(this, world);
#endif

  // CVP partitioning method
  sph_voronoi.Partitioning(this, world);

    time_communication_partition.restart();
  // migrate SPH particles according to the graph
  Migrate_SPH_particle(world);

#if PERI_DIM != 0
  // shift SPH particle position according to
  // the new VD
  sph_voronoi.Shift_SPH_particle_position(this, world);
#endif

#ifdef _MEM_CHECK_
  Check_memory_consumption(world);
#endif
  icount                      = 0;
  need_for_rebuild_graph      = 0;
  need_for_rebuild_local_tree = 0;
  need_for_partition_count    = 0;
  exchange_local_old          = 0;
  error_local                 = 0.;
  error_global                = 0.;

  total_mass_old = Get_total_p_mass_local(world);

#if defined(_WEIGHTED_PARTITION_)

  Real t_force_avg1 = 0.;
  Real t_force_avg2 = 0.;
  Real t_neighbor_avg = 0.;
  all_reduce(world, time_for_refresh_neighbor[world.rank()], t_neighbor_avg,  std::plus<Real>());
  all_reduce(world, time_for_force_calculation[world.rank()], t_force_avg1,  std::plus<Real>());
  all_reduce(world, time_for_density[world.rank()], t_force_avg2,  std::plus<Real>());
  t_neighbor_avg /= world.size();
  t_force_avg1 /= world.size();
  t_force_avg2 /= world.size();

  t_force_old        = t_force_avg1 + t_force_avg2;
  t_neighbor_old     = t_neighbor_avg;

#endif
  total_mass                  = 0.;
  imbalance_local             = 0.;
  imbalance_global            = 0.;
  time_for_communication_partition[world.rank()] += time_communication_partition.elapsed();
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Partitioning finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// update voronoi particles
//--------------------------------------------------
void SPH::Update_vp_position_mean_velo(communicator &world){

  my_set_const (v_avg, 0.0);

  ParallelGetMeanVelocity(total_num_particle, particle, v_avg);
  
  // update voronoi particle position according to 
  // fluid field information
  sph_voronoi.Update_vp_position_mean_velo(this, world);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_vp_position finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// update voronoi particles
//--------------------------------------------------
void SPH::Update_vp_position_mass_center(communicator &world){

  my_set_const( mass_center, 0.);
  my_set_const( glbl_mass_center, 0.);

  ParallelGetMassCenter(total_num_particle, total_mass, particle, mass_center);

  // update voronoi particle position according to 
  // fluid field information
  sph_voronoi.Update_vp_position_mass_center(this, world);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_vp_position finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Output global information
//--------------------------------------------------
void SPH::Output_global_info(int flag, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s","./outdata/global_info.csv");
  ofstream out(filename, ios::app);

  if (flag == 0){
    if (world.rank() == 0){
      out<<"\"run_time\",\"t_graph\",\"t_exchange_graph_topology\",\"t_construct_graph\",\"t_edge_coloring\",\"t_partition\",\"t_communication0\",\"t_communication1\",\"t_communication2\",\"t_communication_partition\",\"t_simulation\",\"t_mapping\",\"t_build_local_map\",\"t_map_to_tree\",\"t_update_map\",\"t_clear_list\",\"t_release_particle\",\"t_reset_neighbor\",\"t_refresh_neighbor\",\"t_force\",\"t_density\",\"t_rest\",\"t_update\",\"t_total\",\"particle_exchange_total\",\"particle_exchange_max\",\"particle_exchange_min\",\"particle_exchange_avg\",\"particle_exchange_error_global\",\"particle_imbalance_total\",\"particle_imbalance_max\",\"particle_imbalance_min\",\"particle_imbalance_avg\",\"particle_imbalance_error_global\",\"particle_migrated_total\",\"particle_migrated_max\",\"particle_migrated_min\",\"particle_migrated_avg\",\"global_force_error\",\"total_edges\",\"maximal_degree\",\"num_colors\",\"num_proc\"\n";
    }
  }else{
    if (world.rank() == 0){
      reduce( world, exchange_local, exchange_total, std::plus<Real>(), 0);
      reduce( world, exchange_local, exchange_max, mpi::maximum<Real>(), 0);
      reduce( world, exchange_local, exchange_min, mpi::minimum<Real>(), 0);
      exchange_avg = exchange_total/glbl_total_num_color;

      reduce( world, imbalance_local, imbalance_total, std::plus<Real>(), 0);
      reduce( world, imbalance_local, imbalance_max, mpi::maximum<Real>(), 0);
      reduce( world, imbalance_local, imbalance_min, mpi::minimum<Real>(), 0);
      imbalance_avg = imbalance_total/glbl_total_num_color;

      reduce( world, migrate_local, migrate_total, std::plus<Real>(), 0);
      reduce( world, migrate_local, migrate_max, mpi::maximum<Real>(), 0);
      reduce( world, migrate_local, migrate_min, mpi::minimum<Real>(), 0);
      migrate_avg = migrate_total/glbl_total_num_color;

      reduce( world, time_for_graph[world.rank()], time_for_graph[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, sph_graph.time_for_exchange_graph_topology[world.rank()], sph_graph.time_for_exchange_graph_topology[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, sph_graph.time_for_construct_graph[world.rank()], sph_graph.time_for_construct_graph[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, sph_graph.time_for_edge_coloring[world.rank()], sph_graph.time_for_edge_coloring[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_partition[world.rank()], time_for_partition[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_communication0[world.rank()], time_for_communication0[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_communication1[world.rank()], time_for_communication1[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_communication2[world.rank()], time_for_communication2[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_communication_partition[world.rank()], time_for_communication_partition[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_simulation[world.rank()], time_for_simulation[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_mapping[world.rank()], time_for_mapping[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_bulid_local_map[world.rank()], time_for_bulid_local_map[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_map_particle_to_tree[world.rank()], time_for_map_particle_to_tree[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_update_every_level_info[world.rank()], time_for_update_every_level_info[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_clear_cell_list[world.rank()], time_for_clear_cell_list[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_release_particle[world.rank()], time_for_release_particle[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_reset_neighbor[world.rank()], time_for_reset_neighbor[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_refresh_neighbor[world.rank()], time_for_refresh_neighbor[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_force_calculation[world.rank()], time_for_force_calculation[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_density[world.rank()], time_for_density[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_rest[world.rank()], time_for_rest[world.size()],  mpi::maximum<Real>(), 0);
      reduce( world, time_for_update[world.rank()], time_for_update[world.size()],  mpi::maximum<Real>(), 0);

      out<<run_time<<","<<time_for_graph[world.size()]<<","
         <<sph_graph.time_for_exchange_graph_topology[world.size()]<<","
         <<sph_graph.time_for_construct_graph[world.size()]<<","
         <<sph_graph.time_for_edge_coloring[world.size()]<<","
         <<time_for_partition[world.size()]<<","
         <<time_for_communication0[world.size()]<<","
         <<time_for_communication1[world.size()]<<","
         <<time_for_communication2[world.size()]<<","
         <<time_for_communication_partition[world.size()]<<","
         <<time_for_simulation[world.size()]<<","
         <<time_for_mapping[world.size()]<<","
         <<time_for_bulid_local_map[world.size()]<<","
         <<time_for_map_particle_to_tree[world.size()]<<","
         <<time_for_update_every_level_info[world.size()]<<","
         <<time_for_clear_cell_list[world.size()]<<","
         <<time_for_release_particle[world.size()]<<","
         <<time_for_reset_neighbor[world.size()]<<","
         <<time_for_refresh_neighbor[world.size()]<<","
         <<time_for_force_calculation[world.size()]<<","
         <<time_for_density[world.size()]<<","
         <<time_for_rest[world.size()]<<","
         <<time_for_update[world.size()]<<","
         <<time_for_total[world.size()]<<","
         <<exchange_total<<","
         <<exchange_max<<","
         <<exchange_min<<","
         <<exchange_avg<<","
         <<error_global<<","
         <<imbalance_total<<","
         <<imbalance_max<<","
         <<imbalance_min<<","
         <<imbalance_avg<<","
         <<imbalance_global<<","
         <<migrate_total<<","
         <<migrate_max<<","
         <<migrate_min<<","
         <<migrate_avg<<","
         <<glbl_force_error<<","
         <<sph_graph.total_edges<<","
         <<sph_graph.maximal_degree<<","
         <<sph_graph.graph_iteration<<","
         <<glbl_total_num_color<<"\n";
    }else{
      reduce( world, exchange_local, std::plus<Real>(), 0);
      reduce( world, exchange_local, mpi::maximum<Real>(), 0);
      reduce( world, exchange_local, mpi::minimum<Real>(), 0);
  
      reduce( world, imbalance_local, std::plus<Real>(), 0);
      reduce( world, imbalance_local, mpi::maximum<Real>(), 0);
      reduce( world, imbalance_local, mpi::minimum<Real>(), 0);

      reduce( world, migrate_local, std::plus<Real>(), 0);
      reduce( world, migrate_local, mpi::maximum<Real>(), 0);
      reduce( world, migrate_local, mpi::minimum<Real>(), 0);

      reduce( world, time_for_graph[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, sph_graph.time_for_exchange_graph_topology[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, sph_graph.time_for_construct_graph[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, sph_graph.time_for_edge_coloring[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_partition[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_communication0[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_communication1[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_communication2[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_communication_partition[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_simulation[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_mapping[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_bulid_local_map[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_map_particle_to_tree[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_update_every_level_info[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_clear_cell_list[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_release_particle[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_reset_neighbor[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_refresh_neighbor[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_force_calculation[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_density[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_rest[world.rank()], mpi::maximum<Real>(), 0);
      reduce( world, time_for_update[world.rank()], mpi::maximum<Real>(), 0);
    }
  }
}
//--------------------------------------------------
// Output color list infomation
//--------------------------------------------------
void SPH::Output_color_list_dat(int n, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s%d%s%d%s","./outdata/simulation_color_list.",n,".",world.rank(),".dat");

  if (world.rank() == 0){
    ofstream out(filename, ios::trunc);
    for (int i = 0; i < num_level; i++){
      level_info[i]->Output_color_list(world, filename, n);
    }
    out.close();
  }
}
#endif
#ifdef _MEM_CHECK_
//--------------------------------------------------
// Check memory usage
//--------------------------------------------------
void SPH::Check_memory_consumption(communicator &world)
{
  local_mem      = 0;
  mem_particle   = 0;
  mem_cell       = 0;
  mem_list       = 0;
  glbl_mem_max   = 0;
  glbl_mem_total = 0;

  long long mem_particle_temp = 0;
  long long mem_cell_temp     = 0;
  long long mem_list_temp     = 0;

  long long count = total_num_particle;
  for (int i = 0; i < total_num_particle; i++){
    count += (long long)(particle[i]->neighbor.capacity() + 1);
  }
  mem_particle_temp += (long long)(particlepool.capacity()*sizeof(Particle));
#ifdef _MPI_
  for (int i = 0; i < particle_buffer.size(); i++){
    count += (long long)(particle_buffer[i]->neighbor.capacity() + 1);
  }
  mem_particle_temp += (long long)(particle_buffer.capacity()*sizeof(Particle));
  count += (long long)(particle_buffer.size()+particle_ghost.size());
#endif
#if SYM_DIM != 0
  for (int i = 0; i < particle_sym.size(); i++){
    count += (long long)(particle_sym[i]->neighbor.capacity() + 1);
  }
  mem_particle_temp += (long long)(particle_sym.capacity()*sizeof(Particle));
  count += (long long)(particle_sym.size());
#endif
#if PERI_DIM != 0
  for (int i = 0; i < particle_peri.size(); i++){
    count += (long long)(particle_peri[i]->neighbor.capacity() + 1);
  }
  mem_particle_temp += (long long)(particle_peri.capacity()*sizeof(Particle));
  count += (long long)(particle_peri.size());
#endif

  mem_cell_temp = (long long)(cell_listpool.capacity()*sizeof(Cell_list));
  for (int i = 0; i < num_level; i++){
    level_info[i]->Check_memory_consumption(world, mem_cell_temp, mem_list_temp);
  }

  mem_particle = Real((long double)(mem_particle_temp + count*sizeof(p_Particle))/1024./1024.);
      mem_cell = Real((long double)mem_cell_temp/1024./1024.);
      mem_list = Real((long double)mem_list_temp/1024./1024.);

  local_mem += (mem_particle + mem_cell + mem_list);

  if (world.rank() == 0){
    reduce(world, local_mem, glbl_mem_max, mpi::maximum<Real>(), 0);
    reduce(world, local_mem, glbl_mem_total, std::plus<int>(), 0);

    cout<<"<<<<<< total memory used(MB): "<<glbl_mem_total<<" | "<<" max memory used in single task(MB): "<<glbl_mem_max<<"\n";
  }else{
    reduce(world, local_mem, mpi::maximum<Real>(), 0);
    reduce(world, local_mem, std::plus<int>(), 0);
  }
}
#endif
//--------------------------------------------------
// Output pov ray infomation
//--------------------------------------------------
void SPH::Output_pov_ray_file(int n, communicator &world)
{
  if (world.rank() == 0){

    char    filename[256];
    sprintf(filename,"%s%d%s%d%s","./tpout/simulation_sp.",n,".",world.rank(),".pov");
    ofstream out(filename, ios::trunc);

    visual.Set_scene(filename, 2);

    for (int i = 0; i < total_num_particle; i++)
      visual.Draw_particle(filename, particle[i]->coord, particle[i]->color,
                           particle[i]->h/CELL_RATIO, 0., Real(glbl_total_num_color), 1.);
    out.close();
  }else{
  
    char    filename[256];
    sprintf(filename,"%s%d%s%d%s","./tpout/simulation_sp.",n,".",world.rank(),".pov");
    ofstream out(filename, ios::trunc);

    for (int i = 0; i < total_num_particle; i++)
      visual.Draw_particle(filename, particle[i]->coord, particle[i]->color,
                           particle[i]->h/CELL_RATIO/2., 0., Real(glbl_total_num_color), 1.);
    out.close();
  }
}
//--------------------------------------------------
// Output level infomation
//--------------------------------------------------
void SPH::Output_every_level_dat(int n, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s%d%s%d%s","./outdata/simulation_level.",n,".",world.rank(),".dat");

  ofstream out(filename, ios::trunc);
  for (int i = 0; i < num_level; i++){
    level_info[i]->Output_level(world, filename, n);
  }

  out.close();
}
//--------------------------------------------------
// Output particle infomation
//--------------------------------------------------
void SPH::Output_plt_file(int n, int flag, communicator &world)
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
  out<<"VARIABLES = \"x\",\"y\",\"z\",\"id\",\"p_mass\",\"color\"\n";
  for (int i = 0; i < total_num_particle; i++){
    Particle *current_particle = particle[i];
    my_real p_coord; my_set_const(p_coord, 0.);
    my_shift_coordinate(p_coord, current_particle->coord, domain, box_l, box_r);
    out<<p_coord.i<<" "
       <<p_coord.j<<" "
       <<p_coord.k<<" "
       <<current_particle->id<<" "
       <<current_particle->p_mass<<" "
       <<current_particle->color<<"\n";
  }
  out.close();
}
//--------------------------------------------------
// Output restart file
//--------------------------------------------------
void SPH::Output_restart_file()
{
}
//--------------------------------------------------
// Read initial particle infomation
//--------------------------------------------------
void SPH::Read_infile(){
  char    filename[256];
  sprintf(filename,"%s","coord.plt");

  ifstream in(filename, ios::in);

  for (int i = 0; i < total_num_particle; i++){
    Particle *current_particle = particle[i];
    in>>current_particle->coord.i>>current_particle->coord.j>>current_particle->coord.k;
  }
  in.close();
}
