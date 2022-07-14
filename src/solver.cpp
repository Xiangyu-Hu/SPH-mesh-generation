#include "solver.h"

/***************************************************/
/*                                                 */
/*       Functions defined in class "Solver"       */
/*                                                 */
/***************************************************/

//--------------------------------------------------
// Solver
//--------------------------------------------------
Solver::Solver(Initialization &Ini, communicator &world):
  current_solver (Ini, world)
{
  total_iterate_num = Ini.num_output;
  output_dt         = Ini.output_dt;
  timestep          = 0.0;
  integral_time     = 0.0;
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Class Solver is initialized\n";
    cout<<"**********************************************************\n";
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<class solver initialized\n";
  out.close();
#endif
}
//--------------------------------------------------
// choose the solver
//--------------------------------------------------
void Solver::Load_solver(Initialization &Ini, communicator &world){

#ifndef _EULARIAN_CAT_ // particle based solver

  #ifdef _MPI_
  current_solver.Output_global_info(0, world);
  #endif
  
  current_solver.Initialize_case(world);
  
  Real Time = Ini.Get_start_time();  //computation time
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Simulation starts at time: "<< Time<< "\n";
  }

  // current_solver.Output_plt_file(0, 1, world);
  // current_solver.Output_paraview_file(0, 1, world);
  current_solver.Post_processing (0, world);
//current_solver.Output_buffer_plt_file(0, world);
  #ifdef _MPI_
  current_solver.sph_voronoi.Output_plt_file(0, 0, world);
  current_solver.sph_voronoi.Output_vtk_file(0, 0, world);
  current_solver.sph_voronoi.Output_pov_ray_file(0,0, world);
  #endif
//  world.barrier();
//  world.abort(-1);
  #ifdef _MEM_CHECK_
  current_solver.Check_memory_consumption(world);
  #endif
  current_solver.time_total.restart();

  for(int n=1; n<=Ini.Get_num_timestep(); n++){
    Set_current_time(Time);

    if (world.rank() == 0) cout<<"Time:  "<<Time<<"\n";
  #ifndef _MESH_GENERATION_
    MRSPH_Solver(Ini, n, world);
  #else
    SPH_Mesh_Generation_Solver(Ini, n, world);
  #endif
    Time += Ini.Get_dt();     //run time for next time step
  }
  current_solver.time_for_total[world.rank()] = current_solver.time_total.elapsed();

  if (world.rank() == 0){
    reduce( world, current_solver.time_for_total[world.rank()], current_solver.time_for_total[world.size()],  mpi::maximum<Real>(), 0);
    int n = 0;
  #ifdef _MPI_
    n = world.size();
  #else
    n = 0;
  #endif
    cout<<"**********************************************************\n";
    cout<<"Time for Simulation:       "<<current_solver.time_for_simulation[n]<<"\n";
    cout<<"Time for Mapping:          "<<current_solver.time_for_mapping[n]<<"\n";
    cout<<"Time for dist_cal:         "<<current_solver.time_for_refresh_neighbor[n]<<"\n";
    cout<<"Time for density:          "<<current_solver.time_for_density[n]<<"\n";
    cout<<"Time for force:            "<<current_solver.time_for_force_calculation[n]<<"\n";
    cout<<"Time for partition:        "<<current_solver.time_for_partition[n]<<"\n";
    cout<<"Time for particle migrate: "<<current_solver.time_for_communication_partition[n]<<"\n";
    cout<<"Time for Communication0:   "<<current_solver.time_for_communication0[n]<<"\n";
    cout<<"Time for Communication1:   "<<current_solver.time_for_communication1[n]<<"\n";
    cout<<"Time for Communication2:   "<<current_solver.time_for_communication2[n]<<"\n";
    cout<<"Time for Graph:            "<<current_solver.time_for_graph[n]<<"\n";
    cout<<"Total time:                "<<current_solver.time_for_total[n]<<"\n";
    cout<<"<<<<< Simulation finished. \n";
    cout<<"**********************************************************\n";
  }else{
    reduce( world, current_solver.time_for_total[world.rank()], mpi::maximum<Real>(), 0);
  }
#else

  //TODO complete the solver
  current_solver.Initialize_case(world);
  
  world.barrier();
  world.abort(-1);

  Real Time = Ini.Get_start_time();  //computation time
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Simulation starts at time: "<< Time<< "\n";
  }

  current_solver.time_total.restart();

  for(int n=1; n<=Ini.Get_num_timestep(); n++){
    Set_current_time(Time);

    if (world.rank() == 0) cout<<"Time:  "<<Time<<"\n";

    Eularian_Cat_Solver(Ini, n, world);

    Time += Ini.Get_dt();     //run time for next time step
  }

  current_solver.time_for_total[world.rank()] = current_solver.time_total.elapsed();

  if (world.rank() == 0){
    reduce( world, current_solver.time_for_total[world.rank()], current_solver.time_for_total[world.size()],  mpi::maximum<Real>(), 0);
    int n = 0;
  #ifdef _MPI_
    n = world.size();
  #else
    n = 0;
  #endif
    cout<<"**********************************************************\n";
    cout<<"Total time:                "<<current_solver.time_for_total[n]<<"\n";
    cout<<"<<<<< Simulation finished. \n";
    cout<<"**********************************************************\n";
  }else{
    reduce( world, current_solver.time_for_total[world.rank()], mpi::maximum<Real>(), 0);
  }

#endif


#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Load_solver finished\n";
  out.close();
#endif
}
//--------------------------------------------------
// Set_current_time
//--------------------------------------------------
void Solver::Set_current_time(Real Time){
  run_time = Time;
}
//--------------------------------------------------
// MRSPH_Solver
//--------------------------------------------------
#if !defined (_EULARIAN_CAT_) && !defined (_MESH_GENERATION_)
void Solver::MRSPH_Solver(Initialization &Ini, int n, communicator &world)
{
  integral_time = 0.0;
  current_solver.i_iter = n;
  while (integral_time < output_dt - 1.0e-10){

      current_solver.time_simulation.restart();
      current_solver.time_rest.restart();
    current_solver.Set_timestep(world);
      current_solver.time_for_rest[world.rank()] += current_solver.time_rest.elapsed();
      current_solver.time_for_simulation[world.rank()] += current_solver.time_simulation.elapsed();

    timestep = current_solver.Get_timestep();
    if (integral_time >= (output_dt - timestep)){
      timestep = output_dt - integral_time;
      current_solver.Reset_timestep(timestep);
    }
    integral_time += timestep;

    current_solver.Set_current_time(run_time+integral_time);

    current_solver.SPH_run_simulation (world);

  #ifdef _MPI_
    if (world.rank() == 0)  cout<<"<<<<<<integral_time: "<<integral_time<<" error "<<current_solver.error_global<<" imbalance "<<current_solver.imbalance_global<<"\n";
  #else
    if (world.rank() == 0)  cout<<"<<<<<<integral_time: "<<integral_time<<"\n";
  #endif
  #ifdef _MEM_CHECK_
    current_solver.Check_memory_consumption(world);
  #endif
  #ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<integral_time "<<integral_time<<" finished\n";
  out<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
  out.close();
  #endif
  }

    current_solver.time_simulation.restart();
    current_solver.time_rest.restart();
  #if SYM_DIM == 0
  // check if the total force is conservative
    current_solver.Check_total_force(world);
  #endif
  // Calculate total energy
  current_solver.Get_energy(world);
//   current_solver.Check_total_mass_energy(world);
    current_solver.time_for_rest[world.rank()] += current_solver.time_rest.elapsed();
    current_solver.time_for_simulation[world.rank()] += current_solver.time_simulation.elapsed();

  // Output computational result
  current_solver.Output_plt_file(n, 1, world);
//   current_solver.Output_fluid_infor(n, world);
//  current_solver.Output_every_level_dat(n, world);
  current_solver.Output_runtime_info(0, world);
  #if PERI_DIM != 0
//    current_solver.Output_PBC_plt_file(n, world);
  #endif
  #if SYM_DIM != 0
//    current_solver.Output_SBC_plt_file(n, world);
  #endif
  #ifdef _MPI_
//  current_solver.Output_pov_ray_file(n, world);
//  current_solver.Output_buffer_vtk_file(n, world);
//  current_solver.Output_buffer_plt_file(n, world);
    current_solver.Output_global_info(1, world);
    current_solver.sph_voronoi.Output_plt_file(0, n, world);
  #endif
  #ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<iteration "<<n<<" finished\n";
  out.close();
  #endif
//  current_solver.Output_restart_file();
}
#endif
//--------------------------------------------------
// particle-based mesh generation
// solver
//--------------------------------------------------
#if !defined (_EULARIAN_CAT_) && defined (_MESH_GENERATION_)
void Solver::SPH_Mesh_Generation_Solver(Initialization &Ini, int n, communicator &world)
{
  integral_time = 0.0;
  current_solver.i_iter = n;
  while (integral_time < output_dt){

    integral_time += 1.;

    current_solver.Set_current_time(run_time+integral_time);

    current_solver.SPH_run_mesh_generation (world);

    timestep = current_solver.Get_timestep();
        
  #ifdef _MPI_
    if (world.rank() == 0)  cout<<"<<<<<<iteration: "<<integral_time<<" current timestep size: "<<timestep<<" error "<<current_solver.error_global<<" imbalance "<<current_solver.imbalance_global<<"\n";
  #else
    if (world.rank() == 0)  cout<<"<<<<<<iteration: "<<integral_time<<" current timestep size: "<<timestep<<"\n";
  #endif
  #ifdef _MEM_CHECK_
    current_solver.Check_memory_consumption(world);
  #endif
  #ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<integral_time "<<integral_time<<" finished\n";
  out<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
  out.close();
  #endif
  }
    current_solver.time_simulation.restart();
    current_solver.time_rest.restart();
  #if SYM_DIM == 0
  // check if the total force is conservative
    // current_solver.Check_total_force(world);
  #else
//     current_solver.Output_SBC_plt_file(n, world);
  #endif
    current_solver.time_for_rest[world.rank()] += current_solver.time_rest.elapsed();
    current_solver.time_for_simulation[world.rank()] += current_solver.time_simulation.elapsed();

  // Output computational result
  // current_solver.Output_plt_file(n, 1, world);
  // current_solver.Output_paraview_file(n, 1, world);
  current_solver.Post_processing (n, world);
  current_solver.Output_simulation_infor(1, world);
//  current_solver.Output_every_level_dat(n, world);
//   current_solver.Output_runtime_info(0, world);
  #ifdef _MPI_
//  current_solver.Output_pov_ray_file(n, world);
//  current_solver.Output_buffer_vtk_file(n, world);
//  current_solver.Output_buffer_plt_file(n, world);
    current_solver.Output_global_info(1, world);
    current_solver.sph_voronoi.Output_vtk_file(0, n, world);
//     current_solver.sph_voronoi.Output_plt_file(0, n, world);
  #endif
  #ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<iteration "<<n<<" finished\n";
  out.close();
  #endif
}
#endif
//--------------------------------------------------
// Eularian Cartisian Solver (structured mesh)
//--------------------------------------------------
#if defined (_EULARIAN_CAT_)
void Solver::Eularian_Cat_Solver(Initialization &Ini, int n, communicator &world)
{
  integral_time = 0.0;
  while (integral_time < output_dt - 1.0e-10){

//     current_solver.Set_timestep(world);

//     timestep = current_solver.Get_timestep();
    if (integral_time >= (output_dt - timestep)){
      timestep = output_dt - integral_time;
//       current_solver.Reset_timestep(timestep);
    }
    integral_time += timestep;

//     current_solver.Set_current_time(run_time+integral_time);

//     current_solver.Eularian_cat_run_simulation (world);

    if (world.rank() == 0)  cout<<"<<<<<<integral_time: "<<integral_time<<"\n";

  #ifdef _TEST_
    char    filename[256];
    sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
    ofstream out(filename, ios::app);
    out<<"<<<<<integral_time "<<integral_time<<" finished\n";
    out<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
    out.close();
  #endif
  }

  #ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<iteration "<<n<<" finished\n";
  out.close();
  #endif
}
#endif
