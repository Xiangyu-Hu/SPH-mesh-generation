#include "sph.h"
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
#include "voronoi.h"

/*********************************************************/
/*                                                       */
/*  Functions defined in class "Voronoi"&"Vo_particle"   */
/*                                                       */
/*********************************************************/

//-------------------------------------------------------
// Important note: if you are using TBB concurrent_vector
// && Vo_particlepool, please don't use
// Vo_particlepool.free() operation. This will cause crash.
// The reason may be the incompatibility of this two Libs.
// If you want to dynamically change the number of VPs,
// please use setups identical to what we do with sph
// particle.
//-------------------------------------------------------

//-------------------------------------------------------
// Get_coord_of_VP_at_iRank
//-------------------------------------------------------
my_real Voronoi::Get_coord_of_VP_at_iRank(int iRank)
{
  if (iRank >= int(vo_particle.size())){
    cout<<"<<<<< ERROR iRank value in getting VP position"<<endl;
    my_real rr; my_set_const (rr, 0.);
    return rr;
  }else
    return vo_particle[iRank]->coord;
}
//-------------------------------------------------------
// Get_h_min_of_VP_at_iRank
//-------------------------------------------------------
Real Voronoi::Get_h_min_of_VP_at_iRank(int iRank)
{
  if (iRank >= int(vo_particle.size())){
    cout<<"<<<<< ERROR iRank value in getting VP position"<<endl;
    Real rr = 0.;
    return rr;
  }else
    return vo_particle[iRank]->h_min;
}
//-------------------------------------------------------
// Set_VP_coords
//-------------------------------------------------------
void Voronoi::Set_VP_coords_and_hmin(std::vector<my_real> &vp_coords, std::vector<Real> &vp_scale, communicator &world)
{
  if (world.rank() == 0){
    if (vp_coords.size() != world.size()){
      cout<<"<<<<< ERROR in vp_coords.size() "<<vp_coords.size()<<endl;
      world.abort(-1);
    }
    for (int i = 0; i < vo_particle.size(); i++){
      my_set_data(vo_particle[i]->coord, vp_coords[i]);
      vo_particle[i]->h_min = vp_scale[i];
    }
  }

  serialization_vector <p_Vo_particle> vp_info;
  vp_info.Vector.clear();
  vp_info.mem_size = int(vo_particle.size());
  vp_info.tag      = 0;

  if (world.rank() == 0){
    for (int i = 0; i < vo_particle.size(); i++){
      vo_particle[i]->tag = 4;
      vp_info.Vector.push_back(vo_particle[i]);
    }
    vp_info.tag      = 1;
    broadcast(world,vp_info,0);

  }else{

    broadcast(world,vp_info,0);
    for (int i = 0; i < vo_particle.size(); i++){
      if (vo_particle[i]->id != vp_info.Vector[i]->id){
        cout<<"VP id is wrong!!!\n"; world.abort(-1);
      }
      my_set_data (vo_particle[i]->coord, vp_info.Vector[i]->coord);
      vo_particle[i]->h_min = vp_info.Vector[i]->h_min;
    }
  }
  if (world.rank() == 0){
    vp_info.Vector.clear();
  }else{
    static affinity_partitioner ap;
    parallel_for( blocked_range<concurrent_vector<p_Vo_particle>::iterator>(vp_info.Vector.begin(), vp_info.Vector.end()),
        [&](const blocked_range<concurrent_vector<p_Vo_particle>::iterator>& r){
      for(concurrent_vector<p_Vo_particle>::iterator it=r.begin(); it!=r.end(); ++it){
        if (NULL != *it){
          delete (*it); // release memory of exchange_particles in heap
          *it = NULL;
        }
      }  
    }, ap);
    vp_info.Vector.clear();
    vp_info.Vector.shrink_to_fit();
    if (int(vp_info.Vector.capacity()) != 0){
      cout<<"Vector memory is not released!!!\n"; world.abort(-1);
    }
  }

  world.barrier();

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Set_VP_coords finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// random the particle distribution
//-------------------------------------------------------
Real Voronoi::Hex_close_packing_inbox(my_real box_l_, my_real box_, int dimX, int dimY, int dimZ, communicator &world)
{
  serialization_vector <p_Vo_particle> vp_info;
  vp_info.Vector.clear();
  vp_info.mem_size = int(vo_particle.size());
  vp_info.tag      = 0;

  Real RX = LONG_MAX;
  Real RY = LONG_MAX;
  Real RZ = LONG_MAX;
  Real R  = LONG_MAX;
  
  RX = DIM_X ? box_.i/(2.*dimX+DIM_Y+DIM_Z) : 0.;
  RY = DIM_Y ? box_.j/(2.*dimY)   : 0.;
  RZ = DIM_Z ? box_.k/(2.*dimZ)   : 0.;
  
  R = DIM_X ? AMIN1 (RX, R) : R;
  R = DIM_Y ? AMIN1 (RY, R) : R;
  R = DIM_Z ? AMIN1 (RZ, R) : R;
    
  if (world.rank() == 0){
    int local_id = 0;
    for (int i = 0; i < dimX; i++)
      for (int j = 0; j < dimY; j++)
	for (int k = 0; k < dimZ; k++){
	  vo_particle[local_id]->coord.i = DIM_X ? (Real(2*i)+1 + Real(j%2) + Real(k%2))*RX + box_l_.i : 0.;
	  vo_particle[local_id]->coord.j = DIM_Y ? (Real(sqrt(3.)*j)+1. + sqrt(3.)/3.*Real(k%2))*RY + box_l_.j : 0.;
	  vo_particle[local_id]->coord.k = DIM_Z ? (2*sqrt(6.)/3.*Real(k)+1.)*RZ + box_l_.k : 0.;
	  local_id ++;
	}

    for (int i = 0; i < vo_particle.size(); i++){
      vo_particle[i]->tag = 1;
      vp_info.Vector.push_back(vo_particle[i]);
    }
    vp_info.tag      = 1;
    broadcast(world,vp_info,0);

  }else{

    broadcast(world,vp_info,0);
    for (int i = 0; i < vo_particle.size(); i++){
      if (vo_particle[i]->id != vp_info.Vector[i]->id){
        cout<<"VP id is wrong!!!\n"; world.abort(-1);
      }
      my_set_data (vo_particle[i]->coord, vp_info.Vector[i]->coord);
    }
  }
  if (world.rank() == 0){
    vp_info.Vector.clear();
  }else{
    static affinity_partitioner ap;
    parallel_for( blocked_range<concurrent_vector<p_Vo_particle>::iterator>(vp_info.Vector.begin(), vp_info.Vector.end()),
        [&](const blocked_range<concurrent_vector<p_Vo_particle>::iterator>& r){
      for(concurrent_vector<p_Vo_particle>::iterator it=r.begin(); it!=r.end(); ++it){
        if (NULL != *it){
          delete (*it); // release memory of exchange_particles in heap
          *it = NULL;
        }
      }  
    }, ap);
    vp_info.Vector.clear();
    vp_info.Vector.shrink_to_fit();
    if (int(vp_info.Vector.capacity()) != 0){
      cout<<"Vector memory is not released!!!\n"; world.abort(-1);
    }
  }
  world.barrier();
  
  return R;
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Hex_close_packing_inbox finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// random the particle distribution
//-------------------------------------------------------
void Voronoi::Random_particle_distribution_inbox(my_real box_l_, my_real box_, communicator &world)
{
  serialization_vector <p_Vo_particle> vp_info;
  vp_info.Vector.clear();
  vp_info.mem_size = int(vo_particle.size());
  vp_info.tag      = 0;

  if (world.rank() == 0){
    // get the random number
    srand((unsigned)time(NULL));

    static affinity_partitioner ap;
    parallel_for( blocked_range<int>(0, int(vo_particle.size())),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        my_real position;
#if (P_DIM_X)
       position.i = rand()/ double(RAND_MAX+1.0);
       vo_particle[i]->coord.i = position.i*box_.i + rand()/ double(RAND_MAX+1.0)*1.e-3 + box_l_.i;
       vo_particle[i]->coord.i = AMAX1(box_l_.i+1.e-20,AMIN1(vo_particle[i]->coord.i, box_.i+box_l_.i-1.e-20));
#endif
#if (P_DIM_Y)
       position.j = rand()/ double(RAND_MAX+1.0);
       vo_particle[i]->coord.j = position.j*box_.j + rand()/ double(RAND_MAX+1.0)*1.e-3 + box_l_.j;
       vo_particle[i]->coord.j = AMAX1(box_l_.j+1.e-20,AMIN1(vo_particle[i]->coord.j, box_.j+box_l_.j-1.e-20));
#endif
#if (P_DIM_Z)
       position.k = rand()/ double(RAND_MAX+1.0);
       vo_particle[i]->coord.k = position.k*box_.k + rand()/ double(RAND_MAX+1.0)*1.e-3 + box_l_.k;
       vo_particle[i]->coord.k = AMAX1(box_l_.k+1.e-20,AMIN1(vo_particle[i]->coord.k, box_.k+box_l_.k-1.e-20));
#endif
      }
    }, ap);

    for (int i = 0; i < vo_particle.size(); i++){
      vo_particle[i]->tag = 1;
      vp_info.Vector.push_back(vo_particle[i]);
    }
    vp_info.tag      = 1;
    broadcast(world,vp_info,0);

  }else{

    broadcast(world,vp_info,0);
    for (int i = 0; i < vo_particle.size(); i++){
      if (vo_particle[i]->id != vp_info.Vector[i]->id){
        cout<<"VP id is wrong!!!\n"; world.abort(-1);
      }
      my_set_data (vo_particle[i]->coord, vp_info.Vector[i]->coord);
    }
  }
  if (world.rank() == 0){
    vp_info.Vector.clear();
  }else{
    static affinity_partitioner ap;
    parallel_for( blocked_range<concurrent_vector<p_Vo_particle>::iterator>(vp_info.Vector.begin(), vp_info.Vector.end()),
        [&](const blocked_range<concurrent_vector<p_Vo_particle>::iterator>& r){
      for(concurrent_vector<p_Vo_particle>::iterator it=r.begin(); it!=r.end(); ++it){
        if (NULL != *it){
          delete (*it); // release memory of exchange_particles in heap
          *it = NULL;
        }
      }  
    }, ap);
    vp_info.Vector.clear();
    vp_info.Vector.shrink_to_fit();
    if (int(vp_info.Vector.capacity()) != 0){
      cout<<"Vector memory is not released!!!\n"; world.abort(-1);
    }
  }
  world.barrier();
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Random_particle_distribution finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// random the particle distribution
//-------------------------------------------------------
void Voronoi::Random_particle_distribution(communicator &world)
{
  serialization_vector <p_Vo_particle> vp_info;
  vp_info.Vector.clear();
  vp_info.mem_size = int(vo_particle.size());
  vp_info.tag      = 0;

  my_real water_box, water_box_l, water_box_r;

//  water_box.i = 2.0;
//  water_box.j = 1.0;
//  water_box.k = 1.0;

//  water_box_l.i = 0.;
//  water_box_l.j = 0.;
//  water_box_l.k = 1.;

//  water_box_r.i = 2.;
//  water_box_r.j = 1.;
//  water_box_r.k = 1.;

  if (world.rank() == 0){
    // get the random number
    srand((unsigned)time(NULL));

    static affinity_partitioner ap;
    parallel_for( blocked_range<int>(0, int(vo_particle.size())),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        my_real position;
//#if (P_DIM_X)
//        position.i = rand()/ double(RAND_MAX+1.0);
//        vo_particle[i]->coord.i = position.i*water_box.i + rand()/ double(RAND_MAX+1.0)*1.e-3 + water_box_l.i;
//        vo_particle[i]->coord.i = AMAX1(water_box_l.i+1.e-20,AMIN1(vo_particle[i]->coord.i, water_box_r.i-1.e-20));
//#endif
//#if (P_DIM_Y)
//        position.j = rand()/ double(RAND_MAX+1.0);
//        vo_particle[i]->coord.j = position.j*water_box.j + rand()/ double(RAND_MAX+1.0)*1.e-3 + water_box_l.j;
//        vo_particle[i]->coord.j = AMAX1(water_box_l.j+1.e-20,AMIN1(vo_particle[i]->coord.j, water_box_r.j-1.e-20));
//#endif
//#if (P_DIM_Z)
//        position.k = rand()/ double(RAND_MAX+1.0);
//        vo_particle[i]->coord.k = position.k*water_box.k + rand()/ double(RAND_MAX+1.0)*1.e-3 + water_box_l.k;
//         vo_particle[i]->coord.k = AMAX1(water_box_l.k+1.e-20,AMIN1(vo_particle[i]->coord.k, water_box_r.k-1.e-20));
//#endif

#if (P_DIM_X)
        position.i = rand()/ double(RAND_MAX+1.0);
        vo_particle[i]->coord.i = position.i*particle_box.i + rand()/ double(RAND_MAX+1.0)*1.e-3 + particle_box_l.i;
        vo_particle[i]->coord.i = AMAX1(particle_box_l.i+1.e-20,AMIN1(vo_particle[i]->coord.i, particle_box_r.i-1.e-20));
#endif
#if (P_DIM_Y)
        position.j = rand()/ double(RAND_MAX+1.0);
        vo_particle[i]->coord.j = position.j*particle_box.j + rand()/ double(RAND_MAX+1.0)*1.e-3 + particle_box_l.j;
        vo_particle[i]->coord.j = AMAX1(particle_box_l.j+1.e-20,AMIN1(vo_particle[i]->coord.j, particle_box_r.j-1.e-20));
#endif
#if (P_DIM_Z)
        position.k = rand()/ double(RAND_MAX+1.0);
        vo_particle[i]->coord.k = position.k*particle_box.k + rand()/ double(RAND_MAX+1.0)*1.e-3 + particle_box_l.k;
         vo_particle[i]->coord.k = AMAX1(particle_box_l.k+1.e-20,AMIN1(vo_particle[i]->coord.k, particle_box_r.k-1.e-20));
#endif
      }
    }, ap);

    for (int i = 0; i < vo_particle.size(); i++){
      vo_particle[i]->tag = 1;
      vp_info.Vector.push_back(vo_particle[i]);
    }
    vp_info.tag      = 1;
    broadcast(world,vp_info,0);

  }else{

    broadcast(world,vp_info,0);
    for (int i = 0; i < vo_particle.size(); i++){
      if (vo_particle[i]->id != vp_info.Vector[i]->id){
        cout<<"VP id is wrong!!!\n"; world.abort(-1);
      }
      my_set_data (vo_particle[i]->coord, vp_info.Vector[i]->coord);
    }
  }
  if (world.rank() == 0){
    vp_info.Vector.clear();
  }else{
    static affinity_partitioner ap;
    parallel_for( blocked_range<concurrent_vector<p_Vo_particle>::iterator>(vp_info.Vector.begin(), vp_info.Vector.end()),
        [&](const blocked_range<concurrent_vector<p_Vo_particle>::iterator>& r){
      for(concurrent_vector<p_Vo_particle>::iterator it=r.begin(); it!=r.end(); ++it){
        if (NULL != *it){
          delete (*it); // release memory of exchange_particles in heap
          *it = NULL;
        }
      }  
    }, ap);
    vp_info.Vector.clear();
    vp_info.Vector.shrink_to_fit();
    if (int(vp_info.Vector.capacity()) != 0){
      cout<<"Vector memory is not released!!!\n"; world.abort(-1);
    }
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Random_particle_distribution finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// Initialization for Voronoi
//-------------------------------------------------------
void Voronoi::Initialize(SPH *sph, int flag, communicator &world){

  sph_particle        = sph->particle;
  num_sph_particle    = sph->total_num_particle;
  num_partition       = sph->glbl_total_num_color;
  convergenced        = 0;
  Max_iterate_number  = 4000;
  error_tolerance     = 0.05;
  scale_coeff         = 0.;
  my_set_data (domain, sph->domain);
  my_set_data (box_l, sph->box_l);
  my_set_data (box_r, sph->box_r);
  my_set_data (particle_box, sph->particle_box);
  my_set_data (particle_box_l, sph->particle_box_l);
  my_set_data (particle_box_r, sph->particle_box_r);
  enlarge_ratio         = pow(2.,20.)/domain.i;
  relax_ratio           = 0.80;
  total_energy          = 0.;
  total_interface_area  = 0.;
  error_max             = 0.;
  error_record          = 0.;

#if PERI_DIM != 0
  cross_boundary.i = 0;
  cross_boundary.j = 0;
  cross_boundary.k = 0;
#endif

  if (flag == 1){
    // generate the partitioning particles
      for(int i=0; i!=num_partition; ++i){
        p_Vo_particle particle_new = Vo_particlepool.malloc();
        particle_new->Initialize(i);
        vo_particle.push_back(particle_new);
      }

    // initialize graph
    vo_graph.Initialize(world);
  }

  // Set up constants for the container geometry
  x_min=-0.5; x_max=0.5;
  y_min=-0.5; y_max=0.5;
  z_min=-0.5; z_max=0.5;
#if P_DIM_X == 1
  x_min=box_l.i, x_max=box_r.i;
#endif
#if P_DIM_Y == 1
  y_min=box_l.j, y_max=box_r.j;
#endif
#if P_DIM_Z == 1
  z_min=box_l.k, z_max=box_r.k;
#endif
  P_X = P_Y = P_Z = false;
#if PERI_DIM_X == 1 && P_DIM_X == 1
  P_X = true;
#endif
#if PERI_DIM_Y == 1 && P_DIM_Y == 1
  P_Y = true;
#endif
#if PERI_DIM_Z == 1 && P_DIM_Z == 1
  P_Z = true;
#endif
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Class voronoi initialized\n";
  out.close();
#endif
}
//-------------------------------------------------------
// get the cpu index of certain position
//-------------------------------------------------------
int Voronoi::Get_CPU_id(my_real p_coord)
{
  Real dist_min = 1.e20;
  int index_record;

  for(int i=0; i<int(vo_particle.size()); i++){
    Real   dist  = get_distance_2p (p_coord, vo_particle[i]->coord);
    if (dist < dist_min){
      dist_min = dist;
      index_record = i;
    }
  }
  return index_record;
}
//-------------------------------------------------------
// Set the target mass for each VP
//-------------------------------------------------------
void Voronoi::Set_target_mass(communicator &world)
{
  static affinity_partitioner ap;
  total_mass = 0.;
  for(int i=0; i!=num_sph_particle; ++i){
    total_mass += sph_particle[i]->p_mass;
  }

  all_reduce(world, total_mass, glbl_total_mass, std::plus<Real>());

  mass_target = glbl_total_mass/num_partition;
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Set_target_mass finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// calcuate VP mass, mass_center, energy and pressure
//-------------------------------------------------------
void Voronoi::Get_VP_info_and_VD_generation(communicator &world)
{
  //   char    filename2[256];
  // sprintf(filename2,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  // ofstream out2(filename2, ios::app);
  
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, num_partition),
           [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      vo_particle[i]->Reset();
    }
  }, ap);

  // Set up the number of blocks that the container is divided into
  int n_x, n_y, n_z;
  double x,y,z;

  particle_order po;

  // Create a pre-container class to import the input file and guess the
  // best computational grid size to use.
  voro::pre_container pcon(x_min,x_max,y_min,y_max,z_min,z_max,P_X,P_Y,P_Z);

  wall_plane plane1(0., 0., 1., fabs(z_max), num_partition + 1);
  wall_plane plane2(0., 1., 0., fabs(y_max), num_partition + 2);
  wall_plane plane3(1., 0., 0., fabs(x_max), num_partition + 3);
  wall_plane plane4(0., 0., -1., fabs(z_min), num_partition + 4);
  wall_plane plane5(0., -1., 0., fabs(y_min), num_partition + 5);
  wall_plane plane6(-1., 0., 0., fabs(x_min), num_partition + 6);

  for(int i=0;i<num_partition;i++){
    x = y = z = 0.;
#if P_DIM_X == 1
    x = vo_particle[i]->coord.i;
#endif
#if P_DIM_Y == 1
    y = vo_particle[i]->coord.j;
#endif
#if P_DIM_Z == 1
    z = vo_particle[i]->coord.k;
#endif
    pcon.put(i,x,y,z);
  }

  pcon.guess_optimal(n_x,n_y,n_z);

  // Set up the container class and import the particles from the
  // pre-container
  voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,P_X,P_Y,P_Z,8);

  // Add wall to the boundary
#if P_DIM_Z == 1 && PERI_DIM_Z == 0
  con.add_wall(plane1);
  con.add_wall(plane4);
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 0
  con.add_wall(plane2);
  con.add_wall(plane5);
#endif
#if P_DIM_X == 1 && PERI_DIM_X == 0
  con.add_wall(plane3);
  con.add_wall(plane6);
#endif

  // add particles into the container
  for(int i=0;i<num_partition;i++){
    x = y = z = 0.;
#if P_DIM_X == 1
    x = vo_particle[i]->coord.i;
#endif
#if P_DIM_Y == 1
    y = vo_particle[i]->coord.j;
#endif
#if P_DIM_Z == 1
    z = vo_particle[i]->coord.k;
#endif
    con.put(po,i,x,y,z);
  }
  
  parallel_for( blocked_range<int>(0, num_sph_particle),
           [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      p_Particle cur_sph_particle = sph_particle[i];
      int index = 0;
      Real rx, ry, rz, xx, yy, zz;
      xx = P_DIM_X==1 ? cur_sph_particle->coord.i : 0;
      yy = P_DIM_Y==1 ? cur_sph_particle->coord.j : 0;
      zz = P_DIM_Z==1 ? cur_sph_particle->coord.k : 0;
      con.find_voronoi_cell(xx,yy,zz,rx,ry,rz,index);
      cur_sph_particle->color = index;
    }
  }, ap);
  
  parallel_for( blocked_range<int>(0, num_partition),
           [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      p_Vo_particle cur_vp = vo_particle[i];
      for (int ip = 0; ip < num_sph_particle; ip++){
        p_Particle cur_sp = sph_particle[ip];
        if (cur_sp->color == cur_vp->color){
          cur_vp->mass        += cur_sp->p_mass;
          my_real        temp  = my_multiply_const(cur_sp->coord, cur_sp->p_mass);
          cur_vp->mass_center  = my_add_data(cur_vp->mass_center, temp);
          Real           dist  = get_distance_2p (cur_vp->coord, cur_sp->coord);
          cur_vp->energy      += cur_sp->p_mass * dist;
        }
      }
    }
  }, ap);
  
  serialization_vector <p_Vo_particle> exchange_vector;
  exchange_vector.Vector.clear();

  // out2<<"OK1 in rank:"<<world.rank()<<endl;

  if (world.rank() == 0){
    serialization_vector <p_Vo_particle> *gather_exchange_vector;
    gather_exchange_vector = new serialization_vector <p_Vo_particle>[world.size()];

    parallel_for( blocked_range<int>(0, world.size()),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        gather_exchange_vector[i].Vector.clear();
        gather_exchange_vector[i].mem_size = 0;
        gather_exchange_vector[i].tag = 0;
      }
    }, ap);

    // out2<<"OK2 in rank:"<<world.rank()<<endl;

    gather(world,exchange_vector,gather_exchange_vector,0);

    // out2<<"OK3 in rank:"<<world.rank()<<endl;

    parallel_for( blocked_range<int>(0, num_partition),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        for (int j = 1; j < world.size(); j++){
          if (vo_particle[i]->id != gather_exchange_vector[j].Vector[i]->id){
            cout<<"VP id is wrong!!!\n"; world.abort(-1);
          }
          vo_particle[i]->mass        +=  gather_exchange_vector[j].Vector[i]->mass;
          vo_particle[i]->mass_center  =  my_add_data(vo_particle[i]->mass_center,
                                          gather_exchange_vector[j].Vector[i]->mass_center);
          vo_particle[i]->energy      +=  gather_exchange_vector[j].Vector[i]->energy;
        }
      }
    }, ap);

    // out2<<"OK4 in rank:"<<world.rank()<<endl;

    parallel_for( blocked_range<int>(0, num_partition),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        p_Vo_particle cur_vp = vo_particle[i];
        cur_vp->mass_center  = my_multiply_const(cur_vp->mass_center, 1./cur_vp->mass);
#if P_DIM_X == 0
        cur_vp->mass_center.i = 0.;
#endif
#if P_DIM_Y == 0
        cur_vp->mass_center.j = 0.;
#endif
#if P_DIM_Z == 0
        cur_vp->mass_center.k = 0.;
#endif
        cur_vp->P            = powern (mass_target/cur_vp->mass, 1);
        cur_vp->error        = fabs(cur_vp->mass - mass_target)/mass_target;
      }
    }, ap);

    // out2<<"OK5 in rank:"<<world.rank()<<endl;

    error_max = 0.;
    total_energy   = 0.;
    for(int i=0; i!=num_partition; ++i){
      p_Vo_particle cur_vp = vo_particle[i];
      error_max      = AMAX1(error_max, cur_vp->error);
      total_energy  += cur_vp->energy;
    }

    for( int i=1; i<world.size(); i++){
      parallel_for( blocked_range<concurrent_vector<p_Vo_particle>::iterator>
                  (gather_exchange_vector[i].Vector.begin(), gather_exchange_vector[i].Vector.end()),
               [&](const blocked_range<concurrent_vector<p_Vo_particle>::iterator>& r){
        for(concurrent_vector<p_Vo_particle>::iterator it=r.begin(); it!=r.end(); ++it){
          if (NULL != *it){
            delete (*it); // release memory of exchange_particles in heap
            *it = NULL;
          }
        }
      }, ap);
      gather_exchange_vector[i].Vector.clear();
      gather_exchange_vector[i].Vector.shrink_to_fit();
      if (int(gather_exchange_vector[i].Vector.capacity()) != 0){
        cout<<"Vector memory is not released!!!\n"; world.abort(-1);
      }
    }
    delete [] gather_exchange_vector;

    // out2<<"OK6 in rank:"<<world.rank()<<endl;

    for( int i=0; i<num_partition; i++){
      p_Vo_particle cur_vp = vo_particle[i];
      my_set_const (cur_vp->f, 0.0);
    }

    int id;
    voronoicell_neighbor c;
    vector<int> neigh;
    vector<int> f_vert;
    vector<double> v;
    vector<double> face;
    vector<double> normal;
    c_loop_all cl(con);

    // out2<<"OK7 in rank:"<<world.rank()<<endl;

    if(cl.start()) do if(con.compute_cell(c,cl)) {
      // Get the position of the current particle under consideration
      cl.pos(x,y,z);
      id=cl.pid();

      // Gather information about the computed Voronoi cell
      c.neighbors(neigh);
      c.face_areas(face);
      c.normals(normal);
      c.face_vertices(f_vert);
      c.vertices(x,y,z,v);

      int     i, j;
      int     num_neighbor = 0;
      Real    scale_record = 0.;
      Real    current_P = vo_particle[id]->P;
      Real    neighbor_P = 0.;
      my_real current_coord;
      my_real neighbor_coord;
      my_set_data (current_coord, vo_particle[id]->coord);

      // iteration along neighbors
      for(i=0, j=0;i<neigh.size();i++) {

        // Get the neighbor vp coord, considering BC
        if (neigh[i] >= 0) { // if the boundary is open, neigh[i]<0
          num_neighbor ++;
          if (neigh[i] < num_partition){
            neighbor_P = vo_particle[neigh[i]]->P;
          }else if (neigh[i] > num_partition){ // neighbor is outside of the domain
            neighbor_P = current_P;
          }

          int l;
          l=3*f_vert[j+1];
          neighbor_coord.i = v[l];
          neighbor_coord.j = v[l+1];
          neighbor_coord.k = v[l+2];

          my_real       dx = my_minus_data (neighbor_coord, current_coord);
          my_real dr;
          dr.i = normal[3*i+0];
          dr.j = normal[3*i+1];
          dr.k = normal[3*i+2];
          Real  dist = fabs(get_dot_product (dx, dr)) * 2.;
          scale_record += dist;
          Real pressure_middle = 0.5*(current_P + neighbor_P);

          vo_particle[id]->f = my_add_data (vo_particle[id]->f, my_multiply_const(dr, -1*pressure_middle*face[i]));
        }
        j+=f_vert[j]+1;
      }

      vo_particle[id]->h = scale_coeff*scale_record/Real(num_neighbor);
    } while (cl.inc());

    // out2<<"OK8 in rank:"<<world.rank()<<endl;

    total_interface_area = c.surface_area();

  }else{
    exchange_vector.mem_size = 0;
    exchange_vector.tag = 0;
  
    for( int i=0; i < num_partition; i++){
      vo_particle[i]->tag = 2;
      exchange_vector.Vector.push_back(vo_particle[i]);

      // out2<<"======================="<<endl;
      // out2<<vo_particle[i]->id<<endl;
      // out2<<vo_particle[i]->mass<<endl;
      // out2<<vo_particle[i]->mass_center.i<<" "<<vo_particle[i]->mass_center.j<<" "<<vo_particle[i]->mass_center.k<<endl;
      // out2<<vo_particle[i]->energy<<endl;
    }
    exchange_vector.mem_size = int(exchange_vector.Vector.size());
    exchange_vector.tag = 1;
    // out2<<"OK2 in rank:"<<world.rank()<<endl;
    gather(world,exchange_vector,0);
    // out2<<"OK3 in rank:"<<world.rank()<<endl;
  }
  con.clear();
  exchange_vector.Vector.clear();
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Get_VP_info_and_VD_generation finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// update the VP positon
//-------------------------------------------------------
void Voronoi::Update_position(communicator &world)
{
  if (world.rank() == 0){
    // update the acceleration and timestep
    static affinity_partitioner ap;
    parallel_for( blocked_range<int>(0, int(vo_particle.size())),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        p_Vo_particle cur_vp = vo_particle[i];
        cur_vp->a = my_multiply_const (cur_vp->f, 1./cur_vp->mass);
        Real accel = get_distance (cur_vp->a);
        cur_vp->timestep = 0.25*sqrt(cur_vp->h/(accel + 1.e-20));
      }
    }, ap);

    glbl_timestep = 1.e20;
    // get the timestep
    for(int i=0; i!=int(vo_particle.size()); ++i){
      glbl_timestep = AMIN1(vo_particle[i]->timestep,glbl_timestep);
    }

    // update VP position
    parallel_for( blocked_range<int>(0, int(vo_particle.size())),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        p_Vo_particle cur_vp = vo_particle[i];
        my_real dx; my_set_const (dx, 0.0);
        dx = my_multiply_const (cur_vp->a, glbl_timestep*glbl_timestep*0.5*relax_ratio);
        cur_vp->coord   = my_add_data (cur_vp->coord, dx);
#if P_DIM_X == 1 && PERI_DIM_X == 0
        cur_vp->coord.i = AMAX1(box_l.i+1.e-20, AMIN1(cur_vp->coord.i, box_r.i-1.e-20));
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 0
        cur_vp->coord.j = AMAX1(box_l.j+1.e-20, AMIN1(cur_vp->coord.j, box_r.j-1.e-20));
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 0
        cur_vp->coord.k = AMAX1(box_l.k+1.e-20, AMIN1(cur_vp->coord.k, box_r.k-1.e-20));
#endif

#if P_DIM_X == 1 && PERI_DIM_X == 1
        if ( cur_vp->coord.i > box_r.i + 1.e-10) cur_vp->coord.i -= domain.i;
        else if ( cur_vp->coord.i < box_l.i - 1.e-10) cur_vp->coord.i += domain.i;
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 1
        if ( cur_vp->coord.j > box_r.j + 1.e-10) cur_vp->coord.j -= domain.j;
        else if ( cur_vp->coord.j < box_l.j - 1.e-10) cur_vp->coord.j += domain.j;
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 1
        if ( cur_vp->coord.k > box_r.k + 1.e-10) cur_vp->coord.k -= domain.k;
        else if ( cur_vp->coord.k < box_l.k - 1.e-10) cur_vp->coord.k += domain.k;
#endif

        cur_vp->CVT_shift     = my_minus_data (cur_vp->mass_center, cur_vp->coord);
        Real shift            =  get_distance (cur_vp->CVT_shift);
        Real shift_constr     = AMIN1(shift, 1./32.*cur_vp->h);
        cur_vp->timestep_CVT  = shift_constr/(shift + 1.e-20);
      }
    }, ap);
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<update the VP positon finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// update the particle position according to 
// the centroidal Voronoi tessellation
//-------------------------------------------------------
void Voronoi::Update_position_CVT(communicator &world)
{
  serialization_vector <p_Vo_particle> vp_info;
  vp_info.Vector.clear();
  vp_info.mem_size = int(vo_particle.size());
  vp_info.tag      = 0;
  if (world.rank() == 0){

    glbl_timestep_CVT = 1.e20;
    // get the timestep_CVT
    for(int i=0; i!=int(vo_particle.size()); ++i){
      glbl_timestep_CVT = AMIN1(vo_particle[i]->timestep_CVT, glbl_timestep_CVT);
    }

    static affinity_partitioner ap;
    parallel_for( blocked_range<int>(0, int(vo_particle.size())),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        p_Vo_particle cur_vp = vo_particle[i];
        cur_vp->CVT_shift    = my_multiply_const (cur_vp->CVT_shift, glbl_timestep_CVT*(1-relax_ratio));
        cur_vp->coord        = my_add_data (cur_vp->coord, cur_vp->CVT_shift);
#if P_DIM_X == 1 && PERI_DIM_X == 0
        cur_vp->coord.i = AMAX1(box_l.i+1.e-20, AMIN1(cur_vp->coord.i, box_r.i-1.e-20));
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 0
        cur_vp->coord.j = AMAX1(box_l.j+1.e-20, AMIN1(cur_vp->coord.j, box_r.j-1.e-20));
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 0
        cur_vp->coord.k = AMAX1(box_l.k+1.e-20, AMIN1(cur_vp->coord.k, box_r.k-1.e-20));
#endif

#if P_DIM_X == 1 && PERI_DIM_X == 1
        if ( cur_vp->coord.i > box_r.i + 1.e-10) cur_vp->coord.i -= domain.i;
        else if ( cur_vp->coord.i < box_l.i - 1.e-10) cur_vp->coord.i += domain.i;
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 1
        if ( cur_vp->coord.j > box_r.j + 1.e-10) cur_vp->coord.j -= domain.j;
        else if ( cur_vp->coord.j < box_l.j - 1.e-10) cur_vp->coord.j += domain.j;
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 1
        if ( cur_vp->coord.k > box_r.k + 1.e-10) cur_vp->coord.k -= domain.k;
        else if ( cur_vp->coord.k < box_l.k - 1.e-10) cur_vp->coord.k += domain.k;
#endif
      }
    }, ap);

    for (int i = 0; i < vo_particle.size(); i++){
      vo_particle[i]->tag = 1;
      vp_info.Vector.push_back(vo_particle[i]);
    }
    vp_info.tag      = 1;
    broadcast(world,vp_info,0);

  }else{
  
    broadcast(world,vp_info,0);
    for (int i = 0; i < vo_particle.size(); i++){
      if (vo_particle[i]->id != vp_info.Vector[i]->id){
        cout<<"VP id is wrong!!!\n"; world.abort(-1);
      }
      my_set_data (vo_particle[i]->coord, vp_info.Vector[i]->coord);
    }
  }

  if (world.rank() == 0){
    vp_info.Vector.clear();
  }else{
    static affinity_partitioner ap;
    parallel_for( blocked_range<concurrent_vector<p_Vo_particle>::iterator>(vp_info.Vector.begin(), vp_info.Vector.end()),
        [&](const blocked_range<concurrent_vector<p_Vo_particle>::iterator>& r){
      for(concurrent_vector<p_Vo_particle>::iterator it=r.begin(); it!=r.end(); ++it){
        if (NULL != *it){
          delete (*it); // release memory of exchange_particles in heap
          *it = NULL;
        }
      }  
    }, ap);
    vp_info.Vector.clear();
    vp_info.Vector.shrink_to_fit();
    if (int(vp_info.Vector.capacity()) != 0){
      cout<<"Vector memory is not released!!!\n"; world.abort(-1);
    }
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_position_CVT finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// update voronoi particle position according to 
// fluid field information 
//-------------------------------------------------------
void Voronoi::Update_vp_position_mean_velo(SPH *sph, communicator &world)
{
  // update vp coord in each processor
  p_Vo_particle cur_vp = vo_particle[world.rank()];
  cur_vp->coord = my_add_data (cur_vp->coord, my_multiply_const(sph->v_avg, sph->glbl_timestep));

#if P_DIM_X == 1 && PERI_DIM_X == 0
  cur_vp->coord.i = AMAX1(box_l.i+1.e-20, AMIN1(cur_vp->coord.i, box_r.i-1.e-20));
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 0
  cur_vp->coord.j = AMAX1(box_l.j+1.e-20, AMIN1(cur_vp->coord.j, box_r.j-1.e-20));
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 0
  cur_vp->coord.k = AMAX1(box_l.k+1.e-20, AMIN1(cur_vp->coord.k, box_r.k-1.e-20));
#endif

#if P_DIM_X == 1 && PERI_DIM_X == 1
  if ( cur_vp->coord.i > box_r.i + 1.e-10) cur_vp->coord.i -= domain.i;
  else if ( cur_vp->coord.i < box_l.i - 1.e-10) cur_vp->coord.i += domain.i;
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 1
  if ( cur_vp->coord.j > box_r.j + 1.e-10) cur_vp->coord.j -= domain.j;
  else if ( cur_vp->coord.j < box_l.j - 1.e-10) cur_vp->coord.j += domain.j;
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 1
  if ( cur_vp->coord.k > box_r.k + 1.e-10) cur_vp->coord.k -= domain.k;
  else if ( cur_vp->coord.k < box_l.k - 1.e-10) cur_vp->coord.k += domain.k;
#endif

  Sync_VP_positions(sph, world);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_vp_position finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// update voronoi particle position according to 
// fluid field information 
//-------------------------------------------------------
void Voronoi::Update_vp_position_mass_center(SPH *sph, communicator &world)
{
  // update vp coord in each processor
  p_Vo_particle cur_vp = vo_particle[world.rank()];

  my_set_data (cur_vp->coord, sph->mass_center);
  cur_vp->mass = sph->total_mass;

#if P_DIM_X == 1 && PERI_DIM_X == 0
  cur_vp->coord.i = AMAX1(box_l.i+1.e-20, AMIN1(cur_vp->coord.i, box_r.i-1.e-20));
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 0
  cur_vp->coord.j = AMAX1(box_l.j+1.e-20, AMIN1(cur_vp->coord.j, box_r.j-1.e-20));
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 0
  cur_vp->coord.k = AMAX1(box_l.k+1.e-20, AMIN1(cur_vp->coord.k, box_r.k-1.e-20));
#endif

#if P_DIM_X == 1 && PERI_DIM_X == 1
  if ( cur_vp->coord.i > box_r.i + 1.e-10) cur_vp->coord.i -= domain.i;
  else if ( cur_vp->coord.i < box_l.i - 1.e-10) cur_vp->coord.i += domain.i;
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 1
  if ( cur_vp->coord.j > box_r.j + 1.e-10) cur_vp->coord.j -= domain.j;
  else if ( cur_vp->coord.j < box_l.j - 1.e-10) cur_vp->coord.j += domain.j;
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 1
  if ( cur_vp->coord.k > box_r.k + 1.e-10) cur_vp->coord.k -= domain.k;
  else if ( cur_vp->coord.k < box_l.k - 1.e-10) cur_vp->coord.k += domain.k;
#endif

  Sync_VP_positions(sph, world);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_vp_position finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// check the partition error
//-------------------------------------------------------
void Voronoi::Shift_coordinate_constrained
(my_int dim, SPH *sph, communicator &world)
{
  p_Vo_particle cur_vp = vo_particle[world.rank()];
#if P_DIM_X == 1 && PERI_DIM_X == 1
  if (dim.i != 0){
    cur_vp->coord.i += dim.i*domain.i;
  } 
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 1
  if (dim.j != 0){
    cur_vp->coord.j += dim.j*domain.j;
  } 
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 1
  if (dim.k != 0){
    cur_vp->coord.k += dim.k*domain.k;
  } 
#endif

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Shift_coordinate_constrained finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// syncronize particle position in all ranks
//-------------------------------------------------------
void Voronoi::Sync_VP_positions(SPH *sph, communicator &world)
{
  serialization_vector <p_Vo_particle> exchange_vector;
  exchange_vector.Vector.clear();
  if (world.rank() == 0){
    // receive all vp information in master node
    serialization_vector <p_Vo_particle> *gather_exchange_vector;
    gather_exchange_vector = new serialization_vector <p_Vo_particle>[world.size()];

    static affinity_partitioner ap;
    parallel_for( blocked_range<int>(0, world.size()),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        gather_exchange_vector[i].Vector.clear();
        gather_exchange_vector[i].mem_size = 0;
        gather_exchange_vector[i].tag = 0;
      }
    }, ap);

    gather(world,exchange_vector,gather_exchange_vector,0);

    for (int i = 1; i < vo_particle.size(); i++){
      my_set_data( vo_particle[i]->coord, gather_exchange_vector[i].Vector[0]->coord);
      vo_particle[i]->mass = gather_exchange_vector[i].Vector[0]->mass;
    }

    for( int i=1; i<world.size(); i++){
      parallel_for( blocked_range<concurrent_vector<p_Vo_particle>::iterator>
                  (gather_exchange_vector[i].Vector.begin(), gather_exchange_vector[i].Vector.end()),
               [&](const blocked_range<concurrent_vector<p_Vo_particle>::iterator>& r){
        for(concurrent_vector<p_Vo_particle>::iterator it=r.begin(); it!=r.end(); ++it){
          if (NULL != *it){
            delete (*it); // release memory of exchange_particles in heap
            *it = NULL;
          }
        }
      }, ap);
      gather_exchange_vector[i].Vector.clear();
      gather_exchange_vector[i].Vector.shrink_to_fit();
      if (int(gather_exchange_vector[i].Vector.capacity()) != 0){
        cout<<"Vector memory is not released!!!\n"; world.abort(-1);
      }
    }

    delete [] gather_exchange_vector;

  }else{
    // send vp information in cuurrent node to master
    exchange_vector.mem_size = 0;
    exchange_vector.tag = 0;

    vo_particle[world.rank()]->tag = 3;
    exchange_vector.Vector.push_back(vo_particle[world.rank()]);
    exchange_vector.mem_size = int(exchange_vector.Vector.size());
    exchange_vector.tag = 1;
    gather(world,exchange_vector,0);
  }
  exchange_vector.Vector.clear();

  serialization_vector <p_Vo_particle> vp_info;
  vp_info.Vector.clear();
  vp_info.mem_size = int(vo_particle.size());
  vp_info.tag      = 0;

  // broadcaset all vp infor to every processor
  if (world.rank() == 0){
    for (int i = 0; i < vo_particle.size(); i++){
      vo_particle[i]->tag = 3;
      vp_info.Vector.push_back(vo_particle[i]);
    }
    vp_info.tag      = 1;
    broadcast(world,vp_info,0);
  }else{
    broadcast(world,vp_info,0);
    for (int i = 0; i < vo_particle.size(); i++){
      if (vo_particle[i]->id != vp_info.Vector[i]->id){
        cout<<"VP id is wrong!!!\n"; world.abort(-1);
      }
      my_set_data (vo_particle[i]->coord, vp_info.Vector[i]->coord);
      vo_particle[i]->mass = vp_info.Vector[i]->mass;
    }
  }
  vp_info.Vector.clear();

  if (world.rank() == 0){
    vp_info.Vector.clear();
  }else{
    static affinity_partitioner ap;
    parallel_for( blocked_range<concurrent_vector<p_Vo_particle>::iterator>(vp_info.Vector.begin(), vp_info.Vector.end()),
        [&](const blocked_range<concurrent_vector<p_Vo_particle>::iterator>& r){
      for(concurrent_vector<p_Vo_particle>::iterator it=r.begin(); it!=r.end(); ++it){
        if (NULL != *it){
          delete (*it); // release memory of exchange_particles in heap
          *it = NULL;
        }
      }  
    }, ap);
    vp_info.Vector.clear();
    vp_info.Vector.shrink_to_fit();
    if (int(vp_info.Vector.capacity()) != 0){
      cout<<"Vector memory is not released!!!\n"; world.abort(-1);
    }
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Sync_VP_positions finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// check the partition error
//-------------------------------------------------------
void Voronoi::Check_partition_error(int &iter, communicator &world)
{
  if (world.rank() == 0){
    if (iter%10 == 0){
      cout<<"Iter: "<<iter<<" max_error: "<<error_max
          <<" Total energy: "<<total_energy
          <<" Interface area: "<<total_interface_area<<"\n";
    }
    if (iter%10 == 0){
      ofstream save("outdata/Error_history.dat", ios::app);
      save<<"Iter: "<<iter<<" max_error: "<<error_max
          <<" Total energy: "<<total_energy
          <<" Interface area: "<<total_interface_area<<endl;
      save.close();
    }
    error_record = error_record + error_max;
    if (iter % 40 == 0){
      error_record /= 40.;
      if (error_max <= error_tolerance && error_record <= error_tolerance)
        convergenced = 1;
      error_record = 0.;
    }
    broadcast(world, convergenced, 0);
  }else{
    broadcast(world, convergenced, 0);
  }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Check_partition_error finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// add edges
//-------------------------------------------------------
void Voronoi::Add_edges(communicator &world)
{
  unsigned int *color_list;
  color_list = new unsigned int [num_partition];
  for (int i = 0; i < num_partition; i++){
    if (i != world.rank())
      color_list[i] = 0;
    else
      color_list[i] = 1;
  }

  static affinity_partitioner ap;
  tbb::mutex  edge_Mutex;
  parallel_for( blocked_range<int>(0, num_sph_particle),
           [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      int cur_color = sph_particle[i]->color;
      if (cur_color != world.rank() && color_list[cur_color] != 1){
        tbb::mutex::scoped_lock lock(edge_Mutex);
        color_list[cur_color] = 1;
      }
    }
  }, ap);

  serialization_vector <unsigned int> exchange_vector;
  exchange_vector.Vector.clear();

  if (world.rank() == 0){

    vo_graph.Remove_all_edges(world);

    serialization_vector <unsigned int> *gather_exchange_vector;
    gather_exchange_vector = new serialization_vector <unsigned int>[world.size()];

    parallel_for( blocked_range<int>(0, world.size()),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        gather_exchange_vector[i].Vector.clear();
        gather_exchange_vector[i].mem_size = 0;
        gather_exchange_vector[i].tag = 0;
      }
    }, ap);

    gather(world,exchange_vector,gather_exchange_vector,0);

    parallel_for( blocked_range<int>(0, world.size()),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        for (int j = 0; j < num_partition; j++){
          int color_i = i;
          int color_j = -1;
          int flag    = 0;
          if (i == 0){
            flag = color_list[j];
          }else{
            flag = gather_exchange_vector[i].Vector[j];
          }
          if (flag == 1){
            color_j = j;
            if (color_i != color_j)
              vo_graph.Add_edge(color_i, color_j);
          }
        }
      }             
    }, ap);
    vo_graph.total_edges = int(num_edges(vo_graph.graph));

    delete [] gather_exchange_vector;

  }else{
    exchange_vector.mem_size = 0;
    exchange_vector.tag = 0;

    for( int i=0; i < num_partition; i++){
      exchange_vector.Vector.push_back(color_list[i]);
    }
    exchange_vector.mem_size = int(exchange_vector.Vector.size());
    exchange_vector.tag = 1;
    gather(world,exchange_vector,0);
  }
  exchange_vector.Vector.clear();
  delete [] color_list;
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Add_edges finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// construct graph to migrate sph particles
//-------------------------------------------------------
void Voronoi::Construct_graph(communicator &world)
{
  // add edges
  Add_edges(world);

  // edge coloring
  vo_graph.Handle_graph_edge_coloring(world);
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Construct_graph finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// CVP Partitioning
//-------------------------------------------------------
void Voronoi::Partitioning(SPH *sph, communicator &world)
{
  sph->time_partition.restart();

  Initialize(sph, 2, world);

  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< CVP is processing\n";
  }

  Set_target_mass(world);

  for(int iter=0; iter<=Max_iterate_number && convergenced == 0; iter++){

//    if (world.rank() == 0)
//      cout<<"<<<<< iter: "<<iter<<"\n";
    // define the scale ratio
    scale_coeff = 0.25 + 1.75 * (1.- AMIN1((iter-1)/500., 1.));

    // calcuate VP mass, mass_center, energy and pressure
    Get_VP_info_and_VD_generation(world);

    // check the partition error
    Check_partition_error(iter, world);

    if (convergenced == 0){

      // update the VP positon
      Update_position(world);

      // update the particle position according to the centroidal Voronoi tessellation
      Update_position_CVT(world);
    }

    // if (iter%1 == 0){
    //  Output_vtk_file(sph->i_iter, iter, world);
//      Output_pov_ray_file(sph->i_iter,iter, world);
//      Output_plt_file(sph->i_iter,iter, world);
//      Output_ghost_vtk_file(sph->i_iter,iter, world);
//      sph->Output_vtk_file(iter, 2, world);
    //  sph->Output_plt_file(iter, 2, world);
    // }

    if (iter%600 == 0)
      relax_ratio = 1.02*relax_ratio;

    relax_ratio = AMIN1(relax_ratio, 0.90);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Partition iter "<<iter<<" finished\n";
  out.close();
#endif
  }

  sph->time_for_partition[world.rank()] += sph->time_partition.elapsed();
  sph->time_graph.restart();

  // construct graph for migrating SPH particles
  Construct_graph(world);

  sph->time_for_graph[world.rank()] += sph->time_graph.elapsed();
  if (world.rank() == 0){
    cout<<"<<<<< CVP partitioning finished\n";
    cout<<"**********************************************************\n";
  }
}
//-------------------------------------------------------
// output VP coord
//-------------------------------------------------------
void Voronoi::Output_vtk_file(int i_iter, int n, communicator &world)
{
  if (world.rank() == 0){
    FILE    *fp;
    char    filename[256];
    sprintf(filename,"%s%d%s%d%s","./outdata/voro_particle.",i_iter,".",n,".vtk");
  
    ofstream out(filename, ios::trunc);
    out<<"# vtk DataFile Version 4.0\nvtk output t: "<<n<<"\n";
    out<<"ASCII\nDATASET POLYDATA\n";

    out<<"POINTS "<<num_partition<<" double\n";
    for (int i = 0; i < num_partition; i++)
      out<<vo_particle[i]->coord.i<<" "
         <<vo_particle[i]->coord.j<<" "
         <<vo_particle[i]->coord.k<<"\n";

    out<<"VERTICES "<<num_partition<<" "<<2*num_partition<<"\n";
    for (int i = 0; i < num_partition; i++)
      out<<"1 "<<i<<"\n";

    out<<"POINT_DATA "<<num_partition<<"\n";
    out<<"FIELD FieldData 7\n";

    //color
    out<<"Index 1 "<<num_partition<<" int\n";
    for (int i = 0; i < num_partition; i++)
      out<<vo_particle[i]->id<<"\n";

    //id
    out<<"Color 1 "<<num_partition<<" int\n";
    for (int i = 0; i < num_partition; i++)
      out<<vo_particle[i]->color<<"\n";

    //mass
    out<<"Mass 1 "<<num_partition<<" double\n";
    for (int i = 0; i < num_partition; i++)
      out<<vo_particle[i]->mass<<"\n";

    //P
    out<<"P 1 "<<num_partition<<" double\n";
    for (int i = 0; i < num_partition; i++)
      out<<vo_particle[i]->P<<"\n";

    //h
    out<<"h 1 "<<num_partition<<" double\n";
    for (int i = 0; i < num_partition; i++)
      out<<vo_particle[i]->h<<"\n";

    //error
    out<<"error 1 "<<num_partition<<" double\n";
    for (int i = 0; i < num_partition; i++)
      out<<vo_particle[i]->error<<"\n";

    //f
    out<<"force 3 "<<num_partition<<" double\n";
    for (int i = 0; i < num_partition; i++)
      out<<vo_particle[i]->f.i<<" "<<vo_particle[i]->f.j<<" "<<vo_particle[i]->f.k<<"\n";
    out.close();
  }
}
//--------------------------------------------------
// Output VP infomation
//--------------------------------------------------
void Voronoi::Output_plt_file(int i_iter, int n, communicator &world)
{
  if (world.rank() == 0){
    FILE    *fp;
    char    filename[256];
    sprintf(filename,"%s%d%s%d%s","./outdata/voro_particle.", i_iter,".",n,".plt");

    ofstream out(filename, ios::trunc);
  out<<"VARIABLES = \"x\",\"y\",\"z\",\"id\",\"color\",\"m\",\"p\",\"h\",\"error\"\n";
    for (int i = 0; i < num_partition; i++){
      Vo_particle *current_particle = vo_particle[i];
      out<<current_particle->coord.i<<" "<<current_particle->coord.j<<"  "<<current_particle->coord.k<<"  "
         <<current_particle->id<<" "<<current_particle->color<<" "
         <<current_particle->mass<<"  "
         <<current_particle->P<<"  "
         <<current_particle->h<<"  "
         <<current_particle->error<<"\n";
    }
    out.close();
  }
}
//-------------------------------------------------------
// output pov ray file
//-------------------------------------------------------
void Voronoi::Output_pov_ray_file(int i_iter, int n, communicator &world)
{
  if (world.rank() == 0){
    double x,y,z;

    // Set up the number of blocks that the container is divided into
    int n_x, n_y, n_z;

    particle_order po;

    // Create a pre-container class to import the input file and guess the
    // best computational grid size to use.
    voro::pre_container pcon(x_min,x_max,y_min,y_max,z_min,z_max,P_X,P_Y,P_Z);

    wall_plane plane1(0., 0., 1., fabs(z_max), num_partition + 1);
    wall_plane plane2(0., 1., 0., fabs(y_max), num_partition + 2);
    wall_plane plane3(1., 0., 0., fabs(x_max), num_partition + 3);
    wall_plane plane4(0., 0., -1., fabs(z_min), num_partition + 4);
    wall_plane plane5(0., -1., 0., fabs(y_min), num_partition + 5);
    wall_plane plane6(-1., 0., 0., fabs(x_min), num_partition + 6);

    for(int i=0;i<num_partition;i++){
      x = y = z = 0.;
#if P_DIM_X == 1
      x = vo_particle[i]->coord.i;
#endif
#if P_DIM_Y == 1
      y = vo_particle[i]->coord.j;
#endif
#if P_DIM_Z == 1
      z = vo_particle[i]->coord.k;
#endif
      pcon.put(i,x,y,z);
    }

    pcon.guess_optimal(n_x,n_y,n_z);

    // Set up the container class and import the particles from the
    // pre-container
    voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,P_X,P_Y,P_Z,8);

    // Add wall to the boundary
#if P_DIM_Z == 1 && PERI_DIM_Z == 0
    con.add_wall(plane1);
    con.add_wall(plane4);
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 0
    con.add_wall(plane2);
    con.add_wall(plane5);
#endif
#if P_DIM_X == 1 && PERI_DIM_X == 0
    con.add_wall(plane3);
    con.add_wall(plane6);
#endif

    // add particles into the container

    for(int i=0;i<num_partition;i++){
      x = y = z = 0.;
#if P_DIM_X == 1
      x = vo_particle[i]->coord.i;
#endif
#if P_DIM_Y == 1
      y = vo_particle[i]->coord.j;
#endif
#if P_DIM_Z == 1
      z = vo_particle[i]->coord.k;
#endif
      con.put(po,i,x,y,z);
    }

    int id;
    voronoicell_neighbor c;
    vector<int> neigh,f_vert;
    vector<double> v;
    c_loop_all cl(con);

    char    filename1[256];
    sprintf(filename1,"%s%d%s%d%s","./tpout/polygons.",i_iter,".",n,".pov");
    char    filename2[256];
    sprintf(filename2,"%s%d%s%d%s","./tpout/particles.",i_iter,".",n,".pov");

    visual.Set_scene(filename1,1);
    visual.Set_scene(filename2,1);

#ifndef _NON_PERI_

    visual.Draw_boundary(filename1, x_min,x_max,y_min,y_max,z_min,z_max);
    visual.Draw_boundary(filename2, x_min,x_max,y_min,y_max,z_min,z_max);

#endif

    if(cl.start()) do if(con.compute_cell(c,cl)) {

      // Get the position of the current particle under consideration
      cl.pos(x,y,z);
      id=cl.pid();

      // Gather information about the computed Voronoi cell
      c.neighbors(neigh);
      c.face_vertices(f_vert);
      c.vertices(x,y,z,v);

      int     i, j;
      Real    rad = 0.;
      rad = vo_particle[id]->mass/glbl_total_mass*domain.i*1.5;
      // iteration along neighbors

      for(i=0, j=0;i<neigh.size();i++) {
        if(neigh[i]>=0) {
          visual.Draw_polygon(filename1, f_vert, v,j, Real(id), 0., Real(num_partition), 1);
          visual.Draw_polygon(filename2, f_vert, v,j, Real(id), 0., Real(num_partition), 2);
        }
        j+=f_vert[j]+1;
      }
      visual.Draw_particle(filename2, vo_particle[id]->coord, vo_particle[id]->id,
                           rad, 0., Real(num_partition), 0.95);

    } while (cl.inc());
    con.clear();
  }
}
//-------------------------------------------------------
// output VP coord
//-------------------------------------------------------
void Voronoi::Output_ghost_vtk_file(int i_iter, int n, communicator &world)
{
  if (world.rank() == 0){
    FILE    *fp;
    char    filename[256];
    sprintf(filename,"%s%d%s%d%s","./outdata/voro_ghost_particle.",i_iter,".",n,".vtk");
  
    ofstream out(filename, ios::trunc);
    out<<"# vtk DataFile Version 4.0\nvtk output t: "<<n<<"\n";
    out<<"ASCII\nDATASET POLYDATA\n";

    out<<"POINTS "<<vo_particle_ghost.size()<<" double\n";
    for (int i = 0; i < vo_particle_ghost.size(); i++)
      out<<vo_particle_ghost[i]->coord.i<<" "
         <<vo_particle_ghost[i]->coord.j<<" "
         <<vo_particle_ghost[i]->coord.k<<"\n";

    out<<"VERTICES "<<vo_particle_ghost.size()<<" "<<2*vo_particle_ghost.size()<<"\n";
    for (int i = 0; i < vo_particle_ghost.size(); i++)
      out<<"1 "<<i<<"\n";

    out<<"POINT_DATA "<<vo_particle_ghost.size()<<"\n";
    out<<"FIELD FieldData 5\n";

    //color
    out<<"Index 1 "<<vo_particle_ghost.size()<<" int\n";
    for (int i = 0; i < vo_particle_ghost.size(); i++)
      out<<vo_particle_ghost[i]->id + vo_particle.size()<<"\n";

    //id
    out<<"Color 1 "<<vo_particle_ghost.size()<<" int\n";
    for (int i = 0; i < vo_particle_ghost.size(); i++)
      out<<vo_particle_ghost[i]->color<<"\n";

    //mass
    out<<"Mass 1 "<<vo_particle_ghost.size()<<" double\n";
    for (int i = 0; i < vo_particle_ghost.size(); i++)
      out<<vo_particle_ghost[i]->mass<<"\n";

    //P
    out<<"P 1 "<<vo_particle_ghost.size()<<" double\n";
    for (int i = 0; i < vo_particle_ghost.size(); i++)
      out<<vo_particle_ghost[i]->P<<"\n";

    //h
    out<<"h 1 "<<vo_particle_ghost.size()<<" double\n";
    for (int i = 0; i < vo_particle_ghost.size(); i++)
      out<<vo_particle_ghost[i]->h<<"\n";

    out.close();
  }
}
#if PERI_DIM != 0
//-------------------------------------------------------
// shift SPH particle position according to
// the new VD
//-------------------------------------------------------
void Voronoi::Shift_SPH_particle_position(SPH *sph, communicator &world)
{

  cross_boundary.i = 0;
  cross_boundary.j = 0;
  cross_boundary.k = 0;

  vector<my_real> neighbor_v;
  neighbor_v.clear();

  sph_particle = sph->particle;

  double x,y,z;

  // Set up the number of blocks that the container is divided into
  int n_x, n_y, n_z;

  particle_order po;

  // Create a pre-container class to import the input file and guess the
  // best computational grid size to use.
  voro::pre_container pcon(x_min,x_max,y_min,y_max,z_min,z_max,P_X,P_Y,P_Z);

  wall_plane plane1(0., 0., 1., fabs(z_max), num_partition + 1);
  wall_plane plane2(0., 1., 0., fabs(y_max), num_partition + 2);
  wall_plane plane3(1., 0., 0., fabs(x_max), num_partition + 3);
  wall_plane plane4(0., 0., -1., fabs(z_min), num_partition + 4);
  wall_plane plane5(0., -1., 0., fabs(y_min), num_partition + 5);
  wall_plane plane6(-1., 0., 0., fabs(x_min), num_partition + 6);

  for(int i=0;i<num_partition;i++){
    x = y = z = 0.;
#if P_DIM_X == 1
    x = vo_particle[i]->coord.i;
#endif
#if P_DIM_Y == 1
    y = vo_particle[i]->coord.j;
#endif
#if P_DIM_Z == 1
    z = vo_particle[i]->coord.k;
#endif
    pcon.put(i,x,y,z);
  }

  pcon.guess_optimal(n_x,n_y,n_z);

  // Set up the container class and import the particles from the
  // pre-container
  voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,P_X,P_Y,P_Z,8);
  // Add wall to the boundary
#if P_DIM_Z == 1 && PERI_DIM_Z == 0
  con.add_wall(plane1);
  con.add_wall(plane4);
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 0
  con.add_wall(plane2);
  con.add_wall(plane5);
#endif
#if P_DIM_X == 1 && PERI_DIM_X == 0
  con.add_wall(plane3);
  con.add_wall(plane6);
#endif

  // add particles into the container
  for(int i=0;i<num_partition;i++){
    x = y = z = 0.;
#if P_DIM_X == 1
    x = vo_particle[i]->coord.i;
#endif
#if P_DIM_Y == 1
    y = vo_particle[i]->coord.j;
#endif
#if P_DIM_Z == 1
    z = vo_particle[i]->coord.k;
#endif
    con.put(po,i,x,y,z);
  }

  int id;
  voronoicell_neighbor c;
  vector<int> neigh,f_vert;
  vector<double> v;
  vector<double> normal;
  c_loop_all cl(con);

  if(cl.start()) do if(con.compute_cell(c,cl)) {

    // Get the position of the current particle under consideration
    cl.pos(x,y,z);
    id=cl.pid();
    if (id == world.rank()){

      // Gather information about the computed Voronoi cell
      c.neighbors(neigh);
      c.face_vertices(f_vert);
      c.vertices(x,y,z,v);
      c.normals(normal);

      my_real p_coord;
      p_coord.i = x; p_coord.j = y; p_coord.k = z;
      neighbor_v.push_back(p_coord);

      my_real current_coord;
      my_real neighbor_coord;
      my_set_data (current_coord, p_coord);

      // iteration along neighbors
      int j = 0;
      for(int i=0; i<neigh.size(); i++) {
        if(neigh[i]>=0){

          int l;
          l=3*f_vert[j+1];
          neighbor_coord.i = v[l];
          neighbor_coord.j = v[l+1];
          neighbor_coord.k = v[l+2];
          my_real       dx = my_minus_data (neighbor_coord, current_coord);
          my_real dr;
          dr.i = normal[3*i+0];
          dr.j = normal[3*i+1];
          dr.k = normal[3*i+2];
          Real  dist = fabs(get_dot_product (dx, dr)) * 2.;
          dr   = my_multiply_const (dr, dist);
          neighbor_coord = my_add_data (current_coord, dr);
          neighbor_v.push_back(neighbor_coord);
          int n=f_vert[j];
          // iterate all the vertex of the cell boundary
          for(int k=0; k<n; k++) {
            my_real v_temp;
            int ll=3*f_vert[j+k+1];
#if P_DIM_X == 1 && PERI_DIM_X == 1
            if (v[ll] < box_l.i ) cross_boundary.i = -1;
            else if (v[ll] > box_r.i) cross_boundary.i = 1;
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 1
            if (v[ll+1] < box_l.j ) cross_boundary.j = -1;
            else if (v[ll+1] > box_r.j) cross_boundary.j = 1;
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 1
            if (v[ll+2] < box_l.k ) cross_boundary.k = -1;
            else if (v[ll+2] > box_r.k) cross_boundary.k = 1;
#endif
          }
        }
        j+=f_vert[j]+1;
      }
    }
  } while (cl.inc());
  con.clear();

  Real need_for_shift = 0;
  my_self_mold(cross_boundary, need_for_shift);
  if (need_for_shift > 1.e-10){
    int n_neighbor = neighbor_v.size();
    for (int i = 0; i < n_neighbor; i++){
      if (cross_boundary.i == -1) neighbor_v[i].i += domain.i;
      if (cross_boundary.j == -1) neighbor_v[i].j += domain.j;
      if (cross_boundary.k == -1) neighbor_v[i].k += domain.k;
    }
    static affinity_partitioner ap;
    parallel_for( blocked_range<int>(0, sph->total_num_particle),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        p_Particle cur_sph_particle = sph_particle[i];
        my_real p_coord, p_coord_0;
        p_coord.i = P_DIM_X==1 ? cur_sph_particle->coord.i : 0.;
        p_coord.j = P_DIM_Y==1 ? cur_sph_particle->coord.j : 0.;
        p_coord.k = P_DIM_Z==1 ? cur_sph_particle->coord.k : 0.;
        my_set_data (p_coord_0, p_coord);

        my_int shift;
        int index_record = -1;
        shift.i = shift.j = shift.k = 0;
        for (int m = 0; m <= 1 && index_record != 0; m++)
          for (int s = 0; s <= 1 && index_record != 0; s++)
            for (int t = 0; t <= 1 && index_record != 0; t++){
              index_record = -1;
              my_set_data (p_coord, p_coord_0);
              p_coord.i += domain.i*m*P_DIM_X;
              p_coord.j += domain.j*s*P_DIM_Y;
              p_coord.k += domain.k*t*P_DIM_Z;

              Real   dist_min = 1.e20;
              for(int l=0; l<n_neighbor; l++){
                Real   dist  = get_distance_2p (p_coord, neighbor_v[l]);
                if (dist < dist_min){
                  dist_min = dist;
                  index_record = l;
                }
              }
            }
        if (index_record == 0){
          cur_sph_particle->coord.i = P_DIM_X==1 ? p_coord.i : cur_sph_particle->coord.i;
          cur_sph_particle->coord.j = P_DIM_Y==1 ? p_coord.j : cur_sph_particle->coord.j;
          cur_sph_particle->coord.k = P_DIM_Z==1 ? p_coord.k : cur_sph_particle->coord.k;
        }
      }
    }, ap);
  }
}
//-------------------------------------------------------
// shift back SPH particle position
//-------------------------------------------------------
void Voronoi::Shift_back_SPH_particle_position(SPH *sph, communicator &world)
{
  sph_particle = sph->particle;
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, sph->total_num_particle),
           [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      p_Particle cur_sph_particle = sph_particle[i];
#if P_DIM_X == 1 && PERI_DIM_X == 1
      if (cur_sph_particle->coord.i > box_r.i+1.e-20)
        cur_sph_particle->coord.i -= domain.i;
      else if (cur_sph_particle->coord.i < box_l.i-1.e-20)
        cur_sph_particle->coord.i += domain.i;
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 1
      if (cur_sph_particle->coord.j > box_r.j+1.e-20)
        cur_sph_particle->coord.j -= domain.j;
      else if (cur_sph_particle->coord.j < box_l.j-1.e-20)
        cur_sph_particle->coord.j += domain.j;
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 1
      if (cur_sph_particle->coord.k > box_r.k+1.e-20)
        cur_sph_particle->coord.k -= domain.k;
      else if (cur_sph_particle->coord.k < box_l.k-1.e-20)
        cur_sph_particle->coord.k += domain.k;
#endif
    }
  }, ap);
}
#endif

