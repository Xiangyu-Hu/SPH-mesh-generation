#ifdef _MESH_GENERATION_
// #include "sph_mesh_generation.h"
#include "particle_mesh_generation.h"
#include "level_set.h"
#include "lset_level_infor.h"
#include "lset_package.h"
#include "lset_cell.h"
#endif
#include "voronoi_lset.h"

/*********************************************************/
/*                                                       */
/*  Functions defined in class "Voronoi_lset"            */
/*                                                       */
/*********************************************************/
// JZ20190111::1st version, only for initial particle distribution

//-------------------------------------------------------
// Get_coord_of_VP_at_iRank
//-------------------------------------------------------
my_real Voronoi_lset::Get_coord_of_VP_at_iRank(int iRank)
{
  if (iRank >= int(vo_particle.size())){
    cout<<"<<<<< ERROR iRank value in getting VP position"<<endl;
    my_real rr; my_set_const (rr, 0.);
    return rr;
  }else
    return vo_particle[iRank]->coord;
}
//-------------------------------------------------------
// random the particle distribution
//-------------------------------------------------------
void Voronoi_lset::Random_particle_distribution(Levelset *level_set, communicator &world)
{
  // serialization_vector <p_Vo_particle> vp_info;
  // vp_info.Vector.clear();
  // vp_info.mem_size = int(vo_particle.size());
  // vp_info.tag      = 0;

  if (world.rank() == 0){
    // get the random number
    srand((unsigned)time(NULL));

    static affinity_partitioner ap;
    parallel_for( blocked_range<int>(0, int(vo_particle.size())),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        my_real position;
        bool    positive = false;

        Real dx = level_set->lset_level_info[0]->dl*5.;

        while (!positive){
          #if (P_DIM_X)
          position.i = rand()/ double(RAND_MAX+1.0);
          vo_particle[i]->coord.i = position.i*particle_box.i + rand()/ double(RAND_MAX+1.0)*1.e-3 + particle_box_l.i;
          vo_particle[i]->coord.i = AMAX1(particle_box_l.i+dx,AMIN1(vo_particle[i]->coord.i, particle_box_r.i-dx));
          #endif
          #if (P_DIM_Y)
          position.j = rand()/ double(RAND_MAX+1.0);
          vo_particle[i]->coord.j = position.j*particle_box.j + rand()/ double(RAND_MAX+1.0)*1.e-3 + particle_box_l.j;
          vo_particle[i]->coord.j = AMAX1(particle_box_l.j+dx,AMIN1(vo_particle[i]->coord.j, particle_box_r.j-dx));
          #endif
          #if (P_DIM_Z)
          position.k = rand()/ double(RAND_MAX+1.0);
          vo_particle[i]->coord.k = position.k*particle_box.k + rand()/ double(RAND_MAX+1.0)*1.e-3 + particle_box_l.k;
          vo_particle[i]->coord.k = AMAX1(particle_box_l.k+dx,AMIN1(vo_particle[i]->coord.k, particle_box_r.k-dx));
          #endif
          Real phi = level_set->Get_phi_at_position(vo_particle[i]->coord);
          if (phi > dx) positive = true;
        }
      }
    }, ap);

    // for (int i = 0; i < vo_particle.size(); i++){
    //   vo_particle[i]->tag = 1;
    //   vp_info.Vector.push_back(vo_particle[i]);
    // }
    // vp_info.tag      = 1;
    // broadcast(world,vp_info,0);

  }
  // else{

  //   broadcast(world,vp_info,0);
  //   for (int i = 0; i < vo_particle.size(); i++){
  //     if (vo_particle[i]->id != vp_info.Vector[i]->id){
  //       cout<<"VP id is wrong!!!\n"; world.abort(-1);
  //     }
  //     my_set_data (vo_particle[i]->coord, vp_info.Vector[i]->coord);
  //   }
  // }
  // if (world.rank() == 0){
  //   vp_info.Vector.clear();
  // }else{
  //   static affinity_partitioner ap;
  //   parallel_for( blocked_range<concurrent_vector<p_Vo_particle>::iterator>(vp_info.Vector.begin(), vp_info.Vector.end()),
  //       [&](const blocked_range<concurrent_vector<p_Vo_particle>::iterator>& r){
  //     for(concurrent_vector<p_Vo_particle>::iterator it=r.begin(); it!=r.end(); ++it){
  //       if (NULL != *it){
  //         delete (*it); // release memory of exchange_particles in heap
  //         *it = NULL;
  //       }
  //     }  
  //   }, ap);
  //   vp_info.Vector.clear();
  //   vp_info.Vector.shrink_to_fit();
  //   if (int(vp_info.Vector.capacity()) != 0){
  //     cout<<"Vector memory is not released!!!\n"; world.abort(-1);
  //   }
  // }
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
void Voronoi_lset::Initialize(Levelset *level_set, int flag, communicator &world)
{
  num_partition       = world.size();
  #if defined (_CVP_LSET_INIT_)
  num_partition       = level_set->num_color_tmp;
  #endif
  convergenced        = 0;
  Max_iterate_number  = 800;
  error_tolerance     = level_set->error_tolerance;
  scale_coeff         = 0.;
  my_set_data (domain        , level_set->domain);
  my_set_data (box_l         , level_set->box_l);
  my_set_data (box_r         , level_set->box_r);
  my_set_data (particle_box  , level_set->domain);
  my_set_data (particle_box_l, level_set->box_l);
  my_set_data (particle_box_r, level_set->box_r);
  my_set_data (dpkg          , level_set->lset_level_info[0]->dpkg);
  relax_ratio           = 0.80;
  total_energy          = 0.;
  total_interface_area  = 0.;
  error_max             = 0.;
  error_record          = 0.;

  contain_cutcell.clear();

  if (flag == 1){
    total_num_positive_pkg = level_set->Get_positive_pkg (lset_positive_pkg);

    _vo_particle = new p_Vo_particle[num_partition];
    // generate the partitioning particles
      for(int i=0; i!=num_partition; ++i){
        _vo_particle[i] = Vo_particlepool.malloc();
        _vo_particle[i]->Initialize(i);
        vo_particle.push_back(_vo_particle[i]);
      }

    Random_particle_distribution (level_set, world);

    static affinity_partitioner ap;
    parallel_for( blocked_range<int>(0, total_num_positive_pkg),
             [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        p_Levelset_package cur_pkg = lset_positive_pkg[i];
        my_real coord_pkg = cur_pkg->get_pkg_position(cur_pkg->index.i, cur_pkg->index.j, cur_pkg->index.k, dpkg, box_l);
        cur_pkg->color = Get_CPU_id(coord_pkg);
      }
    }, ap);
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

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Class Voronoi_lset initialized\n";
  out.close();
#endif
}
//-------------------------------------------------------
// get the cpu index of certain position
//-------------------------------------------------------
int Voronoi_lset::Get_CPU_id(my_real p_coord)
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
void Voronoi_lset::Set_target_mass(Levelset *level_set, communicator &world)
{
  glbl_total_mass = level_set->total_mass;

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
// Get_vp_position_and_hmin
//-------------------------------------------------------
void Voronoi_lset::Get_vp_position_and_hmin (std::vector<my_real> & vp_coords, std::vector<Real> &vp_scale, communicator &world)
{
  for (int i = 0; i < num_partition; ++i){
    vp_coords.push_back(vo_particle[i]->coord);
    vp_scale.push_back(vo_particle[i]->h_min);
  }

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Get_vp_position_and_hmin finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// Ghost_VP_generation
//-------------------------------------------------------
void Voronoi_lset::Ghost_VP_generation(Levelset *level_set, communicator &world)
{
  tbb::atomic<bool> contains_cutcell[num_partition];
  for (int i = 0; i < num_partition; ++i)
    contains_cutcell[i] = false;

  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, total_num_positive_pkg),
           [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      p_Levelset_package cur_pkg = lset_positive_pkg[i];
      if (cur_pkg->total_mass_surface > 0. && cur_pkg->color > 0) contains_cutcell[cur_pkg->color] = true;
    }
  }, ap);

  for (int i = 0; i < num_partition; ++i){
    if(contains_cutcell[i] == true) contain_cutcell.push_back(i);
  }

  int num_ghost = vo_particle_ghost.size();
  for (int i = 0; i < num_ghost; ++i){
    vo_particle_ghost[i]->Reset();
  }
  vo_particle_ghost.clear();

  int num_lset_bound = contain_cutcell.size();
  int index = 0;
  for (int i = 0; i < num_lset_bound; ++i)
  {
    p_Vo_particle mir_vp = vo_particle[contain_cutcell[i]];

    Particle temp;
    my_set_data (temp.coord, mir_vp->coord);
    temp.Calculate_particle_infor(level_set);

    if (temp.phi >= 1.e-6){
      p_Vo_particle new_vp = Vo_particlepool.malloc();
      new_vp->Initialize(num_partition+index);
      index ++;

      Real n_x = temp.norm.i;
      Real n_y = temp.norm.j;
      Real n_z = temp.norm.k;

      new_vp->coord.i = P_DIM_X ? mir_vp->coord.i - 2.*temp.phi*n_x : mir_vp->coord.i;
      new_vp->coord.j = P_DIM_Y ? mir_vp->coord.j - 2.*temp.phi*n_y : mir_vp->coord.j;
      new_vp->coord.k = P_DIM_Z ? mir_vp->coord.k - 2.*temp.phi*n_z : mir_vp->coord.k;

      vo_particle_ghost.push_back(new_vp);

    }else{
      cout<<"<<<<< WARNING: current VP is not in the positive region"<<endl;
    }
  }

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Ghost_VP_generation finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// calcuate VP mass, mass_center, energy and pressure
//-------------------------------------------------------
void Voronoi_lset::Get_VP_info_and_VD_generation(communicator &world)
{ 
  // char    filename2[256];
  // sprintf(filename2,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  // ofstream out2(filename2, ios::app);

  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, num_partition),
           [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      vo_particle[i]->Reset();
    }
  }, ap);

  int num_ghost = vo_particle_ghost.size();
  for (int i = 0; i < num_ghost; ++i){
    vo_particle_ghost[i]->Reset();
  }

  // Set up the number of blocks that the container is divided into
  int n_x, n_y, n_z;
  double x,y,z;

  particle_order po;

  // Create a pre-container class to import the input file and guess the
  // best computational grid size to use.
  voro::pre_container pcon(x_min,x_max,y_min,y_max,z_min,z_max,P_X,P_Y,P_Z);

  wall_plane plane1( 0.,  0.,  1., fabs(z_max), num_partition + num_ghost + 1);
  wall_plane plane2( 0.,  1.,  0., fabs(y_max), num_partition + num_ghost + 2);
  wall_plane plane3( 1.,  0.,  0., fabs(x_max), num_partition + num_ghost + 3);
  wall_plane plane4( 0.,  0., -1., fabs(z_min), num_partition + num_ghost + 4);
  wall_plane plane5( 0., -1.,  0., fabs(y_min), num_partition + num_ghost + 5);
  wall_plane plane6(-1.,  0.,  0., fabs(x_min), num_partition + num_ghost + 6);

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

  for(int i=0;i<num_ghost;i++){
    x = y = z = 0.;
#if P_DIM_X == 1
    x = vo_particle_ghost[i]->coord.i;
#endif
#if P_DIM_Y == 1
    y = vo_particle_ghost[i]->coord.j;
#endif
#if P_DIM_Z == 1
    z = vo_particle_ghost[i]->coord.k;
#endif
    pcon.put(i+num_partition,x,y,z);
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

  for(int i=0;i<num_ghost;i++){
    x = y = z = 0.;
#if P_DIM_X == 1
    x = vo_particle_ghost[i]->coord.i;
#endif
#if P_DIM_Y == 1
    y = vo_particle_ghost[i]->coord.j;
#endif
#if P_DIM_Z == 1
    z = vo_particle_ghost[i]->coord.k;
#endif
    con.put(po,i+num_partition,x,y,z);
  }
  
  // parallel_for( blocked_range<int>(0, total_num_positive_pkg),
  //          [&](const blocked_range<int>& r){
  //   for(int i=r.begin(); i!=r.end(); ++i){
    for(int i=0; i!=total_num_positive_pkg; ++i){
      p_Levelset_package cur_pkg = lset_positive_pkg[i];
      int index = 0;
      Real rx, ry, rz, xx, yy, zz;
      my_real coord_pkg = cur_pkg->get_pkg_position(cur_pkg->index.i, cur_pkg->index.j, cur_pkg->index.k, dpkg, box_l);
      xx = P_DIM_X==1 ? coord_pkg.i : 0;
      yy = P_DIM_Y==1 ? coord_pkg.j : 0;
      zz = P_DIM_Z==1 ? coord_pkg.k : 0;
      con.find_voronoi_cell(xx,yy,zz,rx,ry,rz,index);
      if (index > num_partition || index < 0)
        cur_pkg->color = -1;
      else
        cur_pkg->color = index;
    }
  // }, ap);
  
  parallel_for( blocked_range<int>(0, num_partition),
           [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      p_Vo_particle cur_vp = vo_particle[i];
      for (int ip = 0; ip < total_num_positive_pkg; ip++){
        p_Levelset_package cur_pkg = lset_positive_pkg[ip];
        my_real          coord_pkg = cur_pkg->get_pkg_position(cur_pkg->index.i, cur_pkg->index.j, cur_pkg->index.k, dpkg, box_l);
        if (cur_pkg->color == cur_vp->color){
          cur_vp->mass        += cur_pkg->total_mass;
          my_real        temp  = my_multiply_const(coord_pkg, cur_pkg->total_mass);
          cur_vp->mass_center  = my_add_data(cur_vp->mass_center, temp);
          Real           dist  = get_distance_2p (cur_vp->coord, coord_pkg);
          cur_vp->energy      += cur_pkg->total_mass * dist;
        }
      }
    }
  }, ap);
  
  // serialization_vector <p_Vo_particle> exchange_vector;
  // exchange_vector.Vector.clear();

  if (world.rank() == 0){
    // serialization_vector <p_Vo_particle> *gather_exchange_vector;
    // gather_exchange_vector = new serialization_vector <p_Vo_particle>[world.size()];

    // parallel_for( blocked_range<int>(0, world.size()),
    //          [&](const blocked_range<int>& r){
    //   for(int i=r.begin(); i!=r.end(); ++i){
    //     gather_exchange_vector[i].Vector.clear();
    //     gather_exchange_vector[i].mem_size = 0;
    //     gather_exchange_vector[i].tag = 0;
    //   }
    // }, ap);

    // gather(world,exchange_vector,gather_exchange_vector,0);

    // parallel_for( blocked_range<int>(0, num_partition),
    //          [&](const blocked_range<int>& r){
    //   for(int i=r.begin(); i!=r.end(); ++i){
    //     for (int j = 1; j < world.size(); j++){
    //       if (vo_particle[i]->id != gather_exchange_vector[j].Vector[i]->id){
    //         cout<<"VP id is wrong!!!\n"; world.abort(-1);
    //       }
    //       vo_particle[i]->mass        +=  gather_exchange_vector[j].Vector[i]->mass;
    //       vo_particle[i]->mass_center  =  my_add_data(vo_particle[i]->mass_center,
    //                                       gather_exchange_vector[j].Vector[i]->mass_center);
    //       vo_particle[i]->energy      +=  gather_exchange_vector[j].Vector[i]->energy;
    //     }
    //   }
    // }, ap);

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

    error_max = 0.;
    total_energy   = 0.;
    for(int i=0; i!=num_partition; ++i){
      p_Vo_particle cur_vp = vo_particle[i];
      error_max      = AMAX1(error_max, cur_vp->error);
      total_energy  += cur_vp->energy;
    }

    // for( int i=1; i<world.size(); i++){
    //   parallel_for( blocked_range<concurrent_vector<p_Vo_particle>::iterator>
    //               (gather_exchange_vector[i].Vector.begin(), gather_exchange_vector[i].Vector.end()),
    //            [&](const blocked_range<concurrent_vector<p_Vo_particle>::iterator>& r){
    //     for(concurrent_vector<p_Vo_particle>::iterator it=r.begin(); it!=r.end(); ++it){
    //       if (NULL != *it){
    //         delete (*it); // release memory of exchange_particles in heap
    //         *it = NULL;
    //       }
    //     }
    //   }, ap);
    //   gather_exchange_vector[i].Vector.clear();
    //   gather_exchange_vector[i].Vector.shrink_to_fit();
    //   if (int(gather_exchange_vector[i].Vector.capacity()) != 0){
    //     cout<<"Vector memory is not released!!!\n"; world.abort(-1);
    //   }
    // }
    // delete [] gather_exchange_vector;

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

    if(cl.start()) do if(con.compute_cell(c,cl)) {
      // Get the position of the current particle under consideration
      cl.pos(x,y,z);
      id=cl.pid();

      if(id < num_partition){
        // Gather information about the computed Voronoi cell
        c.neighbors(neigh);
        c.face_areas(face);
        c.normals(normal);
        c.face_vertices(f_vert);
        c.vertices(x,y,z,v);

        total_interface_area += c.surface_area();

        int     i, j;
        int     num_neighbor = 0;
        Real    scale_record = 0.;
        Real    h_min_record = domain.i+domain.j+domain.k;
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
            }else if (neigh[i] > num_partition){ // neighbor is outside of the domain or ghosts
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
            if (neigh[i] < num_partition) h_min_record = AMIN1 (h_min_record, 0.5*dist);
            Real pressure_middle = 0.5*(current_P + neighbor_P);

            vo_particle[id]->f = my_add_data (vo_particle[id]->f, my_multiply_const(dr, -1*pressure_middle*face[i]));
          }
          j+=f_vert[j]+1;
        }

        vo_particle[id]->h     = scale_coeff*scale_record/Real(num_neighbor);
        vo_particle[id]->h_min = h_min_record;
      }
    } while (cl.inc());
  }
  // else{
  //   exchange_vector.mem_size = 0;
  //   exchange_vector.tag = 0;
  
  //   for( int i=0; i < num_partition; i++){
  //     vo_particle[i]->tag = 2;
  //     exchange_vector.Vector.push_back(vo_particle[i]);
  //   }
  //   exchange_vector.mem_size = int(exchange_vector.Vector.size());
  //   exchange_vector.tag = 1;
  //   gather(world,exchange_vector,0);
  // }
  con.clear();
  // exchange_vector.Vector.clear();
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
void Voronoi_lset::Update_position(communicator &world)
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
void Voronoi_lset::Update_position_CVT(communicator &world)
{
  // serialization_vector <p_Vo_particle> vp_info;
  // vp_info.Vector.clear();
  // vp_info.mem_size = int(vo_particle.size());
  // vp_info.tag      = 0;
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
      }
    }, ap);

    // for (int i = 0; i < vo_particle.size(); i++){
    //   vo_particle[i]->tag = 1;
    //   vp_info.Vector.push_back(vo_particle[i]);
    // }
    // vp_info.tag      = 1;
    // broadcast(world,vp_info,0);

  }
  // else{
  
  //   broadcast(world,vp_info,0);
  //   for (int i = 0; i < vo_particle.size(); i++){
  //     if (vo_particle[i]->id != vp_info.Vector[i]->id){
  //       cout<<"VP id is wrong!!!\n"; world.abort(-1);
  //     }
  //     my_set_data (vo_particle[i]->coord, vp_info.Vector[i]->coord);
  //   }
  // }

  // if (world.rank() == 0){
  //   vp_info.Vector.clear();
  // }else{
  //   static affinity_partitioner ap;
  //   parallel_for( blocked_range<concurrent_vector<p_Vo_particle>::iterator>(vp_info.Vector.begin(), vp_info.Vector.end()),
  //       [&](const blocked_range<concurrent_vector<p_Vo_particle>::iterator>& r){
  //     for(concurrent_vector<p_Vo_particle>::iterator it=r.begin(); it!=r.end(); ++it){
  //       if (NULL != *it){
  //         delete (*it); // release memory of exchange_particles in heap
  //         *it = NULL;
  //       }
  //     }  
  //   }, ap);
  //   vp_info.Vector.clear();
  //   vp_info.Vector.shrink_to_fit();
  //   if (int(vp_info.Vector.capacity()) != 0){
  //     cout<<"Vector memory is not released!!!\n"; world.abort(-1);
  //   }
  // }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Update_position_CVT finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// check the partition error
//-------------------------------------------------------
void Voronoi_lset::Check_partition_error(int &iter, communicator &world)
{
  if (world.rank() == 0){
    if (iter%10 == 0){
      cout<<"Iter: "<<iter<<" max_error: "<<error_max
          <<" Total energy: "<<total_energy
          <<" Interface area: "<<total_interface_area<<"\n";
    }
    if (iter%10 == 0){
      ofstream save("outdata/CVP_Lset_Error_history.dat", ios::app);
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
    // broadcast(world, convergenced, 0);
  }
  // else{
  //   broadcast(world, convergenced, 0);
  // }
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Check_partition_error finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// CVP Partitioning for lset
//-------------------------------------------------------
void Voronoi_lset::Partitioning_for_initial_VP_distribution(Levelset *level_set, communicator &world)
{
  Initialize(level_set, 2, world);

  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< CVP for level set is processing\n";
  }

  Set_target_mass(level_set, world);

  for(int iter=0; iter<=Max_iterate_number && convergenced == 0; iter++){

//    if (world.rank() == 0)
//      cout<<"<<<<< iter: "<<iter<<"\n";
    // define the scale ratio
    scale_coeff = 0.25 + 1.75 * (1.- AMIN1((iter-1)/500., 1.));

    // generating ghost particles
    Ghost_VP_generation(level_set, world);

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

    if (iter%100 == 0){
     Output_vtk_file(0, iter, world);
     Output_ghost_vtk_file(0,iter, world);
     level_set->n_out+=iter;
     level_set->Output_level_set_vti(level_set->n_out, world);
    }

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

  level_set->Output_level_set_vti(level_set->n_out++, world);

  if (world.rank() == 0){
    cout<<"<<<<< CVP partitioning for level_set finished\n";
    cout<<"**********************************************************\n";
  }
}
//-------------------------------------------------------
// output VP coord
//-------------------------------------------------------
void Voronoi_lset::Output_vtk_file(int i_iter, int n, communicator &world)
{
  if (world.rank() == 0){
    FILE    *fp;
    char    filename[256];
    sprintf(filename,"%s%d%s%d%s","./outdata/voro_lset_particle.",i_iter,".",n,".vtk");
  
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
//-------------------------------------------------------
// output VP coord
//-------------------------------------------------------
void Voronoi_lset::Output_ghost_vtk_file(int i_iter, int n, communicator &world)
{
  if (world.rank() == 0){
    FILE    *fp;
    char    filename[256];
    sprintf(filename,"%s%d%s%d%s","./outdata/voro_lset_ghost_particle.",i_iter,".",n,".vtk");
  
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
