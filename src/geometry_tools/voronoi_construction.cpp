#include "voronoi_construction.h"
#ifdef _MESH_GENERATION_
//-------------------------------------------------------
// local VD construction for a set of SPH particles
//-------------------------------------------------------
void Voronoi_construction::Local_VD_generation(Particle *sp, Graphmeshcls *mesh)
{
  static affinity_partitioner ap;

  // Set up the number of blocks that the container is divided into
  int n_x, n_y, n_z;
  double x,y,z;

  particle_order po;

  // set an approximate boundary according to the particle scale for local VD
  Real x_min,x_max,y_min,y_max,z_min,z_max;
  my_real coord = sp->coord;
  Real    scale = sp->h*CUT_OFF*1.5;
  
  x_min= DIM_X ? coord.i - scale : -0.5;
  x_max= DIM_X ? coord.i + scale : 0.5;
  y_min= DIM_Y ? coord.j - scale : -0.5;
  y_max= DIM_Y ? coord.j + scale : 0.5;
  z_min= DIM_Z ? coord.k - scale : -0.5;
  z_max= DIM_Z ? coord.k + scale : 0.5;

  int num_neighbor = sp->neighbor.size();
  
  // Create a pre-container class to import the input file and guess the
  // best computational grid size to use.
  voro::pre_container pcon(x_min,x_max,y_min,y_max,z_min,z_max,false,false,false);

  wall_plane plane1(0.,  0.,  1., fabs(z_max), num_neighbor + 1);
  wall_plane plane2(0.,  1.,  0., fabs(y_max), num_neighbor + 2);
  wall_plane plane3(1.,  0.,  0., fabs(x_max), num_neighbor + 3);
  wall_plane plane4(0.,  0., -1., fabs(z_min), num_neighbor + 4);
  wall_plane plane5(0., -1.,  0., fabs(y_min), num_neighbor + 5);
  wall_plane plane6(-1., 0.,  0., fabs(x_min), num_neighbor + 6);

  for(int i=0;i<num_neighbor;i++){
    x = y = z = 0.;
#if DIM_X == 1
    x = sp->neighbor[i]->coord.i;
#endif
#if DIM_Y == 1
    y = sp->neighbor[i]->coord.j;
#endif
#if DIM_Z == 1
    z = sp->neighbor[i]->coord.k;
#endif
    pcon.put(i,x,y,z);
  }

  pcon.guess_optimal(n_x,n_y,n_z);

  // Set up the container class and import the particles from the
  // pre-container
  voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,false,false,false,8);

  // Add wall to the boundary
#if DIM_Z == 1
  con.add_wall(plane1);
  con.add_wall(plane4);
#endif
#if DIM_Y == 1
  con.add_wall(plane2);
  con.add_wall(plane5);
#endif
#if DIM_X == 1
  con.add_wall(plane3);
  con.add_wall(plane6);
#endif

  // add particles into the container
  for(int i=0;i<num_neighbor;i++){
    x = y = z = 0.;
#if DIM_X == 1
    x = sp->neighbor[i]->coord.i;
#endif
#if DIM_Y == 1
    y = sp->neighbor[i]->coord.j;
#endif
#if DIM_Z == 1
    z = sp->neighbor[i]->coord.k;
#endif
    con.put(po,i,x,y,z);
  }

  int id;
  voronoicell_neighbor c;
  vector<int> neigh;
  // vector<int> f_vert;
  // vector<double> v;
  // vector<double> face;
  // vector<double> normal;
  c_loop_all cl(con);

  // TODO: do a small test here
  if(cl.start()) do if(con.compute_cell(c,cl) && sp->neighbor[cl.pid()]->local_id == sp->local_id) {
    // Get the position of the current particle under consideration
    cl.pos(x,y,z);
    id=cl.pid();

    // Gather information about the computed Voronoi cell
    c.neighbors(neigh);
    // iteration along neighbors
    int num_neighbor_cell = neigh.size();
    
    // if (DIM ==2 && num_neighbor_cell <= 5 && num_neighbor_cell == num_neighbor)
    //   cout<<"The recorded triangulation point is wrongh !!!"<<" "<<num_neighbor_cell<<" "<<sp->neighbor.size()<<endl;

    for(int i=0;i<num_neighbor_cell;i++) {

      // Get the neighbor vp coord, considering BC
      if (neigh[i] >= 0 && neigh[i] < num_neighbor) { // if the boundary is open, neigh[i]<0
        //Real dist = get_distance_2p (sp->coord, sp->neighbor[neigh[i]]->coord);
        //if (dist < sp->h)
          mesh->Add_edge(sp->local_id, sp->neighbor[neigh[i]]->local_id);
      }else{
//         cout<<"<<<<< ERROR!:Not a locally compact cell!!!"<<endl;
      }
    }
  } while (cl.inc());

  con.clear();
}
#endif