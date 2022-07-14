#include "boost/multi_array.hpp"
#include "glbfunc.h"
#include "graph_mesh.h"
#ifdef _MESH_GENERATION_
#include "level_set.h"
#include "lset_level_infor.h"
#include "particle_mesh_generation.h"
#include "sph_mesh_generation.h"
#endif

/***************************************************/
/*                                                 */
/*    Functions defined in class "Graphmeshcls"    */
/*                                                 */
/***************************************************/

//-------------------------------------------------------
// Initialization
//-------------------------------------------------------
void Graphmeshcls::Initialize(int nvert, communicator &world)
{
  for(int i=0; i<nvert; i++)
    add_vertex(i,graph);

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<class Graphcls is initialized\n";
  out.close();
#endif
}
#if defined (_MESH_GENERATION_) 
//-------------------------------------------------------
// Reconstruct::local triangulation
//-------------------------------------------------------
void Graphmeshcls::Reconstruct_tris(std::vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  tetgenmesh m;
  // reconstruct the element tris
  tris.clear();
//   tris.shrink_to_fit();

  // reset the output
  int tag_error = 0;
  
  int iRank = world.rank();

  // edges to be deleted
  std::set<int2_graph> edge_delete;

  //iterate the graph
  std::pair<EdgeIterator, EdgeIterator> ei = edges(graph);
  for (EdgeIterator eit = ei.first; eit != ei.second; ++eit)
  { 
    EdgeDescriptor ed = *eit;

    int id_src = int(graph[source(ed,graph)]);
    int id_targ = int(graph[target(ed,graph)]);
    int id_third = -1;
    std::pair<OutEdgeIterator, OutEdgeIterator> src_ei = out_edges(source(ed,graph), graph);
    std::vector<VertexDescriptor> vertex_record;
    vertex_record.clear();
    for (OutEdgeIterator src_eit = src_ei.first; src_eit != src_ei.second; ++src_eit)
    { 
      EdgeDescriptor src = *src_eit;
      std::pair<OutEdgeIterator, OutEdgeIterator> targ_ei = out_edges(target(ed,graph), graph);
      for (OutEdgeIterator targ_eit = targ_ei.first; targ_eit != targ_ei.second; ++targ_eit)
      {             
              EdgeDescriptor targ = *targ_eit;
        if(int(graph[target(src,graph)]) == int(graph[target(targ,graph)]))
        {
                 vertex_record.push_back(target(src,graph));
          id_third = int(graph[target(src,graph)]);
          // construct the tri
          int3_graph tri;
          tri.i = AMAX1(id_src, AMAX1(id_targ, id_third));
          tri.k = AMIN1(id_src, AMIN1(id_targ, id_third));
          if(id_src>tri.k && id_src<tri.i)
                  tri.j = id_src;
          if(id_targ>tri.k && id_targ<tri.i)
                  tri.j = id_targ;
          if(id_third>tri.k && id_third<tri.i)
                  tri.j = id_third;

          if(tri.i == tri.j || tri.i == tri.k || tri.j == tri.k)
                  cout<<"<<<<< ["<<world.rank()<<"] ERROR!!! The particle coincides in pair !!!"<<endl;

          #ifdef _MPI_
          int count_color = 0;
          int color_list[3] = {iRank, iRank, iRank};
          if (particle_total[tri.i]->color != iRank ) {count_color ++; color_list[0] = particle_total[tri.i]->color;}
          if (particle_total[tri.j]->color != iRank ) {count_color ++; color_list[1] = particle_total[tri.j]->color;}
          if (particle_total[tri.k]->color != iRank ) {count_color ++; color_list[2] = particle_total[tri.k]->color;}
          
          if(count_color == 3)
          {
            if (color_list[0] == color_list[1] &&
                color_list[0] == color_list[2] &&
                color_list[1] == color_list[2]){
              int2_graph edge_id;
              edge_id.i      = tri.i;
              edge_id.j      = tri.j;
              edge_delete.insert(edge_id);
              
              edge_id.i      = tri.i;
              edge_id.j      = tri.k;
              edge_delete.insert(edge_id);
              
              edge_id.i      = tri.j;
              edge_id.j      = tri.k;
              edge_delete.insert(edge_id);
            }
          }else if (count_color == 2){
            if (color_list[0] == iRank){
              if (color_list[1] == color_list[2]){
                int2_graph edge_id;
                edge_id.i      = tri.j;
                edge_id.j      = tri.k;
                edge_delete.insert(edge_id);
              }else{
                tris.insert(tri);
              }
            }else if (color_list[1] == iRank){
              if (color_list[0] == color_list[2]){
                int2_graph edge_id;
                edge_id.i      = tri.i;
                edge_id.j      = tri.k;
                edge_delete.insert(edge_id);
              }else{
                tris.insert(tri);
              }
            }else if (color_list[2] == iRank){
              if (color_list[0] == color_list[1]){
                int2_graph edge_id;
                edge_id.i      = tri.i;
                edge_id.j      = tri.j;
                edge_delete.insert(edge_id);
              }else{
                tris.insert(tri);
              }
            }
          }else{
          #endif
            if( particle_total[tri.i]->type != REAL_PARTICLE && 
                particle_total[tri.j]->type != REAL_PARTICLE && 
                particle_total[tri.k]->type != REAL_PARTICLE)
            {
                // judge something
              p_Particle particle_1 = particle_total[tri.i];
              p_Particle particle_2 = particle_total[tri.j];
              p_Particle particle_3 = particle_total[tri.k];
                // construct particle temporally
              Particle particle_temp;
              particle_temp.coord.i = 1./3.*(particle_1->coord.i + particle_2->coord.i + particle_3->coord.i);
              particle_temp.coord.j = 1./3.*(particle_1->coord.j + particle_2->coord.j + particle_3->coord.j);
              particle_temp.coord.k = 1./3.*(particle_1->coord.k + particle_2->coord.k + particle_3->coord.k);
              particle_temp.Calculate_particle_infor(sph);

              Real phi_size         = 1./3.*(particle_1->phi + particle_2->phi + particle_3->phi);

              if(particle_temp.phi > 0.0001*phi_size)
                tris.insert(tri);
              else{

                int2_graph edge_id;
                Real dist;
                my_real vertex_0, vertex_1, vertex_2;

                  // obtain the coordinates
                my_set_data (vertex_0, particle_1->coord);
                my_set_data (vertex_1, particle_2->coord);
                my_set_data (vertex_2, particle_3->coord);

                  // obtain the vectors
                my_real vector_1 = my_minus_data (vertex_1, vertex_0);
                dist             = get_distance  (vector_1);
                edge_id.i        = tri.i;
                edge_id.j        = tri.j;

                my_real vector_2 = my_minus_data (vertex_2, vertex_1);
                Real      dist_2 = get_distance  (vector_2);
                if(dist_2 > dist)
                {
                  dist           = dist_2;
                  edge_id.i      = tri.j;
                  edge_id.j      = tri.k;
                }

                my_real vector_3 = my_minus_data (vertex_0, vertex_2);
                Real      dist_3 = get_distance  (vector_3);
                if(dist_3 > dist)
                {
                  dist           = dist_3;
                  edge_id.i      = tri.i;
                  edge_id.j      = tri.k;
                }

                edge_delete.insert(edge_id);
              }
            }else{
              tris.insert(tri);
            }
          #ifdef _MPI_
          }
          #endif
        }
      }
    }
    // if (DIM == 3){
    //   if (int(vertex_record.size()) >= 2){
    //     my_real vec_21 = my_minus_data (particle_total[id_targ]->coord, particle_total[id_src]->coord);
    //     int vr_size = int(vertex_record.size());
    //     for (int i = 0; i < vr_size-1; i++){
    //       my_real vec_31 = my_minus_data (particle_total[int(graph[vertex_record[i]])]->coord, particle_total[id_src]->coord);
    //       my_real norm   = get_cross_product (vec_21, vec_31);
    //       norm = my_multiply_const (norm, 1./(get_distance(norm)+1.e-15));
    //       for (int j = i+1; j < vr_size; j++){
    //         my_real vec_41 = my_minus_data (particle_total[int(graph[vertex_record[j]])]->coord, particle_total[id_src]->coord);
    //         Real dist = fabs(get_dot_product (norm, vec_41));
    //         // cout<<"dist "<<dist<<endl;

    //         if (dist <= 0.25*sph->level_set.minimum_dl){
    //           std::pair<EdgeDescriptor, bool> test_1 = edge(graph[vertex_record[i]],graph[vertex_record[j]],graph);
    //           std::pair<EdgeDescriptor, bool> test_2 = edge(graph[vertex_record[i]],graph[vertex_record[j]],graph);
    //           if (test_1.second == true || test_2.second == true) {
    //             Real dist_21 = get_distance    (vec_21);
    //             Real dist_31 = get_distance    (vec_31);
    //             Real dist_41 = get_distance    (vec_41);
    //             Real dist_32 = get_distance_2p (particle_total[int(graph[vertex_record[i]])]->coord, particle_total[id_targ]->coord);
    //             Real dist_42 = get_distance_2p (particle_total[int(graph[vertex_record[j]])]->coord, particle_total[id_targ]->coord);
    //             Real dist_43 = get_distance_2p (particle_total[int(graph[vertex_record[j]])]->coord, particle_total[int(graph[vertex_record[i]])]->coord);

    //             Real max_edge = AMAX1(dist_21, AMAX1(dist_31, AMAX1(dist_41,AMAX1(dist_32,AMAX1(dist_42, dist_43)))));

    //             int2_graph edge_id;
    //             if      (fabs(max_edge-dist_21) < 1.e-10){edge_id.i = AMAX1(id_targ, id_src); edge_id.j = AMIN1(id_targ, id_src);}
    //             else if (fabs(max_edge-dist_31) < 1.e-10){edge_id.i = AMAX1(int(graph[vertex_record[i]]), id_src); edge_id.j = AMIN1(int(graph[vertex_record[i]]), id_src);}
    //             else if (fabs(max_edge-dist_41) < 1.e-10){edge_id.i = AMAX1(int(graph[vertex_record[j]]), id_src); edge_id.j = AMIN1(int(graph[vertex_record[j]]), id_src);}
    //             else if (fabs(max_edge-dist_32) < 1.e-10){edge_id.i = AMAX1(int(graph[vertex_record[i]]), id_targ); edge_id.j = AMIN1(int(graph[vertex_record[i]]), id_targ);}
    //             else if (fabs(max_edge-dist_42) < 1.e-10){edge_id.i = AMAX1(int(graph[vertex_record[j]]), id_targ); edge_id.j = AMIN1(int(graph[vertex_record[j]]), id_targ);}
    //             else if (fabs(max_edge-dist_43) < 1.e-10){edge_id.i = AMAX1(int(graph[vertex_record[j]]), int(graph[vertex_record[i]])); edge_id.j = AMIN1(int(graph[vertex_record[j]]), int(graph[vertex_record[i]]));}

    //             edge_delete.insert(edge_id);
    //             // cout<<"<<<<< edge: "<<edge_id.i<<" "<<edge_id.j<<"is deleted (dist: "<<dist<<" )"<<endl;
    //           }
    //         }
    //       }
    //     }
    //   }
    // }
  }

  if(int(edge_delete.size()) != 0)
    cout<<"<<<<< ["<<world.rank()<<"] The number of edges to be deleted : "<<int(edge_delete.size())<<endl;

  // remove all the edges
  Remove_all_edges (world);

  std::set <         int3_graph>::iterator        it;
  // std::pair<std::set<int3_graph>::iterator,bool>  ret;
  std::set <         int2_graph>::iterator        edge_delete_iter;

  for (it=tris.begin(); it!=tris.end(); ++it)
  {
    int vertex_0, vertex_1, vertex_2;

    // obtain the coordinates
    vertex_0 = (*it).i;
    vertex_1 = (*it).j;
    vertex_2 = (*it).k;

    int flag = 0;
    for(edge_delete_iter=edge_delete.begin(); edge_delete_iter!=edge_delete.end(); ++edge_delete_iter)
    {
      if(vertex_0 == (*edge_delete_iter).i && vertex_1 == (*edge_delete_iter).j)
        flag = 1;
      if(vertex_0 == (*edge_delete_iter).i && vertex_2 == (*edge_delete_iter).j)
        flag = 1;
      if(vertex_1 == (*edge_delete_iter).i && vertex_2 == (*edge_delete_iter).j)
        flag = 1;
    }

    if(flag == 0){
      Add_edge(vertex_0, vertex_1);
      Add_edge(vertex_1, vertex_2);
      Add_edge(vertex_0, vertex_2);
    }
  }

  // reset tri elements
  std::set<int3_graph> tris_delete;
  tris.clear();
  tris_delete.clear();

  // reconstruct the mesh graph for the second time
  // iterate the graph
  ei = edges(graph);
  std::vector <p_Particle> error_particle;
  for (EdgeIterator eit = ei.first; eit != ei.second; ++eit)
  { 
    EdgeDescriptor ed = *eit; 

    int id_src = int(graph[source(ed,graph)]);
    int id_targ = int(graph[target(ed,graph)]);
    int id_third;
    std::pair<OutEdgeIterator, OutEdgeIterator> src_ei = out_edges(source(ed,graph), graph);
    std::vector<VertexDescriptor> vertex_record;
    vertex_record.clear();
    for (OutEdgeIterator src_eit = src_ei.first; src_eit != src_ei.second; ++src_eit)
    { 
      EdgeDescriptor src = *src_eit;
      std::pair<OutEdgeIterator, OutEdgeIterator> targ_ei = out_edges(target(ed,graph), graph);
      for (OutEdgeIterator targ_eit = targ_ei.first; targ_eit != targ_ei.second; ++targ_eit)
      {             
        EdgeDescriptor targ = *targ_eit;
        if(int(graph[target(src,graph)]) == int(graph[target(targ,graph)]))
        {
          id_third = int(graph[target(src,graph)]);
          vertex_record.push_back(target(src,graph));
          // construct the tri
          int3_graph tri;
          tri.i = AMAX1(id_src, AMAX1(id_targ, id_third));
          tri.k = AMIN1(id_src, AMIN1(id_targ, id_third));
          if(id_src>tri.k && id_src<tri.i)
            tri.j = id_src;
          if(id_targ>tri.k && id_targ<tri.i)
            tri.j = id_targ;
          if(id_third>tri.k && id_third<tri.i)
            tri.j = id_third;
          tri.vol = 0.;
          tri.aspect = 0.;

          tris.insert(tri);

          if (DIM == 2){
            if(int(vertex_record.size()) > 2)
            {
              // cout<<"Something is wrong about the graph !!! "<<endl;

              error_particle.push_back(particle_total[id_src]);
              error_particle.push_back(particle_total[id_targ]);
              error_particle.push_back(particle_total[id_third]);

              tag_error = 1;
              //exit(0);
            }else if(int(vertex_record.size())== 2)
            {
              std::pair<EdgeDescriptor, bool> test_1 = edge(graph[vertex_record[0]],graph[vertex_record[1]],graph);
              std::pair<EdgeDescriptor, bool> test_2 = edge(graph[vertex_record[1]],graph[vertex_record[0]],graph);
              if (test_1.second == true || test_2.second == true) {
                cout<<"<<<<< ["<<world.rank()<<"] Some edges in the graph intersect !!! "<<endl;

                error_particle.push_back(particle_total[id_src]);
                error_particle.push_back(particle_total[id_targ]);
                error_particle.push_back(particle_total[id_third]);

                tag_error = 1;
                //exit(0);
              }
            }
          }
        }
      }
    }
  }
  if (error_particle.size() > 0){
    cout<<"<<<<< ["<<world.rank()<<"] Num. of error particle for tris: "<<error_particle.size()<<endl;
  }

  //iterate the graph
  std::pair<VertexIterator, VertexIterator> vertice_ei = vertices(graph);
  for (VertexIterator eit = vertice_ei.first; eit != vertice_ei.second; ++eit)
  { 
    VertexDescriptor vt = *eit;
    int id_vt = int(graph[vt]);

    if(particle_total[id_vt]->type == SINGULARITY_PARTICLE)
    {
      if(out_degree(vt, graph) < 2)
      {
        cout<<"<<<< ["<<world.rank()<<"] The characteristic particles are not involved in the final mesh."<<endl;
        tag_error = 1;
        //exit(0);
      }
    }
  }

  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Mesh quality section"<<endl;
  }
  
  // calculate mesh qualities
  Tri_quality_statistics (particle_total, sph, world);
  
  if (world.rank() == 0)
    cout<<"**********************************************************\n";

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Reconstruct finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// Tri_quality_statistics
//-------------------------------------------------------
void Graphmeshcls::Tri_quality_statistics(std::vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  Real dangle = 4.0;
  int num_angletable = int(180./dangle);

  int _theta_30 = 0;
  int *_angle;
  int _num_tris = int(tris.size());

  int local_communication_tris = 0;

  int  _ntet_sm_10 = 0;
  int  _ntet_sm_20 = 0;
  int  _ntet_sm_30 = 0;
  int  _ntet_sm_40 = 0;

  _angle = new int [num_angletable];

  for(int i=0; i<num_angletable; i++)
    _angle[i] = 0;
  
  Real _maximum_angle            = 0.;
  Real _minimum_angle            = 1.e15;
  Real _maximum_tri_aspect_ratio = 0.;
  Real _minimum_angle_avg        = 0.;
  Real _maximum_edge             = 0.;
  Real _minimum_edge             = 1.e15;
  Real _maximum_area             = 0.;
  Real _minimum_area             = 1.e15;
  Real _total_vol                = 0.;
  Real _G_min                    = 1.e20;
  Real _G_avg                    = 0.;

  std::set<int3_graph>::iterator it;
  for (it=tris.begin(); it!=tris.end(); ++it)
  {
    my_real vertex_0, vertex_1, vertex_2;

    // obtain the coordinates
    my_set_data (vertex_0, particle_total[int((*it).i)]->coord);
    my_set_data (vertex_1, particle_total[int((*it).j)]->coord);
    my_set_data (vertex_2, particle_total[int((*it).k)]->coord);

    std::vector<int> color_list;
    int color_count = 1;

    color_list.push_back(particle_total[int((*it).i)]->color);
    color_list.push_back(particle_total[int((*it).j)]->color);
    color_list.push_back(particle_total[int((*it).k)]->color);

    std::sort (color_list.begin(), color_list.end(), Greater());

    for (int i = 0; i < 2; ++i)
    {
      if(color_list[i] != color_list[i+1]) color_count++;
    }
    if (color_count > 1) local_communication_tris++;

    // obtain the vectors
    my_real vector_1 = my_minus_data (vertex_1, vertex_0);
    Real      dist_1 = get_distance  (vector_1);

    my_real vector_2 = my_minus_data (vertex_2, vertex_1);
    Real      dist_2 = get_distance  (vector_2);

    my_real vector_3 = my_minus_data (vertex_0, vertex_2);
    Real      dist_3 = get_distance  (vector_3);

    Real dist_maximum = AMAX1(dist_1, AMAX1(dist_2, dist_3));
    Real dist_minimum = AMIN1(dist_1, AMIN1(dist_2, dist_3));

    _maximum_edge = AMAX1(_maximum_edge, dist_maximum);
    _minimum_edge = AMIN1(_minimum_edge, dist_minimum);

    _maximum_tri_aspect_ratio = AMAX1(_maximum_tri_aspect_ratio, dist_maximum/(dist_minimum + 1.e-20));
    (*it).aspect = dist_maximum/(dist_minimum + 1.e-20);

    Real alpha_1 = 180./acos(-1.) * acos(-(vector_1.i * vector_3.i + vector_1.j * vector_3.j)/(dist_3*dist_1 + 1.e-20));
    Real alpha_2 = 180./acos(-1.) * acos(-(vector_1.i * vector_2.i + vector_1.j * vector_2.j)/(dist_2*dist_1 + 1.e-20));
    Real alpha_3 = 180./acos(-1.) * acos(-(vector_2.i * vector_3.i + vector_2.j * vector_3.j)/(dist_3*dist_2 + 1.e-20));

    alpha_1 = AMAX1(0., AMIN1(alpha_1, 180.));
    alpha_2 = AMAX1(0., AMIN1(alpha_2, 180.));
    alpha_3 = AMAX1(0., AMIN1(alpha_3, 180.));

    _maximum_angle = AMAX1(_maximum_angle, AMAX1(alpha_1, AMAX1(alpha_2, alpha_3)));
    _minimum_angle = AMIN1(_minimum_angle, AMIN1(alpha_1, AMIN1(alpha_2, alpha_3)));

    _angle[int(alpha_1/dangle)] += 1;
    _angle[int(alpha_2/dangle)] += 1;
    _angle[int(alpha_3/dangle)] += 1;

    if(AMIN1(alpha_1, AMIN1(alpha_2, alpha_3)) < 30.)
      _theta_30++;

    Real minmin = AMIN1(alpha_1, AMIN1(alpha_2, alpha_3));
    if (minmin < 10.) _ntet_sm_10++;
    if (minmin < 20.) _ntet_sm_20++;
    if (minmin < 30.) _ntet_sm_30++;
    if (minmin < 40.) _ntet_sm_40++;

    _minimum_angle_avg +=  AMIN1(alpha_1, AMIN1(alpha_2, alpha_3));

    Real area = 0.5*(dist_3*dist_1)*sin(acos(-(vector_1.i * vector_3.i + vector_1.j * vector_3.j)/(dist_3*dist_1  + 1.e-20)));

    (*it).vol = area;
    _total_vol += area;

    _maximum_area = AMAX1(_maximum_area, area);
    _minimum_area = AMIN1(_minimum_area, area);

    Real G = 2. * sqrt(3.) * area / (dist_maximum * (dist_1 + dist_2 + dist_3)/2.);

    _G_min = AMIN1(G, _G_min);
    _G_avg = _G_avg + G;
  }
  
  //JZ20181230::adding for MPI implementation
  if (world.rank() == 0){
    int  theta_30      = 0;
    int  glbl_num_edge = 0;
    Real total_vol     = 0.;

    int *angle;
    angle = new int [num_angletable];
    for(int i=0; i<num_angletable; i++)
      angle[i] = 0;

    glbl_num_tris             = 0;
    maximum_angle             = 0.;
    minimum_angle             = 0.;
    maximum_tri_aspect_ratio  = 0.;
    maximum_edge              = 0.;
    minimum_edge              = 0.;
    maximum_area              = 0.;
    minimum_area              = 0.;
    G_min                     = 0.;
    minimum_angle_avg         = 0.;
    G_avg                     = 0.;
    glbl_communication_volume = 0.;
    glbl_communication_tris   = 0;
    
    reduce( world, _maximum_angle,            maximum_angle,            mpi::maximum<Real>(), 0);
    reduce( world, _minimum_angle,            minimum_angle,            mpi::minimum<Real>(), 0);
    reduce( world, _maximum_tri_aspect_ratio, maximum_tri_aspect_ratio, mpi::maximum<Real>(), 0);
    reduce( world, _maximum_edge,             maximum_edge,             mpi::maximum<Real>(), 0);
    reduce( world, _minimum_edge,             minimum_edge,             mpi::minimum<Real>(), 0);
    reduce( world, _maximum_area,             maximum_area,             mpi::maximum<Real>(), 0);
    reduce( world, _minimum_area,             minimum_area,             mpi::minimum<Real>(), 0);
    reduce( world, _G_min,                    G_min,                    mpi::minimum<Real>(), 0);
    reduce( world, _minimum_angle_avg,        minimum_angle_avg,        std::plus<Real>()   , 0);
    reduce( world, _total_vol,                total_vol,                std::plus<Real>()   , 0);
    reduce( world, _G_avg,                    G_avg,                    std::plus<Real>()   , 0);
    reduce( world, _num_tris,                 glbl_num_tris,            std::plus<int>()    , 0);
    reduce( world, local_communication_tris,  glbl_communication_tris,  std::plus<int>()    , 0);
    reduce( world, int(num_edges(graph)),     glbl_num_edge,            std::plus<int>()    , 0);
    reduce( world, _theta_30,                 theta_30,                 std::plus<int>()    , 0);
    reduce( world,    _angle, num_angletable, angle,                    std::plus<int>()    , 0);
    reduce( world, _ntet_sm_10     ,          ntet_sm_10   ,            std::plus<int>()          , 0);
    reduce( world, _ntet_sm_20     ,          ntet_sm_20   ,            std::plus<int>()          , 0);
    reduce( world, _ntet_sm_30     ,          ntet_sm_30   ,            std::plus<int>()          , 0);
    reduce( world, _ntet_sm_40     ,          ntet_sm_40   ,            std::plus<int>()          , 0);

    G_avg /= glbl_num_tris;
    minimum_angle_avg /= glbl_num_tris;
    glbl_communication_volume = Real(glbl_communication_tris)/Real(glbl_num_tris);
    
    ntri_sm_10 = int(angle[0]);
    ntri_sm_20 = int(angle[1])+ntri_sm_10;
    ntri_sm_30 = int(angle[2])+ntri_sm_20;
    ntri_sm_40 = int(angle[3])+ntri_sm_30;
    
    cout.setf(ios::fixed);
    cout<<"  By current algorithms, the total edges is                    : "<<glbl_num_edge<<endl;
    cout<<"  By current algorithms, the total triangles is                : "<<glbl_num_tris<<endl;
    cout<<"  By current algorithms, the maximum aspect ratio              : "<<setprecision(6)<<maximum_tri_aspect_ratio<<endl;
    cout<<"  By current algorithms, the maximum and minimum angle         : "<<setprecision(6)<<maximum_angle<<" "<<minimum_angle<<endl;
    cout<<"  By current algorithms, the maximum and minimum edge          : "<<setprecision(6)<<maximum_edge<<" "<<minimum_edge<<endl;
    cout<<"  By current algorithms, the maximum and minimum area          : "<<setprecision(6)<<maximum_area<<" "<<minimum_area<<endl;
    cout<<"  By current algorithms, the average and minimum quality       : "<<setprecision(6)<<G_avg<<" "<<G_min<<endl;
    cout<<"  By current algorithms, the minimum angle less than 30 degree : "<<setprecision(6)<<Real(theta_30)/Real(glbl_num_tris)<<endl;
    cout<<"  By current algorithms, the minimum angle average degree      : "<<setprecision(6)<<minimum_angle_avg<<endl;
    cout<<endl;
    
    int total_angle = 0;
    // for(int i=0; i<18; i++)
    // {
    //   cout<<"  the angle between "<<i*10<<" - "<<(i+1)*10<<" degree is: "<<setprecision(6)<<int(angle[i])<<" with percentage "<<Real(angle[i])/(3.*Real(glbl_num_tris))<<endl;
    //   total_angle += angle[i];
    // }

    cout<<" The angle table:"<<endl;
    for (int i = 0; i < num_angletable/2; ++i)
    {
      Real st1, st2, ed1, ed2, pp1, pp2;
      unsigned long ang1, ang2;
      st1 = Real(i                   )*dangle;
      st2 = Real(i  +num_angletable/2)*dangle;
      ed1 = Real(i+1                 )*dangle;
      ed2 = Real(i+1+num_angletable/2)*dangle;
      ang1 = angle[i];
      ang2 = angle[i+num_angletable/2];
      pp1  = Real(ang1)/Real(glbl_num_tris*3)*100.;
      pp2  = Real(ang2)/Real(glbl_num_tris*3)*100.;

      printf("     %4.4g - %4.4g degrees:  %8ld (%03.2f %%)   |    %4.4g - %4.4g degrees:  %8ld (%03.2f %%) \n",
              st1, ed1, ang1, pp1, st2, ed2, ang2, pp2);
      total_angle += (ang1+ang2);
    }

    if(total_angle/glbl_num_tris != 3)
    {
      cout<<"<<<<<ERROR!!! The angle calculation is wrong !!!"<<endl;
      // exit(-1);
    }

    char    ang_histo[256];
    sprintf(ang_histo,"%s%d%s","./outdata/angle_histogram_step_",int(sph->run_time),".csv");
    ofstream ang_histo_out(ang_histo, ios::trunc);

    ang_histo_out<<"\"degree\",\"number\",\"percentage\"\n";
    for (int i = 0; i < num_angletable; ++i)
    {
      Real           st1 = Real(i+1)*dangle;
      unsigned long ang1 = angle[i];
      Real           pp1 = Real(ang1)/Real(glbl_num_tris*3)*100.;
      ang_histo_out<<st1<<","<<ang1<<","<<pp1<<endl;
    }
    ang_histo_out.close();

  }else{
    reduce( world, _maximum_angle,            mpi::maximum<Real>(), 0);
    reduce( world, _minimum_angle,            mpi::minimum<Real>(), 0);
    reduce( world, _maximum_tri_aspect_ratio, mpi::maximum<Real>(), 0);
    reduce( world, _maximum_edge,             mpi::maximum<Real>(), 0);
    reduce( world, _minimum_edge,             mpi::minimum<Real>(), 0);
    reduce( world, _maximum_area,             mpi::maximum<Real>(), 0);
    reduce( world, _minimum_area,             mpi::minimum<Real>(), 0);
    reduce( world, _G_min,                    mpi::minimum<Real>(), 0);
    reduce( world, _minimum_angle_avg,        std::plus<Real>()   , 0);
    reduce( world, _total_vol,                std::plus<Real>()   , 0);
    reduce( world, _G_avg,                    std::plus<Real>()   , 0);
    reduce( world, _num_tris,                 std::plus<int>()    , 0);
    reduce( world, local_communication_tris,  std::plus<int>()    , 0);
    reduce( world, int(num_edges(graph)),     std::plus<int>()    , 0);
    reduce( world, _theta_30,                 std::plus<int>()    , 0);
    reduce( world, _angle,    num_angletable, std::plus<int>()    , 0);
    reduce( world, _ntet_sm_10     ,          std::plus<int>()    , 0);
    reduce( world, _ntet_sm_20     ,          std::plus<int>()    , 0);
    reduce( world, _ntet_sm_30     ,          std::plus<int>()    , 0);
    reduce( world, _ntet_sm_40     ,          std::plus<int>()    , 0);
  }
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Tri_quality_statistics finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// Check_whether_this_tet_should_be_delete
//-------------------------------------------------------
bool Graphmeshcls::Check_whether_this_tet_should_be_delete(int stage, p_Particle particle_1, p_Particle particle_2, p_Particle particle_3, p_Particle particle_4, int iRank, SOLVER *sph)
{
  int count_color = 0;

  if (particle_1->type == DUMMY_PARTICLE) return true;
  if (particle_2->type == DUMMY_PARTICLE) return true;
  if (particle_3->type == DUMMY_PARTICLE) return true;
  if (particle_4->type == DUMMY_PARTICLE) return true;

  #ifdef _MPI_
  int color_list[4] = {iRank, iRank, iRank, iRank};
  if (particle_1->color != iRank ) {count_color ++; color_list[0] = particle_1->color;}
  if (particle_2->color != iRank ) {count_color ++; color_list[1] = particle_2->color;}
  if (particle_3->color != iRank ) {count_color ++; color_list[2] = particle_3->color;}
  if (particle_4->color != iRank ) {count_color ++; color_list[3] = particle_4->color;}

  if (count_color == 4){
    return true;
  }else if (count_color > 0){
    for (int i = 0; i < 4; i++){
      if (color_list[i] > iRank){
       return true;
      }
    }
  }
  #endif
    // Particle particle_temp1;
    // particle_temp1.coord.i = 1./4.*(particle_1->coord.i + particle_2->coord.i + particle_3->coord.i + particle_4->coord.i);
    // particle_temp1.coord.j = 1./4.*(particle_1->coord.j + particle_2->coord.j + particle_3->coord.j + particle_4->coord.j);
    // particle_temp1.coord.k = 1./4.*(particle_1->coord.k + particle_2->coord.k + particle_3->coord.k + particle_4->coord.k);
    // particle_temp1.Calculate_particle_infor(sph);

    my_real  box_l = sph->level_set.box_l;
    my_real  box_r = sph->level_set.box_r;
    Real        dl = sph->level_set.lset_level_info[0]->dl;

    Particle particle_temp1;
    Real s1=get_distance(get_cross_product(my_minus_data(particle_3->coord, particle_2->coord), my_minus_data(particle_4->coord, particle_2->coord)))/2.0; 
    Real s2=get_distance(get_cross_product(my_minus_data(particle_3->coord, particle_1->coord), my_minus_data(particle_4->coord, particle_1->coord)))/2.0; 
    Real s3=get_distance(get_cross_product(my_minus_data(particle_2->coord, particle_1->coord), my_minus_data(particle_4->coord, particle_1->coord)))/2.0; 
    Real s4=get_distance(get_cross_product(my_minus_data(particle_3->coord, particle_1->coord), my_minus_data(particle_2->coord, particle_1->coord)))/2.0; 
    particle_temp1.coord.i=(s1*particle_1->coord.i+s2*particle_2->coord.i+s3*particle_3->coord.i+s4*particle_4->coord.i)/(s1+s2+s3+s4); 
    particle_temp1.coord.j=(s1*particle_1->coord.j+s2*particle_2->coord.j+s3*particle_3->coord.j+s4*particle_4->coord.j)/(s1+s2+s3+s4); 
    particle_temp1.coord.k=(s1*particle_1->coord.k+s2*particle_2->coord.k+s3*particle_3->coord.k+s4*particle_4->coord.k)/(s1+s2+s3+s4);
    #if (DIM_X)
    particle_temp1.coord.i = AMIN1(AMAX1(particle_temp1.coord.i, box_l.i+4.0*dl), box_r.i-4.0*dl);
    #endif
    #if (DIM_Y)
    particle_temp1.coord.j = AMIN1(AMAX1(particle_temp1.coord.j, box_l.j+4.0*dl), box_r.j-4.0*dl);
    #endif
    #if (DIM_Z)
    particle_temp1.coord.k = AMIN1(AMAX1(particle_temp1.coord.k, box_l.k+4.0*dl), box_r.k-4.0*dl);
    #endif

    particle_temp1.Calculate_particle_infor(sph);

    tetgenmesh m;

    tetgenmesh::triface tetloop, neightet;
    Real *p[4];
    Real ccent[3] = {0., 0., 0.};
    Real pp[12];

    pp[ 0] = particle_1->coord.i;
    pp[ 1] = particle_1->coord.j;
    pp[ 2] = particle_1->coord.k;

    pp[ 3] = particle_2->coord.i;
    pp[ 4] = particle_2->coord.j;
    pp[ 5] = particle_2->coord.k;

    pp[ 6] = particle_3->coord.i;
    pp[ 7] = particle_3->coord.j;
    pp[ 8] = particle_3->coord.k;

    pp[ 9] = particle_4->coord.i;
    pp[10] = particle_4->coord.j;
    pp[11] = particle_4->coord.k;

    // Get four vertices: p0, p1, p2, p3.
    for (int i = 0; i < 4; i++) {
      p[i] = &pp[i*3];
    }
    // Get the tet volume.
    m.circumsphere(p[0], p[1], p[2], p[3], ccent, NULL);

    Particle particle_temp2;
    particle_temp2.coord.i = ccent[0];
    particle_temp2.coord.j = ccent[1];
    particle_temp2.coord.k = ccent[2];
    #if (DIM_X)
    particle_temp2.coord.i = AMIN1(AMAX1(particle_temp2.coord.i, box_l.i+4.0*dl), box_r.i-4.0*dl);
    #endif
    #if (DIM_Y)
    particle_temp2.coord.j = AMIN1(AMAX1(particle_temp2.coord.j, box_l.j+4.0*dl), box_r.j-4.0*dl);
    #endif
    #if (DIM_Z)
    particle_temp2.coord.k = AMIN1(AMAX1(particle_temp2.coord.k, box_l.k+4.0*dl), box_r.k-4.0*dl);
    #endif

    particle_temp2.Calculate_particle_infor(sph);

    Real phi_size = 1./4.*(particle_1->phi + particle_2->phi + particle_3->phi + particle_4->phi);
    Real h_avg    = 1./4.*(particle_1->h   + particle_2->h   + particle_3->h   + particle_4->h);
    Real curv_max = AMAX1(particle_4->curv, AMAX1(particle_3->curv, AMAX1(particle_1->curv, particle_2->curv)));
    Real curv_min = AMIN1(particle_4->curv, AMIN1(particle_3->curv, AMIN1(particle_1->curv, particle_2->curv)));
    Real Rmin = AMIN1(0.5*h_avg, 1./(curv_max+1.e-10));
    Real Rmax = AMIN1(1.0*h_avg, 1./(curv_min+1.e-10));
    // if(particle_temp1.phi <= 0.0001*phi_size){
    // if(particle_temp1.phi <= -0.05*h_avg || particle_temp2.phi <= -1.25*h_avg){
    if(particle_temp1.phi <= -Rmin || particle_temp2.phi <= -Rmax){
    // if(particle_temp1.phi <= -0.5*dl || particle_temp2.phi <= -dl){
      // return true;
    }

    // Particle particle_temp1;
    // particle_temp1.coord.i = 1./3.*(particle_1->coord.i + particle_2->coord.i + particle_3->coord.i);
    // particle_temp1.coord.j = 1./3.*(particle_1->coord.j + particle_2->coord.j + particle_3->coord.j);
    // particle_temp1.coord.k = 1./3.*(particle_1->coord.k + particle_2->coord.k + particle_3->coord.k);
    // particle_temp1.Calculate_particle_infor(sph);

    // Real phi_size = 1./3.*(particle_1->phi + particle_2->phi + particle_3->phi);
    // Real h_avg    = 1./3.*(particle_1->h   + particle_2->h   + particle_3->h);

    // // if(particle_temp1.phi <= -0.1*phi_size){
    // if(particle_temp1.phi <= -0.25*h_avg){
    //   return true;
    // }

    // Particle particle_temp2;
    // particle_temp2.coord.i = 1./3.*(particle_2->coord.i + particle_3->coord.i + particle_4->coord.i);
    // particle_temp2.coord.j = 1./3.*(particle_2->coord.j + particle_3->coord.j + particle_4->coord.j);
    // particle_temp2.coord.k = 1./3.*(particle_2->coord.k + particle_3->coord.k + particle_4->coord.k);
    // particle_temp2.Calculate_particle_infor(sph);

    // phi_size = 1./3.*(particle_2->phi + particle_3->phi + particle_4->phi);
    // h_avg    = 1./3.*(particle_2->h   + particle_3->h   + particle_4->h);

    // // if(particle_temp2.phi <= -0.1*phi_size){
    // if(particle_temp2.phi <= -0.25*h_avg){
    //   return true;
    // }

    // Particle particle_temp3;
    // particle_temp3.coord.i = 1./3.*(particle_1->coord.i + particle_2->coord.i + particle_4->coord.i);
    // particle_temp3.coord.j = 1./3.*(particle_1->coord.j + particle_2->coord.j + particle_4->coord.j);
    // particle_temp3.coord.k = 1./3.*(particle_1->coord.k + particle_2->coord.k + particle_4->coord.k);
    // particle_temp3.Calculate_particle_infor(sph);

    // phi_size = 1./3.*(particle_1->phi + particle_2->phi + particle_4->phi);
    // h_avg    = 1./3.*(particle_1->h   + particle_2->h   + particle_4->h);

    // // if(particle_temp3.phi <= -0.1*phi_size){
    // if(particle_temp3.phi <= -0.25*h_avg){
    //   return true;
    // }

    // Particle particle_temp4;
    // particle_temp4.coord.i = 1./3.*(particle_1->coord.i + particle_3->coord.i + particle_4->coord.i);
    // particle_temp4.coord.j = 1./3.*(particle_1->coord.j + particle_3->coord.j + particle_4->coord.j);
    // particle_temp4.coord.k = 1./3.*(particle_1->coord.k + particle_3->coord.k + particle_4->coord.k);
    // particle_temp4.Calculate_particle_infor(sph);

    // phi_size = 1./3.*(particle_1->phi + particle_3->phi + particle_4->phi);
    // h_avg    = 1./3.*(particle_1->h   + particle_3->h   + particle_4->h);

    // // if(particle_temp4.phi <= -0.1*phi_size){
    // if(particle_temp4.phi <= -0.25*h_avg){
    //   return true;
    // }

    if( (particle_1->type != REAL_PARTICLE && 
         particle_2->type != REAL_PARTICLE && 
         particle_3->type != REAL_PARTICLE &&
         particle_4->type != REAL_PARTICLE    )||
        (count_color > 0 && stage == 2))
    {
      tetgenmesh m;

      int  i,j;
      int  indx[4];
      Real pp[12];
      Real edgelength[6];
      Real *p[4];
      Real V[6][3], N[4][3], H[4]; // edge-vectors, face-normals, face-heights.
      Real A[4][4], rhs[4], D;
      Real longlen = 0.;
      Real shortlen = 0.;
      Real insradius = 0.;
      Real minheightinv = 0.;
      Real cirradius = 0.;

      pp   [0*3 + 0] = particle_1->coord.i;
      pp   [0*3 + 1] = particle_1->coord.j;
      pp   [0*3 + 2] = particle_1->coord.k;
      pp   [1*3 + 0] = particle_2->coord.i;
      pp   [1*3 + 1] = particle_2->coord.j;
      pp   [1*3 + 2] = particle_2->coord.k;
      pp   [2*3 + 0] = particle_3->coord.i;
      pp   [2*3 + 1] = particle_3->coord.j;
      pp   [2*3 + 2] = particle_3->coord.k;
      pp   [3*3 + 0] = particle_4->coord.i;
      pp   [3*3 + 1] = particle_4->coord.j;
      pp   [3*3 + 2] = particle_4->coord.k;

      // Get four vertices: p0, p1, p2, p3.
      for (i = 0; i < 4; i++) {
        p[i] = &pp[i*3];
      }

      // Get the tet volume.
      Real tetvol = m.orient3dfast(p[1], p[0], p[2], p[3]) / 6.0;

      // Set the edge vectors: V[0], ..., V[5]
      for (i = 0; i < 3; i++) V[0][i] = p[0][i] - p[3][i]; // V[0]: p3->p0.
      for (i = 0; i < 3; i++) V[1][i] = p[1][i] - p[3][i]; // V[1]: p3->p1.
      for (i = 0; i < 3; i++) V[2][i] = p[2][i] - p[3][i]; // V[2]: p3->p2.
      for (i = 0; i < 3; i++) V[3][i] = p[1][i] - p[0][i]; // V[3]: p0->p1.
      for (i = 0; i < 3; i++) V[4][i] = p[2][i] - p[1][i]; // V[4]: p1->p2.
      for (i = 0; i < 3; i++) V[5][i] = p[0][i] - p[2][i]; // V[5]: p2->p0.

      // Get the squares of the edge lengths.
      for (i = 0; i < 6; i++) edgelength[i] = m.dot(V[i], V[i]);

      // Calculate the longest and shortest edge length.
      for (i = 0; i < 6; i++) {
        if (i == 0) {
          shortlen = longlen = edgelength[i];
        } else {
          shortlen = edgelength[i] < shortlen ? edgelength[i] : shortlen;
          longlen  = edgelength[i] > longlen  ? edgelength[i] : longlen;
        }
      }

      // Set the matrix A = [V[0], V[1], V[2]]^T.
      for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++) A[j][i] = V[j][i];
      }

      // Decompose A just once.
      if (m.lu_decmp(A, 3, indx, &D, 0)) {   
        // Get the three faces normals.
        for (j = 0; j < 3; j++) {
          for (i = 0; i < 3; i++) rhs[i] = 0.0;
          rhs[j] = 1.0;  // Positive means the inside direction
          m.lu_solve(A, 3, indx, rhs, 0);
          for (i = 0; i < 3; i++) N[j][i] = rhs[i];
        }
        // Get the fourth face normal by summing up the first three.
        for (i = 0; i < 3; i++) N[3][i] = - N[0][i] - N[1][i] - N[2][i];
        // Get the radius of the circumsphere.
        for (i = 0; i < 3; i++) rhs[i] = 0.5 * m.dot(V[i], V[i]);
        m.lu_solve(A, 3, indx, rhs, 0);
        cirradius = sqrt(m.dot(rhs, rhs));
        // Normalize the face normals.
        for (i = 0; i < 4; i++) {
          // H[i] is the inverse of height of its corresponding face.
          H[i] = sqrt(m.dot(N[i], N[i]));
          for (j = 0; j < 3; j++) N[i][j] /= H[i];
        }
        // Get the radius of the inscribed sphere.
        insradius = 1.0 / (H[0] + H[1] + H[2] + H[3]);
        // Get the biggest H[i] (corresponding to the smallest height).
        minheightinv = H[0];
        for (i = 1; i < 4; i++) {
          if (H[i] > minheightinv) minheightinv = H[i];
        }
      } else {
        // A nearly degenerated tet.
        // if (tetvol < 0.0) {
        //   cout<<"<<<<< !! Warning:  A inverted tet ["<<p_[0]<<" "<<p_[1]<<" "<<p_[2]<<" "<<p_[3]<<" ]"<<endl;
        //   continue;
        // }else if (tetvol == 0.0){
        //   cout<<"<<<<< !! Warning:  A degenerated tet ["<<p_[0]<<" "<<p_[1]<<" "<<p_[2]<<" "<<p_[3]<<" ]"<<endl;
        // }
        // Calculate the four face normals.
        m.facenormal(p[2], p[1], p[3], N[0], 1, NULL);
        m.facenormal(p[0], p[2], p[3], N[1], 1, NULL);
        m.facenormal(p[1], p[0], p[3], N[2], 1, NULL);
        m.facenormal(p[0], p[1], p[2], N[3], 1, NULL);
        // Normalize the face normals.
        for (i = 0; i < 4; i++) {
          // H[i] is the twice of the area of the face.
          H[i] = sqrt(m.dot(N[i], N[i]));
          for (j = 0; j < 3; j++) N[i][j] /= H[i];
        }
        // Get the biggest H[i] / tetvol (corresponding to the smallest height).
        minheightinv = (H[0] / tetvol);
        for (i = 1; i < 4; i++) {
          if ((H[i] / tetvol) > minheightinv) minheightinv = (H[i] / tetvol);
        }
        // Let the circumradius to be the half of its longest edge length.
        cirradius = 0.5 * sqrt(longlen);
      }

      Real tetaspect = sqrt(longlen) * minheightinv;

      if (tetaspect > sph->tet_delete_thres) 
        return true;
    }

  return false;
}
//-------------------------------------------------------
// Reconstruct::find tetrahedron V0.0.3
//-------------------------------------------------------
void Graphmeshcls::Reconstruct_tets(std::vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< 3D tetrahedralization processing"<<endl;
  }

  #ifdef _GEN_DUMMY_PART_
  int num_dummy = sph->dummy_particle.size();
  #endif
  
  tetgenio pp_in, tet_out;

  // All indices start from 0
  pp_in.firstnumber = 0;
  tet_out.firstnumber = 0;

  int num_total = particle_total.size();
  pp_in.numberofpoints = num_total;

  #ifdef _GEN_DUMMY_PART_
  pp_in.numberofpoints += num_dummy;
  #endif

  pp_in.pointlist = new Real[pp_in.numberofpoints * 3];
  for (int i=0; i < num_total; i++){
    pp_in.pointlist[i * 3 + 0] = particle_total[i]->coord.i;
    pp_in.pointlist[i * 3 + 1] = particle_total[i]->coord.j;
    pp_in.pointlist[i * 3 + 2] = particle_total[i]->coord.k;
  }

  #ifdef _GEN_DUMMY_PART_
  for (int i=num_total; i < num_total+num_dummy; i++){
    pp_in.pointlist[i * 3 + 0] = sph->dummy_particle[i-num_total]->coord.i;
    pp_in.pointlist[i * 3 + 1] = sph->dummy_particle[i-num_total]->coord.j;
    pp_in.pointlist[i * 3 + 2] = sph->dummy_particle[i-num_total]->coord.k;
  }
  #endif

  // do the tetrahedralization 1st time (without face swapping and flipping)
  tetrahedralize("JCnQ", &pp_in, &tet_out);
  // tetrahedralize("JCnlVV", &pp_in, &tet_out);
  // tetrahedralize("CnlQ", &pp_in, &tet_out);
  // tetrahedralize("CnlVVO/2", &pp_in, &tet_out);

  tris.clear();
  tets.clear();

  //load the tets from tetgen the first time
  int   n_tets       = tet_out.numberoftetrahedra;
  int  *tet_list     = tet_out.tetrahedronlist;
  int  *neighborlist = tet_out.neighborlist;
  if (neighborlist == (int *) NULL) cout<<"<<<<< ["<<world.rank()<<"] !! ERROR neighborlist empty"<<endl;

  for (int i=0; i < n_tets; i++){
    std::vector<int> vert_list; vert_list.clear();
    for (int j = 0; j < 4; ++j)
    {
      vert_list.push_back(tet_list[i*4 + j]);
    }
    // std::sort(vert_list.begin(), vert_list.end(), Greater());

    int4_graph tet;
    tet.i      = vert_list[0];
    tet.j      = vert_list[1];
    tet.k      = vert_list[2];
    tet.l      = vert_list[3];
    tet.vol    = 0.;
    tet.aspect = 0.;
    tets.insert (tet);
  }

  if (tet_out.numberofpoints != pp_in.numberofpoints )
    cout<<"<<<<< ["<<world.rank()<<"] ERROR!! number of output points is wrong in 1st tetrahedralize: tet_out.numberofpoints "<<tet_out.numberofpoints<<" pp_in.numberofpoints "<<pp_in.numberofpoints<<endl;

  for (int i=0; i < num_total; i++){
    my_real pp;
    pp.i = tet_out.pointlist[i * 3 + 0];
    pp.j = tet_out.pointlist[i * 3 + 1];
    pp.k = tet_out.pointlist[i * 3 + 2];
    
    Real dist = get_distance_2p(particle_total[i]->coord, pp);
    if (dist> 1.e-12){
      cout<<"<<<<< ["<<world.rank()<<"] ERROR!! output coordinates from teten (the 1st time) not the same, error: "<<dist<<endl;
    }
//     particle_total[i]->coord.i = tet_out.pointlist[i * 3 + 0];
//     particle_total[i]->coord.j = tet_out.pointlist[i * 3 + 1];
//     particle_total[i]->coord.k = tet_out.pointlist[i * 3 + 2];

//     int tag_interface      = NORMAL_CELL;
//     int tag_characteristic = NORMAL_CELL;
//     int idx_characteristic = -1;
//     particle_total[i]->Calculate_particle_infor(sph);
//     particle_total[i]->Get_extended_char_cell_tag_and_idx (sph, tag_interface, tag_characteristic, idx_characteristic);
    // particle_total[i]->Get_current_particle_level_set_cell_tag (sph, tag_interface, tag_characteristic);

//     if ( tag_interface != NORMAL_CELL ){
//       if (tag_characteristic == SINGULARITY_CELL)
//         particle_total[i]->type = SINGULARITY_PARTICLE;
//       else if (tag_characteristic == SEGMENT_CELL)
//         particle_total[i]->type = SEGMENT_PARTICLE;
//       else if (tag_characteristic == NORMAL_CELL && tag_interface == CUT_CELL)
//         particle_total[i]->type = SURFACE_PARTICLE;
//       else
//         particle_total[i]->type = REAL_PARTICLE;
//     }else{
//       particle_total[i]->type = REAL_PARTICLE;
//     }
  }
  
  int  iRank = world.rank();
  std::set<int4_graph>::iterator it;
  std::vector<std::set<int4_graph>::iterator> tet_delete;
  tet_delete.clear();

  for (it=tets.begin(); it!=tets.end(); ++it){
    int p1 = int((*it).i);
    int p2 = int((*it).j);
    int p3 = int((*it).k);
    int p4 = int((*it).l);

    #ifdef _GEN_DUMMY_PART_
    p_Particle particle_1; 
    p_Particle particle_2;
    p_Particle particle_3;
    p_Particle particle_4;

    if (p1 < num_total) particle_1 = particle_total[p1];
    else                particle_1 = sph->dummy_particle[p1-num_total];
    if (p2 < num_total) particle_2 = particle_total[p2];
    else                particle_2 = sph->dummy_particle[p2-num_total];
    if (p3 < num_total) particle_3 = particle_total[p3];
    else                particle_3 = sph->dummy_particle[p3-num_total];
    if (p4 < num_total) particle_4 = particle_total[p4];
    else                particle_4 = sph->dummy_particle[p4-num_total];

    #else
    p_Particle particle_1 = particle_total[p1];
    p_Particle particle_2 = particle_total[p2];
    p_Particle particle_3 = particle_total[p3];
    p_Particle particle_4 = particle_total[p4];
    #endif

    bool error_tet = false;

    error_tet = Check_whether_this_tet_should_be_delete(1, particle_1, particle_2, particle_3, particle_4, iRank, sph);

    if (error_tet) tet_delete.push_back(it);

    // find boundary surfaces if there is any
    if (!error_tet){
      // if (sph->level_set.num_singularity == 0){
        int type[4];
        type[0] = particle_1->type;
        type[1] = particle_2->type;
        type[2] = particle_3->type;
        type[3] = particle_4->type;

        int sum = 0;
        int count_real = 0;
        int count_surf = 0;
        int count_sing = 0;
        int count_segm = 0;

        for (int i = 0; i < 4; ++i)
        {
          if (type[i] == REAL_PARTICLE)           count_real++;
          if (type[i] == SURFACE_PARTICLE)        count_surf++;
          if (type[i] == SINGULARITY_PARTICLE)    count_sing++;
          if (type[i] == SEGMENT_PARTICLE)        count_segm++;
        }

        if (count_real == 1){
          // find a surface tri in the tet
          int3_graph tri;
          if       (type[0] == REAL_PARTICLE){
            tri.i = p2;
            tri.j = p3;
            tri.k = p4;
            tris.insert(tri);
          }else if (type[1] == REAL_PARTICLE){
            tri.i = p1;
            tri.j = p3;
            tri.k = p4;
            tris.insert(tri);
          }else if (type[2] == REAL_PARTICLE){
            tri.i = p1;
            tri.j = p2;
            tri.k = p4;
            tris.insert(tri);
          }else if (type[3] == REAL_PARTICLE){
            tri.i = p1;
            tri.j = p2;
            tri.k = p3;
            tris.insert(tri);
          }
        }else if (count_real == 0 && count_sing == 1){
          int3_graph tri;
          if       (type[0] == SINGULARITY_PARTICLE){
            tri.i = p1;
            tri.j = p2;
            tri.k = p3;
            tris.insert(tri);
            tri.k = p4;
            tris.insert(tri);
            tri.j = p3;
            tri.k = p4;
            tris.insert(tri);
          }else if (type[1] == SINGULARITY_PARTICLE){
            tri.i = p1;
            tri.j = p2;
            tri.k = p3;
            tris.insert(tri);
            tri.k = p4;
            tris.insert(tri);
            tri.i = p2;
            tri.j = p3;
            tris.insert(tri);
          }else if (type[2] == SINGULARITY_PARTICLE){
            tri.i = p1;
            tri.j = p2;
            tri.k = p3;
            tris.insert(tri);
            tri.i = p2;
            tri.j = p3;
            tri.k = p4;
            tris.insert(tri);
            tri.i = p1;
            tris.insert(tri);
          }else if (type[3] == SINGULARITY_PARTICLE){
            tri.i = p1;
            tri.j = p2;
            tri.k = p4;
            tris.insert(tri);
            tri.j = p3;
            tris.insert(tri);
            tri.i = p2;
            tris.insert(tri);
          }
        }else if (count_real == 0 && count_segm == 2){
          int pp[4]; pp[0] = p1; pp[1] = p2; pp[2] = p3; pp[3] = p4;
          for (int i = 0; i < 4; ++i){
            if (type[i] == SEGMENT_PARTICLE){
              for (int j = i; j < 4; ++j){
                if (type[j] == SEGMENT_PARTICLE){
                  for (int k = 0; k < 4; ++k){
                    if (type[k] != SEGMENT_PARTICLE){
                      int3_graph tri;
                      tri.i=pp[i];
                      tri.j=pp[j];
                      tri.k=pp[k];
                    }
                  }
                }            
              }
            }
          }
        }
      // }else{
      //   cout<<"<<<<< ERROR singularities are not able to be considered here!!!"<<endl;
      // }
    }
  }

  if (tet_delete.size()>0){
    int n_delete = tet_delete.size();
    cout<<"<<<<< ["<<world.rank()<<"] Num. of tets to be deleted: "<<tet_delete.size()<<endl;

    for (int i=0; i < n_delete; i++){
      tets.erase(tet_delete[i]);
    }

    if (tets.size() + n_delete != tet_out.numberoftetrahedra){
      cout<<"<<<<< ["<<world.rank()<<"] ERROR: number of tets are wrong!!!"<<endl;
    }
  }

  #if !defined (_READ_TARGET_FIELD_FOR_POST_)
  if (sph->run_time >= sph->v_reini_change){
  #endif
    
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Flipping to improve mesh quality"<<endl;
  }
  
  #ifdef _MPI_
  // JZ20190101::do the flipping currently using a master-slave strategy
  // it is not very easy to perform parallel one
  // TODO::looking for existing parallel code to do the job
  int numberofpoints = sph->total_num_particle;
  int totalnumberofpoints = sph->glbl_total_num_particle;

  int numberoftetrahedra = tets.size();
  int nproc = world.size();
    
  int *local_id_table;
  local_id_table = new int [totalnumberofpoints];
  
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, totalnumberofpoints),
            [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      local_id_table[i] = -1;
    }
  }, ap);
  
  int numberoflocalidpoints = particle_total.size();
  parallel_for( blocked_range<int>(0, numberoflocalidpoints),
            [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      local_id_table[particle_total[i]->id] = particle_total[i]->local_id;
    }
  }, ap);
  
  serialization_vector <Real> out_coord;
  out_coord.Vector.clear();
  serialization_vector <int> out_tets;
  out_tets.Vector.clear();
  
  if (world.rank() == 0){
    serialization_vector <Real> *gather_coord;
    gather_coord = new serialization_vector <Real>[nproc];
    serialization_vector <int> *gather_tets;
    gather_tets = new serialization_vector <int>[nproc];
    
    parallel_for( blocked_range<int>(0, nproc),
              [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        gather_coord[i].Vector.clear();
        gather_coord[i].mem_size = 0;
        gather_coord[i].tag = 0;
        gather_tets [i].Vector.clear();
        gather_tets [i].mem_size = 0;
        gather_tets [i].tag = 0;
      }
    }, ap);
    
    cout<<"[0]-> gathering data"<<endl;
    
    gather(world,out_coord,gather_coord,0);
    gather(world,out_tets ,gather_tets ,0);
    
    cout<<"[0]-> data received"<<endl;
    
    tetgenio::facet *f;
    tetgenio::polygon *p;

    pp_in.deinitialize();
    pp_in.initialize();

    tet_out.deinitialize();
    tet_out.initialize();

    pp_in.firstnumber = 0;
    tet_out.firstnumber = 0;
    
    int off_set[nproc];
    off_set[0] = 0;
    
    for (int i = 1; i < nproc; i++){
      off_set[i] = numberofpoints;
      numberofpoints += int(gather_coord[i].Vector.size()/3);
      numberoftetrahedra += int(gather_tets[i].Vector.size()/4);
      cout<<"[0]->number of points received from proc "<<i<<" is "<<int(gather_coord[i].Vector.size()/3)<<endl;
    }
    
    if (numberofpoints != totalnumberofpoints){
      cout<<"<<<<< ERROR!!! numberofpoints wrong!!! numberofpoints: "<<numberofpoints<<" glbl_total_num_particle: "<<totalnumberofpoints<<endl;
    }
    
    pp_in.numberofpoints = numberofpoints;
    pp_in.pointlist = new Real[pp_in.numberofpoints * 3];
    int idd = 0;
    int numberofpoints_in_master = sph->total_num_particle;
    for (int i=0; i < numberofpoints_in_master; i++){
      pp_in.pointlist[i * 3 + 0] = particle_total[i]->coord.i;
      pp_in.pointlist[i * 3 + 1] = particle_total[i]->coord.j;
      pp_in.pointlist[i * 3 + 2] = particle_total[i]->coord.k;
      idd++;
    }
    for (int i = 1; i < nproc; i++){
      int nj = int(gather_coord[i].Vector.size()/3);
      for(int j = 0; j < nj; j++){
        pp_in.pointlist[idd * 3 + 0] = gather_coord[i].Vector[j * 3 + 0];
        pp_in.pointlist[idd * 3 + 1] = gather_coord[i].Vector[j * 3 + 1];
        pp_in.pointlist[idd * 3 + 2] = gather_coord[i].Vector[j * 3 + 2];
        idd++;
      }
    }
    if (idd != totalnumberofpoints){
      cout<<"<<<<< ERROR!!! idd wrong!!! idd: "<<idd<<" glbl_total_num_particle: "<<totalnumberofpoints<<endl;
    }
    
    pp_in.numberoftetrahedra = numberoftetrahedra;
    pp_in.tetrahedronlist = new int[pp_in.numberoftetrahedra*4];
    int itet = 0;
    for (it=tets.begin(); it!=tets.end(); ++it){
      pp_in.tetrahedronlist[itet*4 + 0] = particle_total[int((*it).i)]->id;
      pp_in.tetrahedronlist[itet*4 + 1] = particle_total[int((*it).j)]->id;
      pp_in.tetrahedronlist[itet*4 + 2] = particle_total[int((*it).k)]->id;
      pp_in.tetrahedronlist[itet*4 + 3] = particle_total[int((*it).l)]->id;
      itet ++;
    }
    for (int i = 1; i < nproc; i++){
      int nj = int(gather_tets[i].Vector.size()/4);
      for(int j = 0; j < nj; j++){
        pp_in.tetrahedronlist[itet*4 + 0] = gather_tets[i].Vector[j * 4 + 0];
        pp_in.tetrahedronlist[itet*4 + 1] = gather_tets[i].Vector[j * 4 + 1];
        pp_in.tetrahedronlist[itet*4 + 2] = gather_tets[i].Vector[j * 4 + 2];
        pp_in.tetrahedronlist[itet*4 + 3] = gather_tets[i].Vector[j * 4 + 3];
        itet ++;
      }
    }
    if (itet != numberoftetrahedra){
      cout<<"<<<<< ERROR!!! itet wrong!!! itet: "<<itet<<" numberoftetrahedra: "<<numberoftetrahedra<<endl;
    }
    
    cout<<"[0]-> tetgen performing flipping"<<endl;
    
    tetrahedralize("JrO/1o/100Q", &pp_in, &tet_out);
    
    cout<<"[0]-> flipping finished"<<endl;
    
    if (tet_out.numberofpoints != totalnumberofpoints)
      cout<<"<<<<< ERROR!! number of output points is wrong in 2nd tetrahedralize: glbl_total_num_particles: "<<totalnumberofpoints<<" tet_out.numberofpoints: "<<tet_out.numberofpoints<<endl;

    cout<<"[0]-> preparing for sending data"<<endl;

    parallel_for( blocked_range<int>(0, nproc),
              [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        gather_tets [i].Vector.clear();
        gather_tets [i].mem_size = 0;
        gather_tets [i].tag = 0;
      }
    }, ap);
    
    n_tets       = tet_out.numberoftetrahedra;
    tet_list     = tet_out.tetrahedronlist;
    tets.clear();
    
    std::vector<int> vert_list;
    for (int i=0; i < n_tets; i++){
      vert_list.clear();
      int color_max = 0;
      for (int j = 0; j < 4; ++j)
      {
        vert_list.push_back(tet_list[i*4 + j]);
        int color_ = Get_color (vert_list[j], off_set, nproc, totalnumberofpoints);
        if (color_ >= 0)
          color_max = AMAX1(color_max, color_);
        else{
          cout<<"[0]-> something goes wrong with the color returned ("<<color_<<")"<<endl;
          world.abort(-1);
        }
      }
      
      if (color_max == 0){
        // JZ20190102::master node only contains local points
        int4_graph tet;
        tet.i      = local_id_table[vert_list[0]];
        tet.j      = local_id_table[vert_list[1]];
        tet.k      = local_id_table[vert_list[2]];
        tet.l      = local_id_table[vert_list[3]];
        tet.vol    = 0.;
        tet.aspect = 0.;
        tets.insert (tet);
      }else if (color_max > 0 && color_max < nproc){
        gather_tets[color_max].Vector.push_back(vert_list[0]);
        gather_tets[color_max].Vector.push_back(vert_list[1]);
        gather_tets[color_max].Vector.push_back(vert_list[2]);
        gather_tets[color_max].Vector.push_back(vert_list[3]);
      }else{
        cout<<"[0]-> something goes wrong with the color_max"<<endl;
      }
    }
    cout<<"[0]-> sending data"<<endl;
    
    mpi::request reqs[nproc];
    mpi::request reqs_2nd_recv[nproc];
    mpi::request reqs_2nd_send_color[nproc];
    mpi::request reqs_2nd_send_coord[nproc];
    
    for(int i=1; i!=nproc; ++i){
      gather_tets[i].mem_size = int(gather_tets[i].Vector.size());
      gather_tets[i].tag = 1;
      //nonblocking communication
      reqs[i] = world.isend(i, i, gather_tets[i]);
    }

    serialization_vector <int> *request_id;
    request_id = new serialization_vector <int>[nproc];
    serialization_vector <int> *color_to_send;
    color_to_send = new serialization_vector <int>[nproc];
    serialization_vector <Real> *coord_to_send;
    coord_to_send = new serialization_vector <Real>[nproc];
    
    parallel_for( blocked_range<int>(0, nproc),
              [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        request_id[i].Vector.clear();
        request_id[i].mem_size = 0;
        request_id[i].tag = 0;
        color_to_send [i].Vector.clear();
        color_to_send [i].mem_size = 0;
        color_to_send [i].tag = 0;
        coord_to_send [i].Vector.clear();
        coord_to_send [i].mem_size = 0;
        coord_to_send [i].tag = 0;
      }
    }, ap);
    
    for(int i=1; i!=nproc; ++i){
      reqs_2nd_recv[i] = world.irecv(i, i, request_id[i]);
    }
    mpi::wait_all(reqs_2nd_recv+1, reqs_2nd_recv+nproc);
    
    for(int i=1; i!=nproc; ++i){
      if(request_id[i].Vector.size() == 1 && request_id[i].Vector[0] == -1) continue;
      else if (request_id[i].Vector.size() >= 1 && request_id[i].Vector[0] >= 0){
        int num_request = request_id[i].Vector.size();
        for (int j = 0; j < num_request; j++){
          int _id = request_id[i].Vector[j];
          coord_to_send [i].Vector.push_back(tet_out.pointlist[_id * 3 + 0]);
          coord_to_send [i].Vector.push_back(tet_out.pointlist[_id * 3 + 1]);
          coord_to_send [i].Vector.push_back(tet_out.pointlist[_id * 3 + 2]);
          
          color_to_send [i].Vector.push_back(Get_color(_id, off_set, nproc, totalnumberofpoints));
        }
        
        coord_to_send[i].mem_size = int(coord_to_send[i].Vector.size());
        coord_to_send[i].tag = 1;
        color_to_send[i].mem_size = int(color_to_send[i].Vector.size());
        color_to_send[i].tag = 1;
        //nonblocking communication
        cout<<"[0]-> sending "<<num_request<<" points to proc "<<i<<endl;
        reqs_2nd_send_color[i] = world.isend(i, 2*i, color_to_send[i]);
        reqs_2nd_send_coord[i] = world.isend(i, 3*i, coord_to_send[i]);
      
      }else {
        cout<<"[0]-> ERROR!!! something wrong with the request extra points from proc "<<i<<endl;
        world.abort(-1);
      }
    }
    
  #else
    // do the tetrahedralization 2nd time (with face swapping and flipping)
    // All indices start from 0
    tetgenio::facet *f;
    tetgenio::polygon *p;

    pp_in.deinitialize();
    pp_in.initialize();

    tet_out.deinitialize();
    tet_out.initialize();

    pp_in.firstnumber = 0;
    tet_out.firstnumber = 0;

    pp_in.numberofpoints = particle_total.size();
    pp_in.pointlist = new Real[pp_in.numberofpoints * 3];
    for (int i=0; i < pp_in.numberofpoints; i++){
      pp_in.pointlist[i * 3 + 0] = particle_total[i]->coord.i;
      pp_in.pointlist[i * 3 + 1] = particle_total[i]->coord.j;
      pp_in.pointlist[i * 3 + 2] = particle_total[i]->coord.k;
    }

    pp_in.numberoftetrahedra = tets.size();
    pp_in.tetrahedronlist = new int[pp_in.numberoftetrahedra*4];
    int itet = 0;
    for (it=tets.begin(); it!=tets.end(); ++it){
      pp_in.tetrahedronlist[itet*4 + 0] = int((*it).i);
      pp_in.tetrahedronlist[itet*4 + 1] = int((*it).j);
      pp_in.tetrahedronlist[itet*4 + 2] = int((*it).k);
      pp_in.tetrahedronlist[itet*4 + 3] = int((*it).l);
      itet ++;
    }
    
    tetrahedralize("JrO/1o/100Q", &pp_in, &tet_out);
  //   tetrahedralize("JrVVO/1o/100Q", &pp_in, &tet_out);
    // tetrahedralize("pdYClVVO/2", &pp_in, &tet_out);

    tets.clear();

    //load the tets from tetgen the first time
    n_tets       = tet_out.numberoftetrahedra;
    tet_list     = tet_out.tetrahedronlist;
    neighborlist = tet_out.neighborlist;

    for (int i=0; i < n_tets; i++){
      std::vector<int> vert_list; vert_list.clear();
      for (int j = 0; j < 4; ++j)
      {
        vert_list.push_back(tet_list[i*4 + j]);
      }
      // std::sort(vert_list.begin(), vert_list.end(), Greater());

      int4_graph tet;
      tet.i      = vert_list[0];
      tet.j      = vert_list[1];
      tet.k      = vert_list[2];
      tet.l      = vert_list[3];
      tet.vol    = 0.;
      tet.aspect = 0.;
      tets.insert (tet);
    }

    if (tet_out.numberofpoints != particle_total.size())
      cout<<"<<<<< ["<<world.rank()<<"] ERROR!! number of output points is wrong in 2nd tetrahedralize: tet_out.numberofpoints "<<tet_out.numberofpoints<<" particle_total.size "<<particle_total.size()<<endl;

    for (int i=0; i < tet_out.numberofpoints; i++){
      my_real pp;
      pp.i = tet_out.pointlist[i * 3 + 0];
      pp.j = tet_out.pointlist[i * 3 + 1];
      pp.k = tet_out.pointlist[i * 3 + 2];
      
      Real dist = get_distance_2p(particle_total[i]->coord, pp);
      if (dist> 1.e-12){
        cout<<"<<<<< ["<<world.rank()<<"] ERROR!! output coordinates from teten (the 2nd time) not the same, error: "<<dist<<endl;
      }
      
  //     particle_total[i]->coord.i = tet_out.pointlist[i * 3 + 0];
  //     particle_total[i]->coord.j = tet_out.pointlist[i * 3 + 1];
  //     particle_total[i]->coord.k = tet_out.pointlist[i * 3 + 2];
  // 
  //     int tag_interface      = NORMAL_CELL;
  //     int tag_characteristic = NORMAL_CELL;
  //     int idx_characteristic = -1;
  //     particle_total[i]->Calculate_particle_infor(sph);
  //     particle_total[i]->Get_extended_char_cell_tag_and_idx (sph, tag_interface, tag_characteristic, idx_characteristic);
  //     // particle_total[i]->Get_current_particle_level_set_cell_tag (sph, tag_interface, tag_characteristic);
  // 
  //     if ( tag_interface != NORMAL_CELL ){
  //       if (tag_characteristic == SINGULARITY_CELL)
  //         particle_total[i]->type = SINGULARITY_PARTICLE;
  //       else if (tag_characteristic == SEGMENT_CELL)
  //         particle_total[i]->type = SEGMENT_PARTICLE;
  //       else if (tag_characteristic == NORMAL_CELL && tag_interface == CUT_CELL)
  //         particle_total[i]->type = SURFACE_PARTICLE;
  //       else
  //         particle_total[i]->type = REAL_PARTICLE;
  //     }else{
  //       particle_total[i]->type = REAL_PARTICLE;
  //     }
    }
  #endif

    tet_delete.clear();

    for (it=tets.begin(); it!=tets.end(); ++it){
      int p1 = int((*it).i);
      int p2 = int((*it).j);
      int p3 = int((*it).k);
      int p4 = int((*it).l);

      p_Particle particle_1 = particle_total[p1];
      p_Particle particle_2 = particle_total[p2];
      p_Particle particle_3 = particle_total[p3];
      p_Particle particle_4 = particle_total[p4];

      bool error_tet = false;

      error_tet = Check_whether_this_tet_should_be_delete(2, particle_1, particle_2, particle_3, particle_4, iRank, sph);

      if (error_tet) tet_delete.push_back(it);

    }

    if (tet_delete.size()>0){
      int n_delete = tet_delete.size();
      cout<<"<<<<< ["<<world.rank()<<"] Num. of tets to be deleted: "<<tet_delete.size()<<endl;

      for (int i=0; i < n_delete; i++){
        tets.erase(tet_delete[i]);
      }

      if (tets.size() + n_delete != tet_out.numberoftetrahedra){
        cout<<"<<<<< ["<<world.rank()<<"] ERROR: number of tets are wrong!!!"<<endl;
      }
    }

  #ifndef _MPI
    pp_in.deinitialize();
    pp_in.initialize();
    tet_out.deinitialize();
    tet_out.initialize();
  #endif
  #ifdef _MPI_
    mpi::wait_all(reqs               +1, reqs                + nproc);
    mpi::wait_all(reqs_2nd_send_color+1, reqs_2nd_send_color + nproc);
    mpi::wait_all(reqs_2nd_send_coord+1, reqs_2nd_send_coord + nproc);
    
    cout<<"[0]-> all data sent"<<endl;
    
    parallel_for( blocked_range<int>(0, nproc),
              [&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        gather_coord  [i].Vector.clear();
        gather_coord  [i].mem_size = 0;
        gather_coord  [i].tag = 0;
        gather_tets   [i].Vector.clear();
        gather_tets   [i].mem_size = 0;
        gather_tets   [i].tag = 0;
        request_id    [i].Vector.clear();
        request_id    [i].mem_size = 0;
        request_id    [i].tag = 0;
        color_to_send [i].Vector.clear();
        color_to_send [i].mem_size = 0;
        color_to_send [i].tag = 0;
        coord_to_send [i].Vector.clear();
        coord_to_send [i].mem_size = 0;
        coord_to_send [i].tag = 0;
      }
    }, ap);
    
    delete [] gather_coord;
    delete [] gather_tets;
    delete [] request_id;
    delete [] color_to_send;
    delete [] coord_to_send;
    
    pp_in.deinitialize();
    pp_in.initialize();
    tet_out.deinitialize();
    tet_out.initialize();
    
  }else{
    out_coord.Vector.clear();
    out_coord.mem_size = 0;
    out_coord.tag = 0;
    out_tets.Vector.clear();
    out_tets.mem_size = 0;
    out_tets.tag = 0;
    
    int num_to_send = sph->total_num_particle;
    for (int i=0; i < num_to_send; i++){
      out_coord.Vector.push_back(particle_total[i]->coord.i);
      out_coord.Vector.push_back(particle_total[i]->coord.j);
      out_coord.Vector.push_back(particle_total[i]->coord.k);
    }

    for (it=tets.begin(); it!=tets.end(); ++it){
      out_tets.Vector.push_back(particle_total[int((*it).i)]->id);
      out_tets.Vector.push_back(particle_total[int((*it).j)]->id);
      out_tets.Vector.push_back(particle_total[int((*it).k)]->id);
      out_tets.Vector.push_back(particle_total[int((*it).l)]->id);
    }
    
    out_coord.mem_size = int(out_coord.Vector.size());
    out_coord.tag = 1;
    out_tets.mem_size = int(out_tets.Vector.size());
    out_tets.tag = 1;

    gather(world,out_coord,0);
    gather(world,out_tets,0);
        
    out_coord.Vector.clear();
    out_coord.mem_size = 0;
    out_coord.tag = 0;
    out_tets.Vector.clear();
    out_tets.mem_size = 0;
    out_tets.tag = 0;
    
    tets.clear();
    
    mpi::request reqs;
    
    reqs = world.irecv(0, iRank, out_tets);
    
    reqs.wait();
    
    cout<<"<<<<< ["<<world.rank()<<"] data received"<<endl;
    
    n_tets = int(out_tets.Vector.size()/4);
    std::vector<int> vert_list;
    
    // JZ20190102 get missing particles
    serialization_vector <int> required_particle_id;
    serialization_vector <int> required_particle_color;
    serialization_vector <Real> required_particle_coord;    
    required_particle_id.Vector.clear();
    required_particle_id.mem_size = 0;
    required_particle_id.tag = 0;
    required_particle_color.Vector.clear();
    required_particle_color.mem_size = 0;
    required_particle_color.tag = 0;
    required_particle_coord.Vector.clear();
    required_particle_coord.mem_size = 0;
    required_particle_coord.tag = 0;
    
    for (int i=0; i < n_tets; i++){
      for (int j = 0; j < 4; ++j)
      {
        int vert = out_tets.Vector[i*4 + j];
        if (local_id_table[vert] == -1){
          local_id_table[vert] = -2;
        }
      }
    }
    
    for (int i=0; i < totalnumberofpoints; i++){
      if (local_id_table[i] == -2)
        required_particle_id.Vector.push_back(i);
    }
    
    mpi::request reqs_2nd[3];
    
    if (required_particle_id.Vector.size() > 0){
      // JZ20180102::temporarily using ghost_particle for the tracking of extra particles sent
      if (sph->ghost_particle.size() > 0){
        for (int i = 0; i < sph->ghost_particle.size(); i++)
          sph->particlepool.free(sph->ghost_particle[i]);
        sph->ghost_particle.clear();
      }
      
      required_particle_id.mem_size = int(required_particle_id.Vector.size());
      required_particle_id.tag = 1;
      
      reqs_2nd[0] = world.isend(0, iRank,   required_particle_id);
      reqs_2nd[1] = world.irecv(0, iRank*2, required_particle_color);
      reqs_2nd[2] = world.irecv(0, iRank*3, required_particle_coord);

      mpi::wait_all(reqs_2nd, reqs_2nd+3);
      
      if (int(required_particle_coord.Vector.size()/3) != required_particle_id.Vector.size()){
        cout<<"<<<<< ["<<world.rank()<<"] ERROR!!! wrong number of received extra points"<<endl;
        world.abort(-1);
      }else
        cout<<"<<<<< ["<<world.rank()<<"] etra data ("<<int(required_particle_coord.Vector.size()/3)<<" point(s)) received"<<endl;
      
      int num_extra_points = required_particle_id.Vector.size();
      int local_id_record = sph->local_id_record;
      for (int i = 0; i < num_extra_points; i++){
        p_Particle p_extra = sph->particlepool.malloc();
        p_extra->coord.i  = required_particle_coord.Vector[i*3 + 0];
        p_extra->coord.j  = required_particle_coord.Vector[i*3 + 1];
        p_extra->coord.k  = required_particle_coord.Vector[i*3 + 2];
        p_extra->color    = required_particle_color.Vector[i];
        p_extra->id       = required_particle_id.Vector[i];
        p_extra->local_id = i+local_id_record;
        
        local_id_table[p_extra->id] = p_extra->local_id;
        
        int tag_interface      = NORMAL_CELL;
        int tag_characteristic = NORMAL_CELL;
        int idx_characteristic = -1;
        
        p_extra->Calculate_particle_infor(sph);
        // p_extra->Get_extended_char_cell_tag_and_idx (sph, tag_interface, tag_characteristic, idx_characteristic);
        p_extra->Get_level_set_char_cell_tag_and_idx (sph, tag_interface, tag_characteristic, idx_characteristic);
    
        if ( tag_interface != NORMAL_CELL ){
          if (tag_characteristic == SINGULARITY_CELL)
            p_extra->type = SINGULARITY_PARTICLE;
          else if (tag_characteristic == SEGMENT_CELL)
            p_extra->type = SEGMENT_PARTICLE;
          else if (tag_characteristic == NORMAL_CELL && tag_interface == CUT_CELL)
            p_extra->type = SURFACE_PARTICLE;
          else
            p_extra->type = REAL_PARTICLE;
        }else{
          p_extra->type = REAL_PARTICLE;
        }
        particle_total.push_back(p_extra);
        sph->ghost_particle.push_back(p_extra);
      }
      
      local_id_record += num_extra_points;
      
      required_particle_id.Vector.clear();
      required_particle_id.mem_size = 0;
      required_particle_id.tag = 0;
      required_particle_color.Vector.clear();
      required_particle_color.mem_size = 0;
      required_particle_color.tag = 0;
      required_particle_coord.Vector.clear();
      required_particle_coord.mem_size = 0;
      required_particle_coord.tag = 0;
      
    }else if (required_particle_id.Vector.size() == 0){
      required_particle_id.Vector.push_back(-1);
      required_particle_id.mem_size = 1;
      required_particle_id.tag = 1;
      
      reqs_2nd[0] = world.isend(0, iRank, required_particle_id);
//       reqs_2nd[0].wait();

      required_particle_id.Vector.clear();
      required_particle_id.mem_size = 0;
      required_particle_id.tag = 0;
    }
    
    for (int i=0; i < n_tets; i++){
      vert_list.clear();
      for (int j = 0; j < 4; ++j)
      {
        vert_list.push_back(out_tets.Vector[i*4 + j]);
      }

      int4_graph tet;
      for (int k = 0; k < 4; k++){
        if (local_id_table[vert_list[k]] == -1){
          cout<<"<<<<< ["<<world.rank()<<"] this global id "<<vert_list[k]<<" is not in the current rank"<<endl;
          world.abort(-1);
        }
      }
      tet.i      = local_id_table[vert_list[0]];
      tet.j      = local_id_table[vert_list[1]];
      tet.k      = local_id_table[vert_list[2]];
      tet.l      = local_id_table[vert_list[3]];
      tet.vol    = 0.;
      tet.aspect = 0.;
      tets.insert (tet);
    }

    out_tets.Vector.clear();
    out_tets.mem_size = 0;
    out_tets.tag = 0;
    
    tet_delete.clear();

    for (it=tets.begin(); it!=tets.end(); ++it){
      int p1 = int((*it).i);
      int p2 = int((*it).j);
      int p3 = int((*it).k);
      int p4 = int((*it).l);

      p_Particle particle_1 = particle_total[p1];
      p_Particle particle_2 = particle_total[p2];
      p_Particle particle_3 = particle_total[p3];
      p_Particle particle_4 = particle_total[p4];

      bool error_tet = false;

      error_tet = Check_whether_this_tet_should_be_delete(2,particle_1, particle_2, particle_3, particle_4, iRank, sph);

      if (error_tet) tet_delete.push_back(it);

    }

    if (tet_delete.size()>0){
      int n_delete = tet_delete.size();
      cout<<"<<<<< ["<<world.rank()<<"] Num. of tets to be deleted: "<<tet_delete.size()<<endl;

      for (int i=0; i < n_delete; i++){
        tets.erase(tet_delete[i]);
      }

      if (tets.size() + n_delete != n_tets){
        cout<<"<<<<< ["<<world.rank()<<"] ERROR: number of tets are wrong!!!"<<endl;
      }
    }
  }
  
  delete [] local_id_table;
  #endif

//   world.barrier();
//   world.abort(-1);
  #if !defined (_READ_TARGET_FIELD_FOR_POST_)
  }
  #endif
  
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Tet mesh quality section"<<endl;
  }
  
  Tet_quality_statistics (particle_total, sph, world);
  
  if (world.rank() == 0)
    cout<<"**********************************************************\n";

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Reconstruct_tets finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// Reconstruct::Tet_quality_statistics
//-------------------------------------------------------
void Graphmeshcls::Tet_quality_statistics(std::vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  tetgenmesh m;

  Real dratio = 0.02;
  Real dangle = 4.0;
  int num_radiustable = int(1.1/dratio)+1;
  int num_dihedangletable = int(180/dangle);
  // int num_faceangletable = int(180/dangle);

  int  _ntet_sm_10 = 0;
  int  _ntet_sm_20 = 0;
  int  _ntet_sm_30 = 0;
  int  _ntet_sm_40 = 0;

  int local_communication_tets = 0;

  tetgenmesh::triface tetloop, neightet;
  Real *p[4];
  char sbuf[128];
  Real *radiusratiotable;
  Real aspectratiotable[12];
  Real A[4][4], rhs[4], D;
  Real V[6][3], N[4][3], H[4]; // edge-vectors, face-normals, face-heights.
  Real edgelength[6], alldihed[6], faceangle[3];
  Real shortest, longest;
  Real smallestvolume, biggestvolume;
  Real smallestratio, biggestratio;
  Real smallestradiusratio, biggestradiusratio; // radius-edge ratio.
  Real smallestdiangle, biggestdiangle;
  // Real smallestfaangle, biggestfaangle;
  Real total_tet_vol, total_tetprism_vol;
  Real total_tet_aspect_radio, total_tet_radius, total_minimum_dihedangle;
  Real tetvol, minaltitude;
  Real cirradius, minheightinv, insradius;
  Real shortlen, longlen;
  Real tetaspect, tetradius;
  Real smalldiangle, bigdiangle;
  Real smallfaangle, bigfaangle;
  unsigned long *radiustable;
  unsigned long aspecttable[16];
  unsigned long *dihedangletable;
  // unsigned long *faceangletable;
  int indx[4];
  int radiusindex;
  int aspectindex;
  int tendegree;
  int i, j;
  // Report the tet which has the biggest radius-edge ratio.
  // triface biggestradiusratiotet;

  shortlen = longlen = 0.0;
  smalldiangle = bigdiangle = 0.0;
  total_tet_vol = 0.0;
  total_tetprism_vol = 0.0;
  total_tet_aspect_radio = 0.0;
  total_tet_radius = 0.0;
  total_minimum_dihedangle = 0.0;

  radiusratiotable = new Real         [num_radiustable];
  radiustable      = new unsigned long[num_radiustable];
  dihedangletable  = new unsigned long[num_dihedangletable];
  // faceangletable   = new unsigned long[num_faceangletable];

  for (int i = 0; i < num_radiustable-1; ++i)
  {
    radiusratiotable[i] = dratio*(i+1);
  }
  radiusratiotable[num_radiustable-1] = 999.;

  aspectratiotable[0]  =      1.5;    aspectratiotable[1]  =     2.0;
  aspectratiotable[2]  =      2.5;    aspectratiotable[3]  =     3.0;
  aspectratiotable[4]  =      4.0;    aspectratiotable[5]  =     6.0;
  aspectratiotable[6]  =     10.0;    aspectratiotable[7]  =    15.0;
  aspectratiotable[8]  =     25.0;    aspectratiotable[9]  =    50.0;
  aspectratiotable[10] =    100.0;    aspectratiotable[11] =     0.0;
  
  for (i = 0; i < num_radiustable    ; i++) radiustable[i]     = 0l;
  for (i = 0; i < 12                 ; i++) aspecttable[i]     = 0l;
  for (i = 0; i < num_dihedangletable; i++) dihedangletable[i] = 0l;
  // for (i = 0; i < num_faceangletable ; i++) faceangletable[i]  = 0l;

  my_self_add (sph->domain, minaltitude);
  minaltitude = minaltitude * minaltitude;
  shortest = minaltitude;
  longest = 0.0;
  smallestvolume = minaltitude;
  biggestvolume = 0.0;
  smallestratio = smallestradiusratio = 1e+16; // minaltitude;
  biggestratio = biggestradiusratio = 0.0;
  smallestdiangle = 180.;
  // smallestfaangle = 180.0;
  biggestdiangle = 0.0;
  // biggestfaangle = 0.0;

  // Loop all elements, calculate quality parameters for each element.

  std::set<int4_graph>::iterator it;
  for (it=tets.begin(); it!=tets.end(); ++it){

    int p_[4];
    p_[0] = int((*it).i);
    p_[1] = int((*it).j);
    p_[2] = int((*it).k);
    p_[3] = int((*it).l);

    p_Particle part_[4];
    Real pp[12];
    for (i = 0; i < 4; i++) {
      part_[i      ] = particle_total[p_[i]];
      pp   [i*3 + 0] = part_[i]->coord.i;
      pp   [i*3 + 1] = part_[i]->coord.j;
      pp   [i*3 + 2] = part_[i]->coord.k;
    }

    #if defined (_READ_TARGET_FIELD_FOR_POST_)
    std::vector<int> color_list;
    int color_count = 1;
    for (i = 0; i < 4; ++i)
    {
      color_list.push_back(part_[i]->color);
    }
    std::sort (color_list.begin(), color_list.end(), Greater());

    for (i = 0; i < 3; ++i)
    {
      if(color_list[i] != color_list[i+1]) color_count++;
    }
    if (color_count > 1) local_communication_tets++;
    #endif

    // Get four vertices: p0, p1, p2, p3.
    for (i = 0; i < 4; i++) {
      p[i] = &pp[i*3];
    }
    // Get the tet volume.
    tetvol = m.orient3dfast(p[1], p[0], p[2], p[3]) / 6.0;
    (*it).vol = tetvol;
    total_tet_vol += tetvol;
    total_tetprism_vol += m.tetprismvol(p[0], p[1], p[2], p[3]);

    // Calculate the largest and smallest volume.
    if (tetvol < smallestvolume) {
      smallestvolume = tetvol;
    } 
    if (tetvol > biggestvolume) {
      biggestvolume = tetvol;
    }

    // Set the edge vectors: V[0], ..., V[5]
    for (i = 0; i < 3; i++) V[0][i] = p[0][i] - p[3][i]; // V[0]: p3->p0.
    for (i = 0; i < 3; i++) V[1][i] = p[1][i] - p[3][i]; // V[1]: p3->p1.
    for (i = 0; i < 3; i++) V[2][i] = p[2][i] - p[3][i]; // V[2]: p3->p2.
    for (i = 0; i < 3; i++) V[3][i] = p[1][i] - p[0][i]; // V[3]: p0->p1.
    for (i = 0; i < 3; i++) V[4][i] = p[2][i] - p[1][i]; // V[4]: p1->p2.
    for (i = 0; i < 3; i++) V[5][i] = p[0][i] - p[2][i]; // V[5]: p2->p0.

    // Get the squares of the edge lengths.
    for (i = 0; i < 6; i++) edgelength[i] = m.dot(V[i], V[i]);

    // Calculate the longest and shortest edge length.
    for (i = 0; i < 6; i++) {
      if (i == 0) {
        shortlen = longlen = edgelength[i];
      } else {
        shortlen = edgelength[i] < shortlen ? edgelength[i] : shortlen;
        longlen  = edgelength[i] > longlen  ? edgelength[i] : longlen;
      }
      if (edgelength[i] > longest) {
        longest = edgelength[i];
      } 
      if (edgelength[i] < shortest) {
        shortest = edgelength[i];
      }
    }

    // Set the matrix A = [V[0], V[1], V[2]]^T.
    for (j = 0; j < 3; j++) {
      for (i = 0; i < 3; i++) A[j][i] = V[j][i];
    }

    // Decompose A just once.
    if (m.lu_decmp(A, 3, indx, &D, 0)) {   
      // Get the three faces normals.
      for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++) rhs[i] = 0.0;
        rhs[j] = 1.0;  // Positive means the inside direction
        m.lu_solve(A, 3, indx, rhs, 0);
        for (i = 0; i < 3; i++) N[j][i] = rhs[i];
      }
      // Get the fourth face normal by summing up the first three.
      for (i = 0; i < 3; i++) N[3][i] = - N[0][i] - N[1][i] - N[2][i];
      // Get the radius of the circumsphere.
      for (i = 0; i < 3; i++) rhs[i] = 0.5 * m.dot(V[i], V[i]);
      m.lu_solve(A, 3, indx, rhs, 0);
      cirradius = sqrt(m.dot(rhs, rhs));
      // Normalize the face normals.
      for (i = 0; i < 4; i++) {
        // H[i] is the inverse of height of its corresponding face.
        H[i] = sqrt(m.dot(N[i], N[i]));
        for (j = 0; j < 3; j++) N[i][j] /= H[i];
      }
      // Get the radius of the inscribed sphere.
      insradius = 1.0 / (H[0] + H[1] + H[2] + H[3]);
      // Get the biggest H[i] (corresponding to the smallest height).
      minheightinv = H[0];
      for (i = 1; i < 4; i++) {
        if (H[i] > minheightinv) minheightinv = H[i];
      }
    } else {
      // A nearly degenerated tet.
      if (tetvol < 0.0) {
        cout<<"<<<<< ["<<world.rank()<<"] !! Warning:  A inverted tet ["<<p_[0]<<" "<<p_[1]<<" "<<p_[2]<<" "<<p_[3]<<" ]"<<endl;
        continue;
      }else if (tetvol == 0.0){
        cout<<"<<<<< ["<<world.rank()<<"] !! Warning:  A degenerated tet ["<<p_[0]<<" "<<p_[1]<<" "<<p_[2]<<" "<<p_[3]<<" ]"<<endl;
      }
      // Calculate the four face normals.
      m.facenormal(p[2], p[1], p[3], N[0], 1, NULL);
      m.facenormal(p[0], p[2], p[3], N[1], 1, NULL);
      m.facenormal(p[1], p[0], p[3], N[2], 1, NULL);
      m.facenormal(p[0], p[1], p[2], N[3], 1, NULL);
      // Normalize the face normals.
      for (i = 0; i < 4; i++) {
        // H[i] is the twice of the area of the face.
        H[i] = sqrt(m.dot(N[i], N[i]));
        for (j = 0; j < 3; j++) N[i][j] /= H[i];
      }
      // Get the biggest H[i] / tetvol (corresponding to the smallest height).
      minheightinv = (H[0] / tetvol);
      for (i = 1; i < 4; i++) {
        if ((H[i] / tetvol) > minheightinv) minheightinv = (H[i] / tetvol);
      }
      // Let the circumradius to be the half of its longest edge length.
      cirradius = 0.5 * sqrt(longlen);
    }

    // Get the dihedrals (in degree) at each edges.
    j = 0;
    for (i = 1; i < 4; i++) {
      alldihed[j] = -m.dot(N[0], N[i]); // Edge cd, bd, bc.
      if (alldihed[j] < -1.0) alldihed[j] = -1; // Rounding.
      else if (alldihed[j] > 1.0) alldihed[j] = 1;
      alldihed[j] = acos(alldihed[j]) / PI * 180.0;
      j++;
    }
    for (i = 2; i < 4; i++) {
      alldihed[j] = -m.dot(N[1], N[i]); // Edge ad, ac.
      if (alldihed[j] < -1.0) alldihed[j] = -1; // Rounding.
      else if (alldihed[j] > 1.0) alldihed[j] = 1;
      alldihed[j] = acos(alldihed[j]) / PI * 180.0;
      j++;
    }
    alldihed[j] = -m.dot(N[2], N[3]); // Edge ab.
    if (alldihed[j] < -1.0) alldihed[j] = -1; // Rounding.
    else if (alldihed[j] > 1.0) alldihed[j] = 1;
    alldihed[j] = acos(alldihed[j]) / PI * 180.0;

    // Calculate the largest and smallest dihedral angles.
    for (i = 0; i < 6; i++) {
      if (i == 0) {
        smalldiangle = bigdiangle = alldihed[i];
      } else {
        smalldiangle = alldihed[i] < smalldiangle ? alldihed[i] : smalldiangle;
        bigdiangle = alldihed[i] > bigdiangle ? alldihed[i] : bigdiangle;
      }
      if (alldihed[i] < smallestdiangle) {
        smallestdiangle = alldihed[i];
      } 
      if (alldihed[i] > biggestdiangle) {
        biggestdiangle = alldihed[i];
      }
      // Accumulate the corresponding number in the dihedral angle histogram.
      tendegree = floor(alldihed[i] / dangle);
      dihedangletable[tendegree]++;

      // if (alldihed[i] < 5.0) {
      //   tendegree = 0;
      // } else if (alldihed[i] >= 5.0 && alldihed[i] < 10.0) {
      //   tendegree = 1;
      // } else if (alldihed[i] >= 80.0 && alldihed[i] < 110.0) {
      //   tendegree = 9; // Angles between 80 to 110 degree are in one entry.
      // } else if (alldihed[i] >= 170.0 && alldihed[i] < 175.0) {
      //   tendegree = 16;
      // } else if (alldihed[i] >= 175.0) {
      //   tendegree = 17;
      // } else {
      //   tendegree = (int) (alldihed[i] / 10.);
      //   if (alldihed[i] < 80.0) {
      //     tendegree++;  // In the left column.
      //   } else {
      //     tendegree--;  // In the right column.
      //   }
      // }
      // dihedangletable[tendegree]++;
    }

    (*it).mindihedangle = smalldiangle;
    if (smalldiangle < 10.) _ntet_sm_10++;
    if (smalldiangle < 20.) _ntet_sm_20++;
    if (smalldiangle < 30.) _ntet_sm_30++;
    if (smalldiangle < 40.) _ntet_sm_40++;

    total_minimum_dihedangle += smalldiangle;

    // int3_graph tri_[3];
    // tri_[0].i = p_[0];
    // tri_[0].j = p_[1];
    // tri_[0].k = p_[2];

    // tri_[1].i = p_[0];
    // tri_[1].j = p_[1];
    // tri_[1].k = p_[3];

    // tri_[2].i = p_[1];
    // tri_[2].j = p_[2];
    // tri_[2].k = p_[3];

    // // Calculate the largest and smallest face angles.
    // for (i = 0; i < 3; i++){
    //   my_real vertex_0, vertex_1, vertex_2;

    //   // obtain the coordinates
    //   my_set_data (vertex_0, particle_total[tri_[i].i]->coord);
    //   my_set_data (vertex_1, particle_total[tri_[i].j]->coord);
    //   my_set_data (vertex_2, particle_total[tri_[i].k]->coord);

    //   // obtain the vectors
    //   my_real vector_1 = my_minus_data (vertex_1, vertex_0);
    //   Real      dist_1 = get_distance  (vector_1);

    //   my_real vector_2 = my_minus_data (vertex_2, vertex_1);
    //   Real      dist_2 = get_distance  (vector_2);

    //   my_real vector_3 = my_minus_data (vertex_0, vertex_2);
    //   Real      dist_3 = get_distance  (vector_3);

    //   Real alpha_1 = 180./acos(-1.) * acos(-(vector_1.i * vector_3.i + vector_1.j * vector_3.j)/(dist_3*dist_1 + 1.e-20));
    //   Real alpha_2 = 180./acos(-1.) * acos(-(vector_1.i * vector_2.i + vector_1.j * vector_2.j)/(dist_2*dist_1 + 1.e-20));
    //   Real alpha_3 = 180./acos(-1.) * acos(-(vector_2.i * vector_3.i + vector_2.j * vector_3.j)/(dist_3*dist_2 + 1.e-20));

    //   alpha_1 = AMAX1(0., AMIN1(alpha_1, 180.));
    //   alpha_2 = AMAX1(0., AMIN1(alpha_2, 180.));
    //   alpha_3 = AMAX1(0., AMIN1(alpha_3, 180.));

    //   tendegree = (int) (alpha_1 / dangle);
    //   faceangletable[tendegree]++;

    //   tendegree = (int) (alpha_2 / dangle);
    //   faceangletable[tendegree]++;

    //   tendegree = (int) (alpha_3 / dangle);
    //   faceangletable[tendegree]++;

    //   biggestfaangle = AMAX1(biggestfaangle, AMAX1(alpha_1, AMAX1(alpha_2, alpha_3)));
    //   smallestfaangle = AMIN1(smallestfaangle, AMIN1(alpha_1, AMIN1(alpha_2, alpha_3)));
    // }

    // Calculate aspect ratio and radius-edge ratio for this element.
    // tetradius = cirradius / sqrt(shortlen);
    tetradius = 3.0* insradius/ cirradius;
    (*it).radius_ratio = tetradius;
    total_tet_radius += tetradius;
    if (tetradius < smallestradiusratio) {
      smallestradiusratio = tetradius;
    }
    if (tetradius > biggestradiusratio) {
      biggestradiusratio = tetradius;
      // biggestradiusratiotet.tet = tetloop.tet;
    }
    // tetaspect = sqrt(longlen) / (2.0 * insradius);
    tetaspect = sqrt(longlen) * minheightinv;
    (*it).aspect = tetaspect;
    total_tet_aspect_radio += tetaspect;
    // Remember the largest and smallest aspect ratio.
    if (tetaspect < smallestratio) {
      smallestratio = tetaspect;
    } 
    if (tetaspect > biggestratio) {
      biggestratio = tetaspect;
    }
    // Accumulate the corresponding number in the aspect ratio histogram.
    aspectindex = 0;
    while ((tetaspect > aspectratiotable[aspectindex]) && (aspectindex < 11)) {
      aspectindex++;
    }
    aspecttable[aspectindex]++;

    int radiusindex = floor (tetradius/dratio);

    if (radiusindex < num_radiustable-1)
      radiustable[radiusindex]++;
    else radiustable[num_radiustable-1]++;

    // radiusindex = 0;
    // while ((tetradius > radiusratiotable[radiusindex]) && (radiusindex < num_radiustable-1)) {
    //   radiusindex++;
    // }
    // radiustable[radiusindex]++;
  }
  
  if (world.rank() == 0){
    Real _shortest              , _longest           ;
    Real _smallestvolume        , _biggestvolume     ;
    Real _smallestratio         , _biggestratio      ;
    Real _smallestradiusratio   , _biggestradiusratio; // radius-edge ratio.
    Real _smallestdiangle       , _biggestdiangle    ;
    // Real _smallestfaangle       , _biggestfaangle    ;
    Real _total_tet_vol         , _total_tetprism_vol;
    Real _total_tet_aspect_radio, _total_tet_radius  , _total_minimum_dihedangle;
    unsigned long *_radiustable;
    unsigned long _aspecttable[16];
    unsigned long *_dihedangletable;
    // unsigned long *_faceangletable;

    _radiustable     = new unsigned long [num_radiustable];
    _dihedangletable = new unsigned long [num_dihedangletable];
    // _faceangletable  = new unsigned long [num_faceangletable];
  
    glbl_num_tets             = 0;
    ntet_sm_10                = 0;
    ntet_sm_20                = 0;
    ntet_sm_30                = 0;
    ntet_sm_40                = 0;
    glbl_communication_tets   = 0;
    glbl_communication_volume = 0.;
    _shortest                 = minaltitude;
    _longest                  = 0.0;
    _total_tet_vol            = 0.0;
    _total_tetprism_vol       = 0.0;
    _total_tet_aspect_radio   = 0.0;
    _total_tet_radius         = 0.0;
    _total_minimum_dihedangle = 0.0;
    _smallestvolume           = minaltitude;
    _biggestvolume            = 0.0;
    
    _smallestratio   = _smallestradiusratio = 1e+16; // minaltitude;
    _biggestratio    = _biggestradiusratio  = 0.0;
    _smallestdiangle = 180.;
    // _smallestfaangle     = 180.0;
    _biggestdiangle  = 0.0;
    // _biggestfaangle      = 0.0;
    
    for (i = 0; i < num_radiustable    ; i++) _radiustable[i]     = 0l;
    for (i = 0; i < 12                 ; i++) _aspecttable[i]     = 0l;
    // for (i = 0; i < 18                 ; i++) _faceangletable[i]  = 0l;
    for (i = 0; i < num_dihedangletable; i++) _dihedangletable[i] = 0l;
    
    reduce( world, total_tet_vol,            _total_tet_vol,                   std::plus<Real>()      , 0);
    reduce( world, total_tet_radius,         _total_tet_radius,                std::plus<Real>()      , 0);
    reduce( world, total_tetprism_vol,       _total_tetprism_vol,              std::plus<Real>()      , 0);
    reduce( world, total_tet_aspect_radio,   _total_tet_aspect_radio,          std::plus<Real>()      , 0);
    reduce( world, total_minimum_dihedangle, _total_minimum_dihedangle,        std::plus<Real>()      , 0);

    reduce( world, shortest,                 _shortest,                     mpi::minimum<Real>()      , 0);
    reduce( world, smallestratio,            _smallestratio,                mpi::minimum<Real>()      , 0);
    reduce( world, smallestvolume,           _smallestvolume,               mpi::minimum<Real>()      , 0);
    reduce( world, smallestdiangle,          _smallestdiangle,              mpi::minimum<Real>()      , 0);
    // reduce( world, smallestfaangle,          _smallestfaangle,              mpi::minimum<Real>()      , 0);
    reduce( world, smallestradiusratio,      _smallestradiusratio,          mpi::minimum<Real>()      , 0);

    reduce( world, longest,                  _longest,                      mpi::maximum<Real>()      , 0);
    reduce( world, biggestratio,             _biggestratio,                 mpi::maximum<Real>()      , 0);
    reduce( world, biggestvolume,            _biggestvolume,                mpi::maximum<Real>()      , 0);
    reduce( world, biggestdiangle,           _biggestdiangle,               mpi::maximum<Real>()      , 0);
    // reduce( world, biggestfaangle,           _biggestfaangle,               mpi::maximum<Real>()      , 0);
    reduce( world, biggestradiusratio,       _biggestradiusratio,           mpi::maximum<Real>()      , 0);
    
    reduce( world, int(tets.size()),          glbl_num_tets,                std::plus<int>()          , 0);
    reduce( world, _ntet_sm_10     ,          ntet_sm_10   ,                std::plus<int>()          , 0);
    reduce( world, _ntet_sm_20     ,          ntet_sm_20   ,                std::plus<int>()          , 0);
    reduce( world, _ntet_sm_30     ,          ntet_sm_30   ,                std::plus<int>()          , 0);
    reduce( world, _ntet_sm_40     ,          ntet_sm_40   ,                std::plus<int>()          , 0);
    reduce( world, local_communication_tets, glbl_communication_tets,       std::plus<int>()          , 0);

    reduce( world, radiustable,      num_radiustable,     _radiustable,     std::plus<unsigned long>(), 0);
    reduce( world, aspecttable,      16,                  _aspecttable,     std::plus<unsigned long>(), 0);
    // reduce( world, faceangletable,   num_faceangletable,  _faceangletable,  std::plus<unsigned long>(), 0);
    reduce( world, dihedangletable,  num_dihedangletable, _dihedangletable, std::plus<unsigned long>(), 0);

    _shortest    = sqrt(_shortest);
    _longest     = sqrt(_longest);

    maximum_dihedangle       = _biggestdiangle;
    minimum_dihedangle       = _smallestdiangle;
    minimum_dihedangle_avg   = _total_minimum_dihedangle/Real(glbl_num_tets);
    maximum_tet_aspect_ratio = _biggestratio;
    minimum_tet_aspect_ratio = _smallestratio;
    average_tet_aspect_ratio = _total_tet_aspect_radio/Real(glbl_num_tets);
    maximum_tetradius        = _biggestradiusratio;
    minimum_tetradius        = _smallestradiusratio;
    average_tetradius        = _total_tet_radius/Real(glbl_num_tets);
    maximum_volume           = _biggestvolume;
    minimum_volume           = _smallestvolume;
    average_volume           = _total_tet_vol/Real(glbl_num_tets);
    #ifdef _READ_TARGET_FIELD_FOR_POST_
    glbl_communication_volume = Real(glbl_communication_tets)/Real(glbl_num_tets);
    #endif
    printf("  Total munber of tets: %16d\n", glbl_num_tets);

    printf("  Smallest volume: %16.5g   |  Largest volume: %16.5g   |  Average volume: %16.5g\n",
          _smallestvolume, _biggestvolume, average_volume);
    printf("  Shortest edge:   %16.5g   |  Longest edge:   %16.5g\n",
          _shortest, _longest);
    printf("  Smallest rad.ratio: %13.5g   |  Largest rad.ratio: %13.5g   |  Average rad.ratio: %13.5g\n",
          _smallestradiusratio, _biggestradiusratio, average_tetradius);
    printf("  Smallest asp.ratio: %13.5g   |  Largest asp.ratio: %13.5g   |  Average asp.ratio: %13.5g\n",
          _smallestratio, _biggestratio, average_tet_aspect_ratio);
    // sprintf(sbuf, "%.17g", _biggestfaangle);
    // if (strlen(sbuf) > 8) {
    //   sbuf[8] = '\0';
    // }
    // printf("  Smallest facangle: %14.5g   |  Largest facangle:       %s\n",
    //       _smallestfaangle, sbuf);
    // sprintf(sbuf, "%.17g", _biggestdiangle);
    // if (strlen(sbuf) > 8) {
    //   sbuf[8] = '\0';
    // }
    printf("  Smallest dihedral: %14.5g   |  Largest dihedral: %14.5g   |  Min. Avg. dihedral: %14.5g\n",
          _smallestdiangle, _biggestdiangle, minimum_dihedangle_avg);

    printf("  tets < 10        : %8d      |  tets < 20       :       %8d\n",
          ntet_sm_10, ntet_sm_20);
    printf("  tets < 30        : %8d      |  tets < 40       :       %8d\n",
          ntet_sm_30, ntet_sm_40);
    printf("\n");

    printf("  Radius ratio histogram:\n");

    for (i = 0; i < num_radiustable/2; i++) {
      Real st1, st2, ed1, ed2, pp1, pp2;
      unsigned long rad1, rad2;
      st1 = (i                    )*dratio;
      st2 = (i  +num_radiustable/2)*dratio;
      ed1 = (i+1                  )*dratio;
      ed2 = (i+1+num_radiustable/2)*dratio;
      rad1 = _radiustable[i                  ];
      rad2 = _radiustable[i+num_radiustable/2];
      pp1  = Real(rad1)/Real(glbl_num_tets)*100.;
      pp2  = Real(rad2)/Real(glbl_num_tets)*100.;

      if (i+1+num_radiustable/2 == num_radiustable) ed2 = 999.;

      printf("  %6.6g - %-6.6g    :  %8ld (%03.2f %%)    | %6.6g - %-6.6g     :  %8ld (%03.2f %%)\n",
             st1, ed1, rad1, pp1, st2, ed2, rad2, pp2);
    }
    if (num_radiustable%2 == 1){
      Real rad = _radiustable[num_radiustable-1];
      Real pp  = Real(_radiustable[num_radiustable-1])/Real(glbl_num_tets)*100.;
      printf("  %6.6g -  999.     :  %8ld (%03.2f %%)\n", 1.1, rad, pp);
    }
    printf("\n");

    // printf("  Radius ratio histogram:\n");
    // printf("         < %-6.6g    :  %8ld      | %6.6g - %-6.6g     :  %8ld\n",
    //       radiusratiotable[0], _radiustable[0], radiusratiotable[5],
    //       radiusratiotable[6], _radiustable[6]);
    // for (i = 1; i < 5; i++) {
    //   printf("  %6.6g - %-6.6g    :  %8ld      | %6.6g - %-6.6g     :  %8ld\n",
    //         radiusratiotable[i - 1], radiusratiotable[i], _radiustable[i],
    //         radiusratiotable[i + 5], radiusratiotable[i + 6],
    //         _radiustable[i + 6]);
    // }
    // printf("  %6.6g - %-6.6g    :  %8ld      | %6.6g -            :  %8ld\n",
    //       radiusratiotable[4], radiusratiotable[5], _radiustable[5],
    //       radiusratiotable[10], _radiustable[11]);

    printf("  Aspect ratio histogram:\n");
    printf("         < %-6.6g    :  %8ld      | %6.6g - %-6.6g     :  %8ld\n",
          aspectratiotable[0], _aspecttable[0], aspectratiotable[5],
          aspectratiotable[6], _aspecttable[6]);
    for (i = 1; i < 5; i++) {
      printf("  %6.6g - %-6.6g    :  %8ld      | %6.6g - %-6.6g     :  %8ld\n",
            aspectratiotable[i - 1], aspectratiotable[i], _aspecttable[i],
            aspectratiotable[i + 5], aspectratiotable[i + 6],
            _aspecttable[i + 6]);
    }
    printf("  %6.6g - %-6.6g    :  %8ld      | %6.6g -            :  %8ld\n",
          aspectratiotable[4], aspectratiotable[5], _aspecttable[5],
          aspectratiotable[10], _aspecttable[11]);
    printf("  (A tetrahedron's aspect ratio is its longest edge length");
    printf(" divided by its\n");
    printf("    smallest side height)\n");

    // printf("  Face angle histogram:\n");
    // for (int i = 0; i < num_faceangletable/2; ++i)
    // {
    //   Real st1, st2, ed1, ed2, pp1, pp2;
    //   unsigned long ang1, ang2;
    //   st1 = Real(i                       )*dangle;
    //   st2 = Real(i  +num_faceangletable/2)*dangle;
    //   ed1 = Real(i+1                     )*dangle;
    //   ed2 = Real(i+1+num_faceangletable/2)*dangle;
    //   ang1 = _faceangletable[i                  ];
    //   ang2 = _faceangletable[i+num_faceangletable/2];
    //   pp1  = Real(ang1)/Real(glbl_num_tets*6)*100.;
    //   pp2  = Real(ang2)/Real(glbl_num_tets*6)*100.;

    //   printf("     %4.4g - %4.4g degrees:  %8ld (%03.2f %%)   |    %4.4g - %4.4g degrees:  %8ld (%03.2f %%) \n",
    //           st1, ed1, ang1, pp1, st2, ed2, ang2, pp2);
    // }
    // for (i = 0; i < 9; i++) {
    //   printf("    %3d - %3d degrees:  %8ld      |    %3d - %3d degrees:  %8ld\n",
    //         i * 10, i * 10 + 10, _faceangletable[i],
    //         i * 10 + 90, i * 10 + 100, _faceangletable[i + 9]);
    // }
    printf("\n");

    printf("  Dihedral angle histogram:\n");
    // Print the three two rows:
    for (int i = 0; i < num_dihedangletable/2; ++i)
    {
      Real st1, st2, ed1, ed2, pp1, pp2;
      unsigned long ang1, ang2;
      st1 = Real(i                        )*dangle;
      st2 = Real(i  +num_dihedangletable/2)*dangle;
      ed1 = Real(i+1                      )*dangle;
      ed2 = Real(i+1+num_dihedangletable/2)*dangle;
      ang1 = _dihedangletable[i                  ];
      ang2 = _dihedangletable[i+num_dihedangletable/2];
      pp1  = Real(ang1)/Real(glbl_num_tets*6)*100.;
      pp2  = Real(ang2)/Real(glbl_num_tets*6)*100.;

      printf("     %4.4g - %4.4g degrees:  %8ld (%03.2f %%)   |    %4.4g - %4.4g degrees:  %8ld (%03.2f %%) \n",
              st1, ed1, ang1, pp1, st2, ed2, ang2, pp2);
    }

    // // Print the three two rows:
    // printf("     %3d - %2d degrees:  %8ld      |    %3d - %3d degrees:  %8ld\n",
    //       0, 5, _dihedangletable[0], 80, 110, _dihedangletable[9]);
    // printf("     %3d - %2d degrees:  %8ld      |    %3d - %3d degrees:  %8ld\n",
    //       5, 10, _dihedangletable[1], 110, 120, _dihedangletable[10]);
    // // Print the third to seventh rows.
    // for (i = 2; i < 7; i++) {
    //   printf("     %3d - %2d degrees:  %8ld      |    %3d - %3d degrees:  %8ld\n",
    //         (i - 1) * 10, (i - 1) * 10 + 10, _dihedangletable[i],
    //         (i - 1) * 10 + 110, (i - 1) * 10 + 120, _dihedangletable[i + 9]);
    // }
    // // Print the last two rows.
    // printf("     %3d - %2d degrees:  %8ld      |    %3d - %3d degrees:  %8ld\n",
    //       60, 70, _dihedangletable[7], 170, 175, _dihedangletable[16]);
    // printf("     %3d - %2d degrees:  %8ld      |    %3d - %3d degrees:  %8ld\n",
    //       70, 80, _dihedangletable[8], 175, 180, _dihedangletable[17]);
    printf("\n");

    printf("\n");

    char    dihang_histo[256];
    sprintf(dihang_histo,"%s%d%s","./outdata/dihedangle_histogram_step_",int(sph->run_time),".csv");
    ofstream dihang_histo_out(dihang_histo, ios::trunc);

    dihang_histo_out<<"\"degree\",\"number\",\"percentage\"\n";
    for (int i = 0; i < num_dihedangletable; ++i)
    {
      Real           st1 = Real(i+1)*dangle;
      unsigned long ang1 = _dihedangletable[i];
      Real           pp1 = Real(ang1)/Real(glbl_num_tets*6)*100.;
      dihang_histo_out<<st1<<","<<ang1<<","<<pp1<<endl;
    }
    dihang_histo_out.close();

    char    radius_histo[256];
    sprintf(radius_histo,"%s%d%s","./outdata/radiousratio_histogram_step_",int(sph->run_time),".csv");
    ofstream radius_histo_out(radius_histo, ios::trunc);

    radius_histo_out<<"\"ratio\",\"number\",\"percentage\"\n";
    for (int i = 0; i < num_radiustable-1; ++i)
    {
      Real           st1 = Real(i+1)*dratio;
      unsigned long rad1 = _radiustable[i];
      Real           pp1 = Real(rad1)/Real(glbl_num_tets)*100.;
      radius_histo_out<<st1<<","<<rad1<<","<<pp1<<endl;
    }
    Real pp1 = Real(_radiustable[num_radiustable-1])/Real(glbl_num_tets)*100.;
    radius_histo_out<<999<<","<<_radiustable[num_radiustable-1]<<","<<pp1<<endl;
    radius_histo_out.close();

  }else{
    reduce( world, total_tet_vol,            std::plus<Real>()         , 0);
    reduce( world, total_tet_radius,         std::plus<Real>()         , 0);
    reduce( world, total_tetprism_vol,       std::plus<Real>()         , 0);
    reduce( world, total_tet_aspect_radio,   std::plus<Real>()         , 0);
    reduce( world, total_minimum_dihedangle, std::plus<Real>()         , 0);

    reduce( world, shortest,                 mpi::minimum<Real>()      , 0);
    reduce( world, smallestratio,            mpi::minimum<Real>()      , 0);
    reduce( world, smallestvolume,           mpi::minimum<Real>()      , 0);
    reduce( world, smallestdiangle,          mpi::minimum<Real>()      , 0);
    // reduce( world, smallestfaangle,          mpi::minimum<Real>()      , 0);
    reduce( world, smallestradiusratio,      mpi::minimum<Real>()      , 0);

    reduce( world, longest,                  mpi::maximum<Real>()      , 0);
    reduce( world, biggestratio,             mpi::maximum<Real>()      , 0);
    reduce( world, biggestvolume,            mpi::maximum<Real>()      , 0);
    reduce( world, biggestdiangle,           mpi::maximum<Real>()      , 0);
    // reduce( world, biggestfaangle,           mpi::maximum<Real>()      , 0);
    reduce( world, biggestradiusratio,       mpi::maximum<Real>()      , 0);
    
    reduce( world, int(tets.size()),         std::plus<int>()          , 0);
    reduce( world, _ntet_sm_10     ,         std::plus<int>()          , 0);
    reduce( world, _ntet_sm_20     ,         std::plus<int>()          , 0);
    reduce( world, _ntet_sm_30     ,         std::plus<int>()          , 0);
    reduce( world, _ntet_sm_40     ,         std::plus<int>()          , 0);
    reduce( world, local_communication_tets, std::plus<int>()          , 0);

    reduce( world, radiustable,      num_radiustable,     std::plus<unsigned long>(), 0);
    reduce( world, aspecttable,      16,                  std::plus<unsigned long>(), 0);
    // reduce( world, faceangletable,   num_faceangletable,  std::plus<unsigned long>(), 0);
    reduce( world, dihedangletable,  num_dihedangletable, std::plus<unsigned long>(), 0);
  }

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Tet_quality_statistics finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
// Reconstruct::Write_tri_quality_plt
//-------------------------------------------------------
void Graphmeshcls::Write_tri_quality_plt  (int n, Real run_time, communicator &world)
{
  if (world.rank() == 0){
    char    filename[256];
    sprintf(filename,"%s","./outdata/tri_infor.csv");
    ofstream out(filename, ios::app);
    if (n == 0){
      out<<"\"run_time\",\"n_tris\",\"n_communicate_tris\",\"communication_vol\",\"maximum_angle\",\"minimum_angle\",\"minimum_angle_avg\",\"ntri_sm_10\",\"ntri_sm_20\",\"ntri_sm_30\",\"ntri_sm_40\",\"maximum_tri_aspect_ratio\",\"maximum_edge\",\"minimum_edge\",\"maximum_area\",\"minimum_area\",\"G_min\",\"G_avg\"\n";
    }

    out<<run_time<<","<<glbl_num_tris<<","<<glbl_communication_tris<<","<<glbl_communication_volume<<","<<maximum_angle<<","<<minimum_angle<<","<<minimum_angle_avg<<","<<ntri_sm_10<<","<<ntri_sm_20<<","<<ntri_sm_30<<","<<ntri_sm_40<<","<<maximum_tri_aspect_ratio<<","<<maximum_edge<<","<<minimum_edge<<","<<maximum_area<<","<<minimum_area<<","<<G_min<<","<<G_avg<<"\n";
  }
}
//-------------------------------------------------------
// Reconstruct::Write_tet_quality_plt
//-------------------------------------------------------
void Graphmeshcls::Write_tet_quality_plt  (int n, Real run_time, communicator &world)
{
  if (world.rank() == 0){
    char    filename[256];
    sprintf(filename,"%s","./outdata/tet_infor.csv");
    ofstream out(filename, ios::app);
    if (n == 0){
      out<<"\"run_time\",\"n_tets\",\"n_communicate_tets\",\"communication_vol\",\"maximum_dihedangle\",\"minimum_dihedangle\",\"minimum_dihedangle_avg\",\"ntet_sm_10\",\"ntet_sm_20\",\"ntet_sm_30\",\"ntet_sm_40\",\"minimum_tet_aspect_ratio\",\"maximum_tet_aspect_ratio\",\"average_tet_aspect_ratio\",\"minimum_tetradius\",\"maximum_tetradius\",\"average_tetradius\",\"maximum_volume\",\"minimum_volume\",\"average_volume\"\n";
    }

    out<<run_time<<","<<glbl_num_tets<<","<<glbl_communication_tets<<","<<glbl_communication_volume<<","<<maximum_dihedangle<<","<<minimum_dihedangle<<","<<minimum_dihedangle_avg<<","<<ntet_sm_10<<","<<ntet_sm_20<<","<<ntet_sm_30<<","<<ntet_sm_40<<","<<minimum_tet_aspect_ratio<<","<<maximum_tet_aspect_ratio<<","<<average_tet_aspect_ratio<<","<<minimum_tetradius<<","<<maximum_tetradius<<","<<average_tetradius<<","<<maximum_volume<<","<<minimum_volume<<","<<average_volume<<"\n";
  }
}
//-------------------------------------------------------
// Reconstruct::find tetrahedron V0.0.2
//-------------------------------------------------------
// void Graphmeshcls::Reconstruct_tets(std::vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
// {
//   tetgenio pp_in, tet_out;
//   tetgenio::facet *f;
//   tetgenio::polygon *p;

//   // All indices start from 0
//   pp_in.firstnumber = 0;
//   tet_out.firstnumber = 0;

//   pp_in.numberofpoints = particle_total.size();
//   pp_in.pointlist = new Real[pp_in.numberofpoints * 3];
//   for (int i=0; i < pp_in.numberofpoints; i++){
//     pp_in.pointlist[i * 3 + 0] = particle_total[i]->coord.i;
//     pp_in.pointlist[i * 3 + 1] = particle_total[i]->coord.j;
//     pp_in.pointlist[i * 3 + 2] = particle_total[i]->coord.k;
//   }

//   // tetrahedralize("CnlV", &pp_in, &tet_out);
//   tetrahedralize("CnlVVO/2", &pp_in, &tet_out);

//   tris.clear();
//   tets.clear();
//   tets_error.clear();

//   //load the tets from tetgen the first time
//   int   n_tets       = tet_out.numberoftetrahedra;
//   int  *tet_list     = tet_out.tetrahedronlist;
//   int  *neighborlist = tet_out.neighborlist;
//   if (neighborlist == (int *) NULL) cout<<" !! ERROR neighborlist empty"<<endl;

//   for (int i=0; i < n_tets; i++){

//     int4_graph tet;
//     tet.i      = tet_list[i*4 + 0];
//     tet.j      = tet_list[i*4 + 1];
//     tet.k      = tet_list[i*4 + 2];
//     tet.l      = tet_list[i*4 + 3];
//     tet.vol    = 0.;
//     tet.aspect = 0.;
//     tets.insert (tet);
//   }

//   if (tet_out.numberofpoints != particle_total.size())
//     cout<<"<<<<< ERROR!! number of output points is wrong"<<endl;

//   for (int i=0; i < tet_out.numberofpoints; i++){
//     particle_total[i]->coord.i = tet_out.pointlist[i * 3 + 0];
//     particle_total[i]->coord.j = tet_out.pointlist[i * 3 + 1];
//     particle_total[i]->coord.k = tet_out.pointlist[i * 3 + 2];
//   }

//   std::set<int4_graph>::iterator it;
//   std::vector<std::set<int4_graph>::iterator> tet_delete;
//   tet_delete.clear();

//   for (it=tets.begin(); it!=tets.end(); ++it){
//     int p1 = int((*it).i);
//     int p2 = int((*it).j);
//     int p3 = int((*it).k);
//     int p4 = int((*it).l);

//     p_Particle particle_1 = particle_total[p1];
//     p_Particle particle_2 = particle_total[p2];
//     p_Particle particle_3 = particle_total[p3];
//     p_Particle particle_4 = particle_total[p4];

//     if( particle_1->type != REAL_PARTICLE && 
//         particle_2->type != REAL_PARTICLE && 
//         particle_3->type != REAL_PARTICLE &&
//         particle_4->type != REAL_PARTICLE)
//     {
//       tet_delete.push_back(it);
//       tets_error.insert(*(it));
//     }else if(particle_1->type != REAL_PARTICLE || 
//              particle_2->type != REAL_PARTICLE || 
//              particle_3->type != REAL_PARTICLE ||
//              particle_4->type != REAL_PARTICLE){
//         // judge something
//         // construct particle temporally
//       Particle particle_temp;
//       particle_temp.coord.i = 1./4.*(particle_1->coord.i + particle_2->coord.i + particle_3->coord.i + particle_4->coord.i);
//       particle_temp.coord.j = 1./4.*(particle_1->coord.j + particle_2->coord.j + particle_3->coord.j + particle_4->coord.j);
//       particle_temp.coord.k = 1./4.*(particle_1->coord.k + particle_2->coord.k + particle_3->coord.k + particle_4->coord.k);
//       particle_temp.Calculate_particle_infor(sph);

//       Real phi_size         = 1./4.*(particle_1->phi + particle_2->phi + particle_3->phi + particle_4->phi);

//       if(particle_temp.phi <= 0.0001*phi_size){
//         tet_delete.push_back(it);
//         tets_error.insert(*(it));
//       }
//     }
//   }

//   if (tet_delete.size()>0){
//     int n_delete = tet_delete.size();
//     cout<<"<<<<< Num. of tets to be deleted: "<<tet_delete.size()<<endl;

//     for (int i=0; i < n_delete; i++){
//       tets.erase(tet_delete[i]);
//     }

//     if (tets.size() + n_delete != tet_out.numberoftetrahedra){
//       cout<<"<<<<< ERROR: number of tets are wrong!!!"<<endl;
//     }
//   }

//   Tet_quality_statistics (particle_total, sph, world);

// #ifdef _TEST_
//   char    filename[256];
//   sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
//   ofstream out(filename, ios::app);
//   out<<"<<<<<Reconstruct_tets finished\n";
//   out.close();
// #endif
// }
//-------------------------------------------------------
// Reconstruct::find tetrahedron V0.0.1
//-------------------------------------------------------
// void Graphmeshcls::Reconstruct_tets(std::vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
// {
//   tets.clear();

//   std::vector <p_Particle> error_particle;

//   std::pair<EdgeIterator, EdgeIterator> ei = edges(graph);

//   // iterate all the triangles to find the tets
//   std::set<int3_graph>::iterator it;
//   for (it=tris.begin(); it!=tris.end(); ++it){
//     VertexDescriptor p1 = VertexDescriptor((*it).i);
//     VertexDescriptor p2 = VertexDescriptor((*it).j);
//     VertexDescriptor p3 = VertexDescriptor((*it).k);

//     // cout<<"<<<< in face: "<<(*it).i<<" "<<(*it).j<<" "<<(*it).k<<endl;

//     std::vector<VertexDescriptor> vertex_record;
//     vertex_record.clear();

//     AdjacencyIterator adjvert_p1_i, adjvert_p1_i_end, p1_next;
//     boost::tie( adjvert_p1_i, adjvert_p1_i_end) = adjacent_vertices(p1, graph);
//     for (p1_next = adjvert_p1_i; adjvert_p1_i != adjvert_p1_i_end; adjvert_p1_i = p1_next) {
//       // cout<<"<<<< p1_next "<<*p1_next<<endl;
//       if (*p1_next != p2 && *p1_next != p3){

//         AdjacencyIterator adjvert_p2_i, adjvert_p2_i_end, p2_next;
//         boost::tie( adjvert_p2_i, adjvert_p2_i_end) = adjacent_vertices(p2, graph);
//         for (p2_next = adjvert_p2_i; adjvert_p2_i != adjvert_p2_i_end; adjvert_p2_i = p2_next) {
//           // cout<<"<<<< p2_next "<<*p2_next<<endl;
//           if (*p2_next != p1 && *p2_next != p3){

//             AdjacencyIterator adjvert_p3_i, adjvert_p3_i_end, p3_next;
//             boost::tie( adjvert_p3_i, adjvert_p3_i_end) = adjacent_vertices(p3, graph);
//             for (p3_next = adjvert_p3_i; adjvert_p3_i != adjvert_p3_i_end; adjvert_p3_i = p3_next) {
//               // cout<<"<<<< p3_next "<<*p3_next<<endl;
//               if (*p3_next != p1 && *p3_next != p2){

//                 if (*p1_next==*p2_next && *p1_next==*p3_next){
//                   int id_fourth = int (*p1_next);
//                   vertex_record.push_back (*p1_next);

//                   int4_graph tet;

//                   std::vector<int> vlist; 
//                   vlist.push_back((*it).i);
//                   vlist.push_back((*it).j);
//                   vlist.push_back((*it).k);
//                   vlist.push_back(id_fourth);
//                   std::sort (vlist.begin(), vlist.end(), Greater());

//                   tet.i = vlist[0];
//                   tet.j = vlist[1];
//                   tet.k = vlist[2];
//                   tet.l = vlist[3];

//                   tets.insert (tet);

//                   // cout<<"<<<< tet generated: "<<tet.i<<" "<<tet.j<<" "<<tet.k<<" "<<tet.l<<endl;

//                   if(int(vertex_record.size()) > 2)
//                   {
//                     cout<<"Something is wrong about the tet created !!! "<<endl;

//                     error_particle.push_back(particle_total[tet.i]);
//                     error_particle.push_back(particle_total[tet.j]);
//                     error_particle.push_back(particle_total[tet.k]);
//                     error_particle.push_back(particle_total[tet.l]);

//                     //exit(0);
//                   }else if(int(vertex_record.size())== 2)
//                   {
//                     std::pair<EdgeDescriptor, bool> test_1 = edge(graph[vertex_record[0]],graph[vertex_record[1]],graph);
//                     std::pair<EdgeDescriptor, bool> test_2 = edge(graph[vertex_record[1]],graph[vertex_record[0]],graph);
//                     if (test_1.second == true || test_2.second == true) {
//                       cout<<"Some edges in the graph intersect !!! "<<endl;

//                       error_particle.push_back(particle_total[tet.i]);
//                       error_particle.push_back(particle_total[tet.j]);
//                       error_particle.push_back(particle_total[tet.k]);
//                       error_particle.push_back(particle_total[tet.l]);

//                         //exit(0);
//                     }
//                   }
//                 }
//               }
//               p3_next++;
//             }
//           }
//           p2_next++;
//         }
//       }
//       p1_next++;
//     }
//   }

//   if (error_particle.size() > 0){
//     cout<<"<<<<< Num. of error particle for tets: "<<error_particle.size()<<endl;
//   }
// #ifdef _TEST_
//   char    filename[256];
//   sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
//   ofstream out(filename, ios::app);
//   out<<"<<<<<Reconstruct_tets finished\n";
//   out.close();
// #endif
// }
#endif
