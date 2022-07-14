#include "mesh_binary_writer.h"
#include "sph_mesh_generation.h"
#include "particle_mesh_generation.h"
#include "graph_mesh.h"

/***************************************************/
/*                                                 */
/* Functions defined in class "Mesh_binary_writer" */
/*                                                 */
/***************************************************/

//-------------------------------------------------------
// construction
//-------------------------------------------------------
Mesh_binary_writer::Mesh_binary_writer()
{
  print_info["id"           ] = false;
  print_info["local_id"     ] = false;
  print_info["type"         ] = true;
  print_info["phi"          ] = false;
  print_info["curv"         ] = false;
  print_info["scale"        ] = false;
  print_info["h"            ] = true;
  print_info["NN"           ] = false;
  print_info["color"        ] = true;
  print_info["norm"         ] = false;
  print_info["vol"          ] = false;
  print_info["aspect"       ] = false;
  print_info["radius_ratio" ] = true;
  print_info["mindihedangle"] = true;
}

//-------------------------------------------------------
// Output_mesh
//-------------------------------------------------------
void Mesh_binary_writer::Output_mesh(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  // tris
  world.barrier();
  
  if (DIM == 2){
    Output_vtu_file_tris (n_out, particle_total, sph, world);
    
    world.barrier();
    if (world.rank() == 0){
      Output_pvtu_file_tris (n_out, particle_total, sph, world);
      Output_pvd_file_tris (n_out, particle_total, sph, world);
    }
    world.barrier();
  }
  // tets
  world.barrier();

  if (DIM == 3){
    Output_vtu_file_tets (n_out, particle_total, sph, world);
    
    world.barrier();
    if (world.rank() == 0){
      Output_pvtu_file_tets (n_out, particle_total, sph, world);
      Output_pvd_file_tets (n_out, particle_total, sph, world);
    }
    world.barrier();
  }
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Output_mesh finished\n";
  out.close();
#endif
}

//-------------------------------------------------------
// Output_pvd_file
//-------------------------------------------------------
void Mesh_binary_writer::Output_pvd_file_tris(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s","./outdata/mesh_tris.pvd");
  std::string line;
  std::vector<std::string> lines;

  if (sph->n_post == 0){
    ofstream pvdfile(filename, ios::trunc);

    pvdfile<<"<VTKFile type=\"Collection\" version=\"0.1\">"<<endl;
    pvdfile<<"  <Collection>"<<endl;
    pvdfile<<"    <DataSet timestep=\""<<sph->run_time<<"\" group=\"\" part=\"0\" file=\"./mesh_tris_step_"<<n_out<<".pvtu\"/>"<<endl;
    pvdfile<<"  </Collection>"<<endl;
    pvdfile<<"</VTKFile>"<<endl;
    pvdfile.close();
  }else{
    ifstream pvdfile(filename);
    while (std::getline(pvdfile, line)){
      lines.push_back(line);
    }
    pvdfile.close();

    ofstream pvdfile_new(filename, ios::trunc);

    for (int i=0; i<lines.size()-2; i++){
      pvdfile_new<<lines.at(i)<<endl;
    }
    pvdfile_new<<"  <DataSet timestep=\""<<sph->run_time<<"\" group=\"\" part=\"0\" file=\"./mesh_tris_step_"<<n_out<<".pvtu\"/>"<<endl;
	pvdfile_new<<"  </Collection>"<<endl;
	pvdfile_new<<"</VTKFile>"<<endl;
	pvdfile_new.close();
  }
}
//-------------------------------------------------------
// Output_pvd_file
//-------------------------------------------------------
void Mesh_binary_writer::Output_pvd_file_tets(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s","./outdata/mesh_tets.pvd");
  std::string line;
  std::vector<std::string> lines;

  if (sph->n_post == 0){
    ofstream pvdfile(filename, ios::trunc);

    pvdfile<<"<VTKFile type=\"Collection\" version=\"0.1\">"<<endl;
    pvdfile<<"  <Collection>"<<endl;
    pvdfile<<"    <DataSet timestep=\""<<sph->run_time<<"\" group=\"\" part=\"0\" file=\"./mesh_tets_step_"<<n_out<<".pvtu\"/>"<<endl;
    pvdfile<<"  </Collection>"<<endl;
    pvdfile<<"</VTKFile>"<<endl;
    pvdfile.close();
  }else{
    ifstream pvdfile(filename);
    while (std::getline(pvdfile, line)){
      lines.push_back(line);
    }
    pvdfile.close();

    ofstream pvdfile_new(filename, ios::trunc);

    for (int i=0; i<lines.size()-2; i++){
      pvdfile_new<<lines.at(i)<<endl;
    }
    pvdfile_new<<"  <DataSet timestep=\""<<sph->run_time<<"\" group=\"\" part=\"0\" file=\"./mesh_tets_step_"<<n_out<<".pvtu\"/>"<<endl;
  pvdfile_new<<"  </Collection>"<<endl;
  pvdfile_new<<"</VTKFile>"<<endl;
  pvdfile_new.close();
  }
}
//-------------------------------------------------------
//Output_vtu_file
//-------------------------------------------------------
void Mesh_binary_writer::Output_vtu_file_tris(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s%d%s%d%s","./outdata/mesh_tris_rank_",world.rank(),"_step_",n_out,".vtu");

  int np    = particle_total.size();
  int ncell = sph->mesh.tris.size();

  ofstream     vtufile(filename, ios::binary);
  stringstream dataStream;
  int          dataOffset = 0;
  unsigned int num_bytes  = 0;

  vtufile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  vtufile<<" <UnstructuredGrid>"<<endl;
  vtufile<<"  <Piece NumberOfPoints=\""<<np<<"\" NumberOfCells=\""<<ncell<<"\">"<<endl;
  vtufile<<"    <PointData Vectors=\"vector\">"<<endl;
  
  if (print_info["id"]){
    vtufile<<"      <DataArray type=\"Int32\" Name=\"id\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          int var = int(particle_total[i]->id);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }
  
  if (print_info["local_id"]){
    vtufile<<"      <DataArray type=\"Int32\" Name=\"local_id\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          int var = int(particle_total[i]->local_id);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }
  
  if (print_info["type"]){
    vtufile<<"      <DataArray type=\"Int32\" Name=\"type\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          int var = int(particle_total[i]->type);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }
  
  if (print_info["phi"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"phi\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var = float(particle_total[i]->phi);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }

  if (print_info["h"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"h\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var = float(particle_total[i]->h);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }

  if (print_info["scale"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"scale\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var = float(particle_total[i]->scale);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }

  if (print_info["curv"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"curv\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var = float(particle_total[i]->curv);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }

  if (print_info["NN"]){
    vtufile<<"      <DataArray type=\"Int32\" Name=\"NN\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          int var = int(particle_total[i]->neighbor.size());
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }
  
  if (print_info["color"]){
    vtufile<<"      <DataArray type=\"Int32\" Name=\"color\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          int var = int(particle_total[i]->color);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }
  
  if (print_info["norm"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"norm\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*3*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var1 = float(particle_total[i]->norm.i);
          float var2 = float(particle_total[i]->norm.j);
          float var3 = float(particle_total[i]->norm.k);
          dataStream.write(reinterpret_cast<char*>(&var1), sizeof(var1));
          dataStream.write(reinterpret_cast<char*>(&var2), sizeof(var2));
          dataStream.write(reinterpret_cast<char*>(&var3), sizeof(var3));
    }
  }
  
  vtufile<<"    </PointData>"<<endl;
  vtufile<<"    <CellData Vectors=\"vector\">"<<endl;
  if (print_info["vol"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"vol\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(ncell*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    std::set<int3_graph>::iterator it;
    for (it=sph->mesh.tris.begin(); it!=sph->mesh.tris.end(); ++it){
      float vol = float((*it).vol);
      dataStream.write(reinterpret_cast<char*>(&vol), sizeof(vol));
    }
  }
  if (print_info["aspect"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"aspect\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(ncell*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    std::set<int3_graph>::iterator it;
    for (it=sph->mesh.tris.begin(); it!=sph->mesh.tris.end(); ++it){
      float aspect = float((*it).aspect);
      dataStream.write(reinterpret_cast<char*>(&aspect), sizeof(aspect));
    }
  }
  vtufile<<"    </CellData>"<<endl;
  vtufile<<"    <Points>"<<endl;
  vtufile<<"      <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
  num_bytes = (unsigned int)(np*3*sizeof(float));
  dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
  dataOffset += num_bytes + sizeof(num_bytes);
  
  for (int i=0; i<np; i++){
  	float coord_i = float(particle_total[i]->coord.i);
  	float coord_j = float(particle_total[i]->coord.j);
  	float coord_k = float(particle_total[i]->coord.k);

  	dataStream.write(reinterpret_cast<char*>(&coord_i), sizeof(coord_i));
  	dataStream.write(reinterpret_cast<char*>(&coord_j), sizeof(coord_j));
  	dataStream.write(reinterpret_cast<char*>(&coord_k), sizeof(coord_k));
  }
//   vtufile<<"      </DataArray>"<<endl;
  vtufile<<"    </Points>"<<endl;
  vtufile<<"    <Cells>"<<endl;
  vtufile<<"      <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
  num_bytes = (unsigned int)(ncell*3*sizeof(int));
  dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
  dataOffset += num_bytes + sizeof(num_bytes);
  
  std::set<int3_graph>::iterator it;
  for (it=sph->mesh.tris.begin(); it!=sph->mesh.tris.end(); ++it){
  	int p1 = int((*it).i);
  	int p2 = int((*it).j);
  	int p3 = int((*it).k);
  	dataStream.write(reinterpret_cast<char*>(&p1), sizeof(p1));
  	dataStream.write(reinterpret_cast<char*>(&p2), sizeof(p2));
  	dataStream.write(reinterpret_cast<char*>(&p3), sizeof(p3));
  }
//   vtufile<<"      </DataArray>"<<endl;
  vtufile<<"      <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
  num_bytes = (unsigned int)(ncell*1*sizeof(int));
  dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
  dataOffset += num_bytes + sizeof(num_bytes);
  
  int offset = 3;
  for (int i=0; i<ncell; i++){
        int _offset=offset;
  	dataStream.write(reinterpret_cast<char*>(&_offset), sizeof(_offset));
  	offset += 3;
  }
//   vtufile<<"      </DataArray>"<<endl;
  vtufile<<"      <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
  num_bytes = (unsigned int)(ncell*1*sizeof(uint8_t));
  dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
  dataOffset += num_bytes + sizeof(num_bytes);
  
  for (int i=0; i<ncell; i++){
        uint8_t type = 5;
  	dataStream.write(reinterpret_cast<char*>(&type), sizeof(type));
  }
//   vtufile<<"      </DataArray>"<<endl;
  vtufile<<"    </Cells>"<<endl;
  vtufile<<"  </Piece>"<<endl;
  vtufile<<" </UnstructuredGrid>"<<endl;
  vtufile<<" <AppendedData encoding=\"raw\">"<<endl;
  vtufile<<" _"<<dataStream.rdbuf()<<endl;
  vtufile<<" </AppendedData>"<<endl;
  vtufile<<"</VTKFile>"<<endl;

  vtufile.flush();
  vtufile.close();
}
//-------------------------------------------------------
//Output_vtu_file
//-------------------------------------------------------
void Mesh_binary_writer::Output_vtu_file_tets(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s%d%s%d%s","./outdata/mesh_tets_rank_",world.rank(),"_step_",n_out,".vtu");

  int np    = particle_total.size();
  int ncell = sph->mesh.tets.size();

  ofstream     vtufile(filename, ios::binary);
  stringstream dataStream;
  int          dataOffset = 0;
  unsigned int num_bytes  = 0;

  vtufile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  vtufile<<" <UnstructuredGrid>"<<endl;
  vtufile<<"  <Piece NumberOfPoints=\""<<np<<"\" NumberOfCells=\""<<ncell<<"\">"<<endl;
  vtufile<<"    <PointData Vectors=\"vector\">"<<endl;
  
  if (print_info["id"]){
    vtufile<<"      <DataArray type=\"Int32\" Name=\"id\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          int var = int(particle_total[i]->id);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }
  
  if (print_info["local_id"]){
    vtufile<<"      <DataArray type=\"Int32\" Name=\"local_id\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          int var = int(particle_total[i]->local_id);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }
  
  if (print_info["type"]){
    vtufile<<"      <DataArray type=\"Int32\" Name=\"type\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          int var = int(particle_total[i]->type);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }
  
  if (print_info["phi"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"phi\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var = float(particle_total[i]->phi);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }

  if (print_info["h"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"h\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var = float(particle_total[i]->h);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }

  if (print_info["scale"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"scale\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var = float(particle_total[i]->scale);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }

  if (print_info["curv"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"curv\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var = float(particle_total[i]->curv);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }

  if (print_info["NN"]){
    vtufile<<"      <DataArray type=\"Int32\" Name=\"NN\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          int var = int(particle_total[i]->neighbor.size());
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }
  
  if (print_info["color"]){
    vtufile<<"      <DataArray type=\"Int32\" Name=\"color\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          int var = int(particle_total[i]->color);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }
  
  if (print_info["norm"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"norm\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*3*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var1 = float(particle_total[i]->norm.i);
          float var2 = float(particle_total[i]->norm.j);
          float var3 = float(particle_total[i]->norm.k);
          dataStream.write(reinterpret_cast<char*>(&var1), sizeof(var1));
          dataStream.write(reinterpret_cast<char*>(&var2), sizeof(var2));
          dataStream.write(reinterpret_cast<char*>(&var3), sizeof(var3));
    }
  }
  vtufile<<"    </PointData>"<<endl;
  vtufile<<"    <CellData Vectors=\"vector\">"<<endl;
  if (print_info["vol"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"vol\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(ncell*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    std::set<int4_graph>::iterator it;
    for (it=sph->mesh.tets.begin(); it!=sph->mesh.tets.end(); ++it){
      float vol = float((*it).vol);
      dataStream.write(reinterpret_cast<char*>(&vol), sizeof(vol));
    }
  }
  if (print_info["aspect"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"aspect\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(ncell*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    std::set<int4_graph>::iterator it;
    for (it=sph->mesh.tets.begin(); it!=sph->mesh.tets.end(); ++it){
      float aspect = float((*it).aspect);
      dataStream.write(reinterpret_cast<char*>(&aspect), sizeof(aspect));
    }
  }
  if (print_info["radius_ratio"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"radius_ratio\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(ncell*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    std::set<int4_graph>::iterator it;
    for (it=sph->mesh.tets.begin(); it!=sph->mesh.tets.end(); ++it){
      float radius_ratio = float((*it).radius_ratio);
      dataStream.write(reinterpret_cast<char*>(&radius_ratio), sizeof(radius_ratio));
    }
  }
  if (print_info["mindihedangle"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"mindihedangle\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(ncell*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    std::set<int4_graph>::iterator it;
    for (it=sph->mesh.tets.begin(); it!=sph->mesh.tets.end(); ++it){
      float mindihedangle = float((*it).mindihedangle);
      dataStream.write(reinterpret_cast<char*>(&mindihedangle), sizeof(mindihedangle));
    }
  }
  vtufile<<"    </CellData>"<<endl;
  vtufile<<"    <Points>"<<endl;
  vtufile<<"      <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
  num_bytes = (unsigned int)(np*3*sizeof(float));
  dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
  dataOffset += num_bytes + sizeof(num_bytes);
  
  for (int i=0; i<np; i++){
    float coord_i = float(particle_total[i]->coord.i);
    float coord_j = float(particle_total[i]->coord.j);
    float coord_k = float(particle_total[i]->coord.k);

    dataStream.write(reinterpret_cast<char*>(&coord_i), sizeof(coord_i));
    dataStream.write(reinterpret_cast<char*>(&coord_j), sizeof(coord_j));
    dataStream.write(reinterpret_cast<char*>(&coord_k), sizeof(coord_k));
  }
//   vtufile<<"      </DataArray>"<<endl;
  vtufile<<"    </Points>"<<endl;
  vtufile<<"    <Cells>"<<endl;
  vtufile<<"      <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
  num_bytes = (unsigned int)(ncell*4*sizeof(int));
  dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
  dataOffset += num_bytes + sizeof(num_bytes);
  
  std::set<int4_graph>::iterator it;
  for (it=sph->mesh.tets.begin(); it!=sph->mesh.tets.end(); ++it){
    int p1 = int((*it).i);
    int p2 = int((*it).j);
    int p3 = int((*it).k);
    int p4 = int((*it).l);
    dataStream.write(reinterpret_cast<char*>(&p1), sizeof(p1));
    dataStream.write(reinterpret_cast<char*>(&p2), sizeof(p2));
    dataStream.write(reinterpret_cast<char*>(&p3), sizeof(p3));
    dataStream.write(reinterpret_cast<char*>(&p4), sizeof(p4));
  }
//   vtufile<<"      </DataArray>"<<endl;
  vtufile<<"      <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
  num_bytes = (unsigned int)(ncell*1*sizeof(int));
  dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
  dataOffset += num_bytes + sizeof(num_bytes);
  
  int offset = 4;
  for (int i=0; i<ncell; i++){
        int _offset=offset;
    dataStream.write(reinterpret_cast<char*>(&_offset), sizeof(_offset));
    offset += 4;
  }
//   vtufile<<"      </DataArray>"<<endl;
  vtufile<<"      <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
  num_bytes = (unsigned int)(ncell*1*sizeof(uint8_t));
  dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
  dataOffset += num_bytes + sizeof(num_bytes);
  
  for (int i=0; i<ncell; i++){
        uint8_t type = 10;
    dataStream.write(reinterpret_cast<char*>(&type), sizeof(type));
  }
//   vtufile<<"      </DataArray>"<<endl;
  vtufile<<"    </Cells>"<<endl;
  vtufile<<"  </Piece>"<<endl;
  vtufile<<" </UnstructuredGrid>"<<endl;
  vtufile<<" <AppendedData encoding=\"raw\">"<<endl;
  vtufile<<" _"<<dataStream.rdbuf()<<endl;
  vtufile<<" </AppendedData>"<<endl;
  vtufile<<"</VTKFile>"<<endl;

  vtufile.flush();
  vtufile.close();
}
// void Mesh_binary_writer::Output_vtu_file(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
// {
//   FILE    *fp;
//   char    filename[256];
//   sprintf(filename,"%s%d%s%d%s","./outdata/mesh_rank_",world.rank(),"_step_",n_out,".vtu");
// 
//   int np    = particle_total.size();
//   int ncell = sph->mesh.tris.size();
// 
//   ofstream     vtufile(filename, ios::trunc);
//   stringstream dataStream;
// 
//   vtufile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
//   vtufile<<" <UnstructuredGrid>"<<endl;
//   vtufile<<"  <Piece NumberOfPoints=\""<<np<<"\" NumberOfCells=\""<<ncell<<"\">"<<endl;
//   vtufile<<"    <PointData Vectors=\"vector\">"<<endl;
//   vtufile<<"    </PointData>"<<endl;
//   vtufile<<"    <Points>"<<endl;
//   vtufile<<"      <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
//   for (int i=0; i<np; i++){
//         float coord_i = float(particle_total[i]->coord.i);
//         float coord_j = float(particle_total[i]->coord.j);
//         float coord_k = float(particle_total[i]->coord.k);
// 
//         vtufile<<coord_i<<" "<<coord_j<<" "<<coord_k<<endl;
//   }
//   vtufile<<"      </DataArray>"<<endl;
//   vtufile<<"    </Points>"<<endl;
//   vtufile<<"    <Cells>"<<endl;
//   vtufile<<"      <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<endl;
//   std::set<int3_graph>::iterator it;
//   for (it=sph->mesh.tris.begin(); it!=sph->mesh.tris.end(); ++it){
//         int p1 = int((*it).i);
//         int p2 = int((*it).j);
//         int p3 = int((*it).k);
//         
//         vtufile<<p1<<" "<<p2<<" "<<p3<<endl;
//   }
//   vtufile<<"      </DataArray>"<<endl;
//   vtufile<<"      <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<endl;
//   int offset = 3;
//   for (int i=0; i<ncell; i++){
//         vtufile<<offset<<endl;
//         offset += 3;
//   }
//   vtufile<<"      </DataArray>"<<endl;
//   vtufile<<"      <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<endl;
//   uint8_t type = 5;
//   for (int i=0; i<ncell; i++){
//         vtufile<<5<<endl;
//   }
//   vtufile<<"      </DataArray>"<<endl;
//   vtufile<<"    </Cells>"<<endl;
//   vtufile<<"  </Piece>"<<endl;
//   vtufile<<" </UnstructuredGrid>"<<endl;
//   vtufile<<"</VTKFile>"<<endl;
// 
//   vtufile.close();
// }
//-------------------------------------------------------
// Output_pvtu_file
//-------------------------------------------------------
void Mesh_binary_writer::Output_pvtu_file_tris(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s%d%s","./outdata/mesh_tris_step_",n_out,".pvtu");

  int np    = particle_total.size();
  int ncell = sph->mesh.tris.size();

  ofstream     pvtufile(filename, ios::trunc);

  pvtufile<<"<?xml version=\"1.0\"?>"<<endl;
  pvtufile<<"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  pvtufile<<" <PUnstructuredGrid GhostLevel=\"0\">"<<endl;
  pvtufile<<"  <PPointData Vectors=\"vector\">"<<endl;
  if (print_info["id"      ]) pvtufile<<"    <PDataArray type=\"Int32\" Name=\"id\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["local_id"]) pvtufile<<"    <PDataArray type=\"Int32\" Name=\"local_id\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["type"    ]) pvtufile<<"    <PDataArray type=\"Int32\" Name=\"type\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["phi"     ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"phi\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["h"       ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"h\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["scale"   ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"scale\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["curv"    ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"curv\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["NN"      ]) pvtufile<<"    <PDataArray type=\"Int32\" Name=\"NN\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["color"   ]) pvtufile<<"    <PDataArray type=\"Int32\" Name=\"color\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["norm"    ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"norm\" NumberOfComponents=\"3\"/>"<<endl;
  pvtufile<<"  </PPointData>"<<endl;
  pvtufile<<"  <PCellData Vectors=\"vector\">"<<endl;
  if (print_info["vol"    ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"vol\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["aspect" ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"aspect\" NumberOfComponents=\"1\"/>"<<endl;
  pvtufile<<"  </PCellData>"<<endl;
  pvtufile<<"  <PPoints>"<<endl;
  pvtufile<<"    <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>"<<endl;
  pvtufile<<"  </PPoints>"<<endl;
  for (int i=0; i < world.size(); i++){
    pvtufile<<"    <Piece Source=\"./mesh_tris_rank_"<<i<<"_step_"<<n_out<<".vtu\"/>"<<endl;
  }
  pvtufile<<" </PUnstructuredGrid>"<<endl;
  pvtufile<<"</VTKFile>"<<endl;

  pvtufile.close();
}
//-------------------------------------------------------
// Output_pvtu_file
//-------------------------------------------------------
void Mesh_binary_writer::Output_pvtu_file_tets(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s%d%s","./outdata/mesh_tets_step_",n_out,".pvtu");

  int np    = particle_total.size();
  int ncell = sph->mesh.tris.size();

  ofstream     pvtufile(filename, ios::trunc);

  pvtufile<<"<?xml version=\"1.0\"?>"<<endl;
  pvtufile<<"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  pvtufile<<" <PUnstructuredGrid GhostLevel=\"0\">"<<endl;
  pvtufile<<"  <PPointData Vectors=\"vector\">"<<endl;
  if (print_info["id"           ]) pvtufile<<"    <PDataArray type=\"Int32\" Name=\"id\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["local_id"     ]) pvtufile<<"    <PDataArray type=\"Int32\" Name=\"local_id\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["type"         ]) pvtufile<<"    <PDataArray type=\"Int32\" Name=\"type\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["phi"          ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"phi\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["h"            ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"h\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["scale"        ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"scale\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["curv"         ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"curv\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["NN"           ]) pvtufile<<"    <PDataArray type=\"Int32\" Name=\"NN\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["color"        ]) pvtufile<<"    <PDataArray type=\"Int32\" Name=\"color\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["norm"         ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"norm\" NumberOfComponents=\"3\"/>"<<endl;
  pvtufile<<"  </PPointData>"<<endl;
  pvtufile<<"  <PCellData Vectors=\"vector\">"<<endl;
  if (print_info["vol"          ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"vol\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["aspect"       ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"aspect\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["radius_ratio" ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"radius_ratio\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["mindihedangle"]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"mindihedangle\" NumberOfComponents=\"1\"/>"<<endl;
  pvtufile<<"  </PCellData>"<<endl;
  pvtufile<<"  <PPoints>"<<endl;
  pvtufile<<"    <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>"<<endl;
  pvtufile<<"  </PPoints>"<<endl;
  for (int i=0; i < world.size(); i++){
    pvtufile<<"    <Piece Source=\"./mesh_tets_rank_"<<i<<"_step_"<<n_out<<".vtu\"/>"<<endl;
  }
  pvtufile<<" </PUnstructuredGrid>"<<endl;
  pvtufile<<"</VTKFile>"<<endl;

  pvtufile.close();
}
