#include "particle_binary_writer.h"
#include "sph_mesh_generation.h"
#include "particle_mesh_generation.h"

/***************************************************/
/*                                                 */
/* Functions defined in class "Particle_binary_writer" */
/*                                                 */
/***************************************************/

//-------------------------------------------------------
// construction
//-------------------------------------------------------
Particle_binary_writer::Particle_binary_writer()
{
  print_info["id"]       = true;
  print_info["local_id"] = true;
  print_info["type"]     = true;
  print_info["phi"]      = false;
  print_info["curv"]     = false;
  print_info["shaperd"]  = false;
  print_info["timestep"] = false;
  print_info["scale"]    = false;
  print_info["h"]        = true;
  print_info["NN"]       = true;
  print_info["color"]    = true;
  print_info["norm"]     = true;
  print_info["accel"]    = false;
  print_info["vel"]      = false;
}

//-------------------------------------------------------
// Output_particle
//-------------------------------------------------------
void Particle_binary_writer::Output_particle(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  Output_vtu_file (n_out, particle_total, sph, world);

  world.barrier();

  if (world.rank() == 0){
    Output_pvtu_file (n_out, particle_total, sph, world);
    Output_pvd_file (n_out, particle_total, sph, world);
  }
  
  world.barrier();

#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Output_particle finished\n";
  out.close();
#endif
}

//-------------------------------------------------------
// Output_pvd_file
//-------------------------------------------------------
void Particle_binary_writer::Output_pvd_file(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s","./outdata/particle.pvd");
  std::string line;
  std::vector<std::string> lines;

  if (n_out == 0){
    ofstream pvdfile(filename, ios::trunc);

    pvdfile<<"<VTKFile type=\"Collection\" version=\"0.1\">"<<endl;
    pvdfile<<"  <Collection>"<<endl;
    pvdfile<<"    <DataSet timestep=\""<<sph->run_time<<"\" group=\"\" part=\"0\" file=\"./particle_step_"<<n_out<<".pvtu\"/>"<<endl;
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
    pvdfile_new<<"  <DataSet timestep=\""<<sph->run_time<<"\" group=\"\" part=\"0\" file=\"./particle_step_"<<n_out<<".pvtu\"/>"<<endl;
	pvdfile_new<<"  </Collection>"<<endl;
	pvdfile_new<<"</VTKFile>"<<endl;
	pvdfile_new.close();
  }
}
//-------------------------------------------------------
//Output_vtu_file
//-------------------------------------------------------
void Particle_binary_writer::Output_vtu_file(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s%d%s%d%s","./outdata/particle_rank_",world.rank(),"_step_",n_out,".vtu");

  int np    = particle_total.size();

  ofstream     vtufile(filename, ios::binary);
  stringstream dataStream;
  int          dataOffset = 0;
  unsigned int num_bytes  = 0;

  vtufile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  vtufile<<" <UnstructuredGrid>"<<endl;
  vtufile<<"  <Piece NumberOfPoints=\""<<np<<"\" NumberOfCells=\""<<np<<"\">"<<endl;
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
  
  if (print_info["accel"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"accel\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*3*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var1 = DIM_X ? float(particle_total[i]->a.i) : 0.;
          float var2 = DIM_Y ? float(particle_total[i]->a.j) : 0.;
          float var3 = DIM_Z ? float(particle_total[i]->a.k) : 0.;
          dataStream.write(reinterpret_cast<char*>(&var1), sizeof(var1));
          dataStream.write(reinterpret_cast<char*>(&var2), sizeof(var2));
          dataStream.write(reinterpret_cast<char*>(&var3), sizeof(var3));
    }
  }
  
  if (print_info["vel"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"vel\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*3*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var1 = DIM_X ? float(particle_total[i]->v.i) : 0.;
          float var2 = DIM_Y ? float(particle_total[i]->v.j) : 0.;
          float var3 = DIM_Z ? float(particle_total[i]->v.k) : 0.;
          dataStream.write(reinterpret_cast<char*>(&var1), sizeof(var1));
          dataStream.write(reinterpret_cast<char*>(&var2), sizeof(var2));
          dataStream.write(reinterpret_cast<char*>(&var3), sizeof(var3));
    }
  }

  if (print_info["shaperd"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"shaperd\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var = float(particle_total[i]->sum);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }

  if (print_info["timestep"]){
    vtufile<<"      <DataArray type=\"Float32\" Name=\"timestep\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int i=0; i<np; i++){
          float var = float(particle_total[i]->timestep);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
    }
  }

  vtufile<<"    </PointData>"<<endl;
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
  num_bytes = (unsigned int)(np*1*sizeof(int));
  dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
  dataOffset += num_bytes + sizeof(num_bytes);

  int pp = 0;
  for (int i=0; i<np; i++){
    dataStream.write(reinterpret_cast<char*>(&pp), sizeof(pp));
    pp ++;
  }

//   vtufile<<"      </DataArray>"<<endl;
  vtufile<<"      <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
  num_bytes = (unsigned int)(np*1*sizeof(int));
  dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
  dataOffset += num_bytes + sizeof(num_bytes);
  
  int offset = 1;
  for (int i=0; i<np; i++){
    dataStream.write(reinterpret_cast<char*>(&offset), sizeof(offset));
    offset ++;
  }

//   vtufile<<"      </DataArray>"<<endl;
  vtufile<<"      <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
  num_bytes = (unsigned int)(np*1*sizeof(uint8_t));
  dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
  dataOffset += num_bytes + sizeof(num_bytes);
  
  for (int i=0; i<np; i++){
    uint8_t type = 1;
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
// Output_pvtu_file
//-------------------------------------------------------
void Particle_binary_writer::Output_pvtu_file(int n_out, vector<p_Particle> particle_total, SOLVER *sph, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s%d%s","./outdata/particle_step_",n_out,".pvtu");

  int np    = particle_total.size();

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
  if (print_info["accel"   ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"accel\" NumberOfComponents=\"3\"/>"<<endl;
  if (print_info["vel"     ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"vel\" NumberOfComponents=\"3\"/>"<<endl;
  if (print_info["shaperd" ]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"shaperd\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["timestep"]) pvtufile<<"    <PDataArray type=\"Float32\" Name=\"timestep\" NumberOfComponents=\"1\"/>"<<endl;
  pvtufile<<"  </PPointData>"<<endl;
  pvtufile<<"  <PPoints>"<<endl;
  pvtufile<<"    <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>"<<endl;
  pvtufile<<"  </PPoints>"<<endl;
  for (int i=0; i < world.size(); i++){
    pvtufile<<"    <Piece Source=\"./particle_rank_"<<i<<"_step_"<<n_out<<".vtu\"/>"<<endl;
  }
  pvtufile<<" </PUnstructuredGrid>"<<endl;
  pvtufile<<"</VTKFile>"<<endl;

  pvtufile.close();
}