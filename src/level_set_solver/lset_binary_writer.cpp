#include "lset_binary_writer.h"
#include "level_set.h"
#include "lset_level_infor.h"
#include "lset_package.h"
#include "lset_cell.h"

/********************************************************/
/*                                                      */
/* Functions defined in class "Level_set_binary_writer" */
/*                                                      */
/********************************************************/

//-------------------------------------------------------
// construction
//-------------------------------------------------------
Level_set_binary_writer::Level_set_binary_writer()
{
  print_info["phi"       ] = true;
  print_info["curv"      ] = true;
  print_info["scale"     ] = true;
  print_info["norm"      ] = false;
  print_info["psi"       ] = false;
  print_info["vol"       ] = true;
  print_info["area"      ] = true;
  print_info["interface" ] = true;
  print_info["char"      ] = true;
  #ifdef _MPI_
  print_info["color"     ] = true;
  #endif
}

//-------------------------------------------------------
// Output_mesh
//-------------------------------------------------------
void Level_set_binary_writer::Output_lset(int n_out, p_Levelset_levelinfo lset_level_info, communicator &world)
{
  // tris
  // world.barrier();
  
  Output_vti_file_lset (n_out, lset_level_info, world);
    
  if (world.rank() == 0){
    Output_pvti_file_lset (n_out, lset_level_info, world);
    Output_pvd_file_lset (n_out, lset_level_info, world);
  }
  
  // world.barrier();
  
#ifdef _TEST_
  char    filename[256];
  sprintf(filename,"%s%d%s","./test_log/test_logfile_PID_",world.rank(),".log");
  ofstream out(filename, ios::app);
  out<<"<<<<<Output_lset finished\n";
  out.close();
#endif
}
//-------------------------------------------------------
//Output_vtu_file
//-------------------------------------------------------
void Level_set_binary_writer::Output_vti_file_lset(int n_out, p_Levelset_levelinfo lset_level_info, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s%d%s%d%s","./outdata/mesh_lset_rank_",world.rank(),"_step_",n_out,".vti");

  int         x1 = 0;
  int         y1 = 0;
  int         z1 = 0;
  int         x2 = lset_level_info->num_pkg.i*ICPX;
  int         y2 = lset_level_info->num_pkg.j*ICPY;
  int         z2 = lset_level_info->num_pkg.k*ICPZ;
  int         np = (x2+1-x1)*(y2+1-y1)*(z2+1-z1);
  int      ncell = (x2-x1)*(y2-y1)*(z2-z1);
  Real   spacing = lset_level_info->dl;
  my_real origin = lset_level_info->box_l;
  
  ofstream     vtifile(filename, ios::binary);
  stringstream dataStream;
  int          dataOffset = 0;
  unsigned int num_bytes  = 0;

  vtifile<<"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  vtifile<<" <ImageData WholeExtent=\" "<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<" "<<z1<<" "<<z2<<"\" Origin=\" "<<origin.i<<" "<<origin.j<<" "<<origin.k<<"\" Spacing=\" "<<spacing<<" "<<spacing<<" "<<spacing<<"\">"<<endl;
  vtifile<<"  <Piece Extent=\" "<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<" "<<z1<<" "<<z2<<"\">"<<endl;
  vtifile<<"    <PointData Vectors=\"vector\">"<<endl;
  
  if (print_info["phi"]){
    vtifile<<"      <DataArray type=\"Float32\" Name=\"phi\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);

    for (int k = z1; k <= z2; k++){
      for (int j = y1; j <= y2; j++){
        for (int i = x1; i <= x2; i++){
          int pkg_i = floor(i/ICPX);
          int pkg_j = floor(j/ICPY);
          int pkg_k = floor(k/ICPZ);
          int cell_i = i - pkg_i*ICPX;
          int cell_j = j - pkg_j*ICPY;
          int cell_k = k - pkg_k*ICPZ;
          
          if(k == z2) {pkg_k--; cell_k= DIM_Z ? ICPZ : 0;}
          if(j == y2) {pkg_j--; cell_j= DIM_Y ? ICPY : 0;}
          if(i == x2) {pkg_i--; cell_i= DIM_X ? ICPX : 0;}
                    
          p_Levelset_package current_pkg = lset_level_info->table_lset_pkg_list[pkg_i][pkg_j][pkg_k];

          int  im = DIM_X ? -1 : 0;
          int  jm = DIM_Y ? -1 : 0;
          int  km = DIM_Z ? -1 : 0;

          Real Target      = 0.;
          
          for(int r=im; r<=0; r++)
            for(int s=jm; s<=0; s++)
              for(int t=km; t<=0; t++){
                Target += current_pkg->p_cell[r+cell_i][s+cell_j][t+cell_k]->phi;
              }
              
          Target *= pow(0.5, DIM);

          float var = float(Target);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
        }
      }
    }
  }

  if (print_info["curv"]){
    vtifile<<"      <DataArray type=\"Float32\" Name=\"curv\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int k = z1; k <= z2; k++){
      for (int j = y1; j <= y2; j++){
        for (int i = x1; i <= x2; i++){
          int pkg_i = floor(i/ICPX);
          int pkg_j = floor(j/ICPY);
          int pkg_k = floor(k/ICPZ);
          int cell_i = i - pkg_i*ICPX;
          int cell_j = j - pkg_j*ICPY;
          int cell_k = k - pkg_k*ICPZ;
          
          if(k == z2) {pkg_k--; cell_k= DIM_Z ? ICPZ : 0;}
          if(j == y2) {pkg_j--; cell_j= DIM_Y ? ICPY : 0;}
          if(i == x2) {pkg_i--; cell_i= DIM_X ? ICPX : 0;}
                    
          p_Levelset_package current_pkg = lset_level_info->table_lset_pkg_list[pkg_i][pkg_j][pkg_k];

          int  im = DIM_X ? -1 : 0;
          int  jm = DIM_Y ? -1 : 0;
          int  km = DIM_Z ? -1 : 0;

          Real Target      = 0.;
          
          for(int r=im; r<=0; r++)
            for(int s=jm; s<=0; s++)
              for(int t=km; t<=0; t++){
                Target += current_pkg->p_cell[r+cell_i][s+cell_j][t+cell_k]->curv;
              }
              
          Target *= pow(0.5, DIM);

          float var = float(Target);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
        }
      }
    }
  }
  
  if (print_info["scale"]){
    vtifile<<"      <DataArray type=\"Float32\" Name=\"scale\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int k = z1; k <= z2; k++){
      for (int j = y1; j <= y2; j++){
        for (int i = x1; i <= x2; i++){
          int pkg_i = floor(i/ICPX);
          int pkg_j = floor(j/ICPY);
          int pkg_k = floor(k/ICPZ);
          int cell_i = i - pkg_i*ICPX;
          int cell_j = j - pkg_j*ICPY;
          int cell_k = k - pkg_k*ICPZ;
          
          if(k == z2) {pkg_k--; cell_k= DIM_Z ? ICPZ : 0;}
          if(j == y2) {pkg_j--; cell_j= DIM_Y ? ICPY : 0;}
          if(i == x2) {pkg_i--; cell_i= DIM_X ? ICPX : 0;}
                    
          p_Levelset_package current_pkg = lset_level_info->table_lset_pkg_list[pkg_i][pkg_j][pkg_k];

          int  im = DIM_X ? -1 : 0;
          int  jm = DIM_Y ? -1 : 0;
          int  km = DIM_Z ? -1 : 0;

          Real Target      = 0.;
          
          for(int r=im; r<=0; r++)
            for(int s=jm; s<=0; s++)
              for(int t=km; t<=0; t++){
                Target += current_pkg->p_cell[r+cell_i][s+cell_j][t+cell_k]->scale;
              }
              
          Target *= pow(0.5, DIM);

          float var = float(Target);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
        }
      }
    }
  }
  
  if (print_info["norm"]){
    vtifile<<"      <DataArray type=\"Float32\" Name=\"norm\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*3*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int k = z1; k <= z2; k++){
      for (int j = y1; j <= y2; j++){
        for (int i = x1; i <= x2; i++){
          int pkg_i = floor(i/ICPX);
          int pkg_j = floor(j/ICPY);
          int pkg_k = floor(k/ICPZ);
          int cell_i = i - pkg_i*ICPX;
          int cell_j = j - pkg_j*ICPY;
          int cell_k = k - pkg_k*ICPZ;
          
          if(k == z2) {pkg_k--; cell_k= DIM_Z ? ICPZ : 0;}
          if(j == y2) {pkg_j--; cell_j= DIM_Y ? ICPY : 0;}
          if(i == x2) {pkg_i--; cell_i= DIM_X ? ICPX : 0;}
                    
          p_Levelset_package current_pkg = lset_level_info->table_lset_pkg_list[pkg_i][pkg_j][pkg_k];

          int  im = DIM_X ? -1 : 0;
          int  jm = DIM_Y ? -1 : 0;
          int  km = DIM_Z ? -1 : 0;

          Real Target1      = 0.;
          Real Target2      = 0.;
          Real Target3      = 0.;
          
          for(int r=im; r<=0; r++)
            for(int s=jm; s<=0; s++)
              for(int t=km; t<=0; t++){
                Target1 += current_pkg->p_cell[r+cell_i][s+cell_j][t+cell_k]->n_x;
                Target2 += current_pkg->p_cell[r+cell_i][s+cell_j][t+cell_k]->n_y;
                Target3 += current_pkg->p_cell[r+cell_i][s+cell_j][t+cell_k]->n_z;
              }
              
          Target1 *= pow(0.5, DIM);
          Target2 *= pow(0.5, DIM);
          Target3 *= pow(0.5, DIM);

          float var1 = float(Target1);
          float var2 = float(Target2);
          float var3 = float(Target3);
          dataStream.write(reinterpret_cast<char*>(&var1), sizeof(var1));
          dataStream.write(reinterpret_cast<char*>(&var2), sizeof(var2));
          dataStream.write(reinterpret_cast<char*>(&var3), sizeof(var3));
        }
      }
    }
  }
  
  if (print_info["psi"]){
    vtifile<<"      <DataArray type=\"Float32\" Name=\"psi\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(np*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int k = z1; k <= z2; k++){
      for (int j = y1; j <= y2; j++){
        for (int i = x1; i <= x2; i++){
          int pkg_i = floor(i/ICPX);
          int pkg_j = floor(j/ICPY);
          int pkg_k = floor(k/ICPZ);
          int cell_i = i - pkg_i*ICPX;
          int cell_j = j - pkg_j*ICPY;
          int cell_k = k - pkg_k*ICPZ;
          
          if(k == z2) {pkg_k--; cell_k= DIM_Z ? ICPZ : 0;}
          if(j == y2) {pkg_j--; cell_j= DIM_Y ? ICPY : 0;}
          if(i == x2) {pkg_i--; cell_i= DIM_X ? ICPX : 0;}
                    
          p_Levelset_package current_pkg = lset_level_info->table_lset_pkg_list[pkg_i][pkg_j][pkg_k];

          int  im = DIM_X ? -1 : 0;
          int  jm = DIM_Y ? -1 : 0;
          int  km = DIM_Z ? -1 : 0;

          Real Target      = 0.;
          
          for(int r=im; r<=0; r++)
            for(int s=jm; s<=0; s++)
              for(int t=km; t<=0; t++){
                Target += current_pkg->p_cell[r+cell_i][s+cell_j][t+cell_k]->psi;
              }
              
          Target *= pow(0.5, DIM);

          float var = float(Target);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
        }
      }
    }
  }
  
  vtifile<<"    </PointData>"<<endl;
  vtifile<<"    <CellData Vectors=\"vector\">"<<endl;
  
  if (print_info["interface"]){
    vtifile<<"      <DataArray type=\"Int32\" Name=\"interface\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(ncell*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int k = z1; k < z2; k++){
      for (int j = y1; j < y2; j++){
        for (int i = x1; i < x2; i++){
          int pkg_i = int(i/ICPX);
          int pkg_j = int(j/ICPY);
          int pkg_k = int(k/ICPZ);
          int cell_i = i - pkg_i*ICPX;
          int cell_j = j - pkg_j*ICPY;
          int cell_k = k - pkg_k*ICPZ;
          int var = int(lset_level_info->table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].tag_interface);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
        }
      }
    }
  }
  
  if (print_info["char"]){
    vtifile<<"      <DataArray type=\"Int32\" Name=\"char\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(ncell*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int k = z1; k < z2; k++){
      for (int j = y1; j < y2; j++){
        for (int i = x1; i < x2; i++){
          int pkg_i = int(i/ICPX);
          int pkg_j = int(j/ICPY);
          int pkg_k = int(k/ICPZ);
          int cell_i = i - pkg_i*ICPX;
          int cell_j = j - pkg_j*ICPY;
          int cell_k = k - pkg_k*ICPZ;
          int var = int(lset_level_info->table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].tag_characteristic);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
        }
      }
    }
  }

  #ifdef _MPI_
  if (print_info["color"]){
    vtifile<<"      <DataArray type=\"Int32\" Name=\"color\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(ncell*1*sizeof(int));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int k = z1; k < z2; k++){
      for (int j = y1; j < y2; j++){
        for (int i = x1; i < x2; i++){
          int pkg_i = int(i/ICPX);
          int pkg_j = int(j/ICPY);
          int pkg_k = int(k/ICPZ);
          int var = int(lset_level_info->table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->color);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
        }
      }
    }
  }
  #endif
  
  if (print_info["vol"]){
    vtifile<<"      <DataArray type=\"Float32\" Name=\"vol\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(ncell*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int k = z1; k < z2; k++){
      for (int j = y1; j < y2; j++){
        for (int i = x1; i < x2; i++){
          int pkg_i = int(i/ICPX);
          int pkg_j = int(j/ICPY);
          int pkg_k = int(k/ICPZ);
          int cell_i = i - pkg_i*ICPX;
          int cell_j = j - pkg_j*ICPY;
          int cell_k = k - pkg_k*ICPZ;
          float var = float(lset_level_info->table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].vol);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
        }
      }
    }
  }
  if (print_info["area"]){
    vtifile<<"      <DataArray type=\"Float32\" Name=\"area\" NumberOfComponents=\"1\" format=\"appended\" offset=\""<<dataOffset<<"\"/>"<<endl;
    
    num_bytes = (unsigned int)(ncell*1*sizeof(float));
    dataStream.write(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
    dataOffset += num_bytes + sizeof(num_bytes);
    for (int k = z1; k < z2; k++){
      for (int j = y1; j < y2; j++){
        for (int i = x1; i < x2; i++){
          int pkg_i = int(i/ICPX);
          int pkg_j = int(j/ICPY);
          int pkg_k = int(k/ICPZ);
          int cell_i = i - pkg_i*ICPX;
          int cell_j = j - pkg_j*ICPY;
          int cell_k = k - pkg_k*ICPZ;
          float var = float(lset_level_info->table_lset_pkg_list[pkg_i][pkg_j][pkg_k]->lset_cell[cell_i][cell_j][cell_k].area);
          dataStream.write(reinterpret_cast<char*>(&var), sizeof(var));
        }
      }
    }
  }
  vtifile<<"    </CellData>"<<endl;
  vtifile<<"  </Piece>"<<endl;
  vtifile<<" </ImageData>"<<endl;
  vtifile<<" <AppendedData encoding=\"raw\">"<<endl;
  vtifile<<" _"<<dataStream.rdbuf()<<endl;
  vtifile<<" </AppendedData>"<<endl;
  vtifile<<"</VTKFile>"<<endl;

  vtifile.flush();
  vtifile.close();
}
//-------------------------------------------------------
// Output_pvd_file
//-------------------------------------------------------
void Level_set_binary_writer::Output_pvd_file_lset(int n_out, p_Levelset_levelinfo lset_level_info, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s","./outdata/mesh_lset.pvd");
  std::string line;
  std::vector<std::string> lines;

  if (n_out == 0){
    ofstream pvdfile(filename, ios::trunc);

    pvdfile<<"<VTKFile type=\"Collection\" version=\"0.1\">"<<endl;
    pvdfile<<"  <Collection>"<<endl;
    pvdfile<<"    <DataSet timestep=\""<<n_out<<"\" group=\"\" part=\"0\" file=\"./mesh_lset_step_"<<n_out<<".pvti\"/>"<<endl;
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
    pvdfile_new<<"  <DataSet timestep=\""<<n_out<<"\" group=\"\" part=\"0\" file=\"./mesh_lset_step_"<<n_out<<".pvti\"/>"<<endl;
  pvdfile_new<<"  </Collection>"<<endl;
  pvdfile_new<<"</VTKFile>"<<endl;
  pvdfile_new.close();
  }
}
//-------------------------------------------------------
// Output_pvtu_file
//-------------------------------------------------------
void Level_set_binary_writer::Output_pvti_file_lset(int n_out, p_Levelset_levelinfo lset_level_info, communicator &world)
{
  FILE    *fp;
  char    filename[256];
  sprintf(filename,"%s%d%s","./outdata/mesh_lset_step_",n_out,".pvti");

  int         x1 = 0;
  int         y1 = 0;
  int         z1 = 0;
  int         x2 = lset_level_info->num_pkg.i*ICPX;
  int         y2 = lset_level_info->num_pkg.j*ICPY;
  int         z2 = lset_level_info->num_pkg.k*ICPZ;
  int         np = (x2+1-x1)*(y2+1-y1)*(z2+1-z1);
  int      ncell = (x2-x1)*(y2-y1)*(z2-z1);
  Real   spacing = lset_level_info->dl;
  my_real origin = lset_level_info->box_l;

  ofstream     pvtifile(filename, ios::trunc);

  pvtifile<<"<?xml version=\"1.0\"?>"<<endl;
  pvtifile<<"<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  pvtifile<<" <PImageData WholeExtent=\" "<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<" "<<z1<<" "<<z2<<"\" GhostLevel=\"0\" Origin=\" "<<origin.i<<" "<<origin.j<<" "<<origin.k<<"\" Spacing=\" "<<spacing<<" "<<spacing<<" "<<spacing<<"\">"<<endl;
  pvtifile<<"  <PPointData Vectors=\"vector\">"<<endl;
  if (print_info["phi"      ]) pvtifile<<"    <PDataArray type=\"Float32\" Name=\"phi\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["curv"     ]) pvtifile<<"    <PDataArray type=\"Float32\" Name=\"curv\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["scale"    ]) pvtifile<<"    <PDataArray type=\"Float32\" Name=\"scale\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["norm"     ]) pvtifile<<"    <PDataArray type=\"Float32\" Name=\"norm\" NumberOfComponents=\"3\"/>"<<endl;
  if (print_info["psi"      ]) pvtifile<<"    <PDataArray type=\"Float32\" Name=\"psi\" NumberOfComponents=\"1\"/>"<<endl;
  pvtifile<<"  </PPointData>"<<endl;
  pvtifile<<"  <PCellData Vectors=\"vector\">"<<endl;
  if (print_info["interface"]) pvtifile<<"    <PDataArray type=\"Int32\" Name=\"interface\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["char"     ]) pvtifile<<"    <PDataArray type=\"Int32\" Name=\"char\" NumberOfComponents=\"1\"/>"<<endl;
  #ifdef _MPI_
  if (print_info["color"     ]) pvtifile<<"    <PDataArray type=\"Int32\" Name=\"color\" NumberOfComponents=\"1\"/>"<<endl;
  #endif
  if (print_info["vol"      ]) pvtifile<<"    <PDataArray type=\"Float32\" Name=\"vol\" NumberOfComponents=\"1\"/>"<<endl;
  if (print_info["area"     ]) pvtifile<<"    <PDataArray type=\"Float32\" Name=\"area\" NumberOfComponents=\"1\"/>"<<endl;
  pvtifile<<"  </PCellData>"<<endl;
  for (int i=0; i < world.size(); i++){
    pvtifile<<"    <Piece Extent=\" "<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<" "<<z1<<" "<<z2<<"\" Source=\"./mesh_lset_rank_"<<i<<"_step_"<<n_out<<".vti\"/>"<<endl;
  }
  pvtifile<<" </PImageData>"<<endl;
  pvtifile<<"</VTKFile>"<<endl;

  pvtifile.close();
}
