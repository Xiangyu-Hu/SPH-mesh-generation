#include "geometry_calculator.h"
#include "glbfunc.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <functional>
#include <iostream>

//------------------------------------------------------------
// Get_corner
//------------------------------------------------------------
void GeometryCalculator::Get_corner(Real phii0[1+2*DIM_X][1+2*DIM_Y][1+2*DIM_Z], Real corner[1+2*DIM_X][1+2*DIM_Y][1+2*DIM_Z])
{
  corner[DIM_X][DIM_Y][DIM_Z] = phii0[DIM_X][DIM_Y][DIM_Z];

  int r2 = 1+DIM_X;
  int s2 = 1+DIM_Y;
  int t2 = 1+DIM_Z;

  for(int r=0; r<r2; r++)
    for(int s=0; s<s2; s++)
      for(int t=0; t<t2; t++){
        corner[2*r][2*s][2*t] = 0.125*(  phii0[r][s      ][t      ] + phii0[r+DIM_X][s      ][t      ] 
                                       + phii0[r][s+DIM_Y][t      ] + phii0[r+DIM_X][s+DIM_Y][t      ]
                                       + phii0[r][s      ][t+DIM_Z] + phii0[r+DIM_X][s      ][t+DIM_Z] 
                                       + phii0[r][s+DIM_Y][t+DIM_Z] + phii0[r+DIM_X][s+DIM_Y][t+DIM_Z]);
      }

  #if (DIM_Z==1 && DIM_X+DIM_Y>=1)
  {
    for(int r=0; r<r2; r++)
      for(int s=0; s<s2; s++){
        corner[2*r][2*s][1] = 0.25*(  phii0[r][s      ][1] + phii0[r+DIM_X][s      ][1] 
                                    + phii0[r][s+DIM_Y][1] + phii0[r+DIM_X][s+DIM_Y][1]);
      }    
  }
  #endif

  #if (DIM_X==1 && DIM_Z+DIM_Y>=1)
  {
    for(int s=0; s<s2; s++)
      for(int t=0; t<t2; t++){
        corner[1][2*s][2*t] = 0.25*(  phii0[1][s      ][t] + phii0[1][s      ][t+DIM_Z] 
                                    + phii0[1][s+DIM_Y][t] + phii0[1][s+DIM_Y][t+DIM_Z]);
      }    
  }
  #endif

  #if (DIM_Y==1 && DIM_Z+DIM_X>=1)
  {    
    for(int r=0; r<r2; r++)
      for(int t=0; t<t2; t++){
        corner[2*r][1][2*t] = 0.25*(  phii0[r][1][t      ] + phii0[r+DIM_X][1][t      ] 
                                    + phii0[r][1][t+DIM_Z] + phii0[r+DIM_X][1][t+DIM_Z]);
      }
  }
  #endif

  #if (DIM_Y+DIM_Z+DIM_X==3)
  {    
    for(int r=0; r<2; r++){
      corner[2*r][1  ][1  ] = 0.5*(phii0[r][1][1] + phii0[r+1][  1][  1]);
      corner[1  ][2*r][1  ] = 0.5*(phii0[1][r][1] + phii0[  1][r+1][  1]);
      corner[1  ][1  ][2*r] = 0.5*(phii0[1][1][r] + phii0[  1][  1][r+1]);
    }
  }
  #endif
}

//------------------------------------------------------------
// Get_surface_area using marching cube method
//------------------------------------------------------------
Real GeometryCalculator::Get_surface_area(Real var[1+2*DIM_X][1+2*DIM_Y][1+2*DIM_Z], Real dx, Real dy, Real dz, Real dl)
{
  Real area = 0.;

  #if (DIM_Y+DIM_Z+DIM_X==3)
  {
    int marker;
    Real face[3][3];
    marker =  mark_volume(var);

    if(marker == -1){
      area = get_area_MQ(var, dx, dy, dz);
    }
    else{
      area = 0.0;
    }
  }
  #endif

  #if (DIM_Y+DIM_Z+DIM_X==2)
  {
    int marker;
    Real face[3][3];
    Real var_[3][3][3];
    Real dx_, dy_, dz_;
    dx_ = dx;
    dy_ = dy;
    dz_ = dz;

    #if(DIM_X==0)
    {
      dx_ = 1.;
      for(int r=0; r<3; r++){
        for(int s=0; s<3; s++){
          Real varr = var[0][r][s];
          face[r][s] = varr;
          var_[0][r][s] = varr;
          var_[1][r][s] = varr;
          var_[2][r][s] = varr;
        }
      }
    }
    #endif

    #if(DIM_Y==0)
    {
      dy_ = 1.;
      for(int r=0; r<3; r++){
        for(int s=0; s<3; s++){
          Real varr = var[r][0][s];
          face[r][s] =  varr;
          var_[r][0][s] = varr;
          var_[r][1][s] = varr;
          var_[r][2][s] = varr;
        }
      }
    }
    #endif

    #if(DIM_Z==0)
    {
      dz_ = 1.;
      for(int r=0; r<3; r++){
        for(int s=0; s<3; s++){
          Real varr = var[r][s][0];
          face[r][s] = varr;
          var_[r][s][0] = varr;
          var_[r][s][1] = varr;
          var_[r][s][2] = varr;
        }
      }
    }
    #endif

    //volume fraction
    marker =  mark_face(face);

    if(marker == -1){         
      //record interface node
      area = get_area_MQ(var_, dx_, dy_, dz_);
    }
    else{
      area = 0.0;
    }
  }
  #endif

  return area;
}

//------------------------------------------------------------
// Get_volume_fraction
//------------------------------------------------------------
Real GeometryCalculator::Get_cutcell_volume(Real var[1+2*DIM_X][1+2*DIM_Y][1+2*DIM_Z], Real dx, Real dy, Real dz, Real dl)
{
  Real H = 0.;
  #if (DIM_Y+DIM_Z+DIM_X==3)
  {
    int marker;
    Real face[3][3];
    //volume fraction
    marker =  mark_volume(var);

    if(marker == -1){         
      H = get_volume(var);
    }
    else if(marker == 1){
      H = 1.0;
    }
    else{
      H = 0.0;
    }
  }
  #endif

  #if (DIM_Y+DIM_Z+DIM_X==2)
  {
    int marker;
    Real face[3][3];
    #if(DIM_X==0)
    {
      for(int r=0; r<3; r++){
        for(int s=0; s<3; s++){
          face[r][s] = var[0][r][s];
        }
      }
    }
    #endif

    #if(DIM_Y==0)
    {
      for(int r=0; r<3; r++){
        for(int s=0; s<3; s++){
          face[r][s] =  var[r][0][s];
        }
      }
    }
    #endif

    #if(DIM_Z==0)
    {
      for(int r=0; r<3; r++){
        for(int s=0; s<3; s++){
          face[r][s] = var[r][s][0];
        }
      }
    }
    #endif

    //volume fraction
    marker =  mark_face(face);

    if(marker == -1){         
      //record interface node
      H = get_area(face);
    }
    else if(marker == 1){
      H = 1.0;
    }
    else{
      H = 0.0;
    }
  }
  #endif

  Real volume = (DIM_X==1 ? dx : 1.)*(DIM_Y==1 ? dy : 1.)*(DIM_Z==1 ? dz : 1.);

  return H*volume;
}
//-------------------------------------------------------------------------------------------------
//  Get volume fraction of a 2D cell from level set on cell corners
//-------------------------------------------------------------------------------------------------
Real GeometryCalculator::get_area(Real phii[3][3])
{
  Real piece = 0.25;
  Real area = 0.0;
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
    {
      Real sub_cor[2][2];
        for(int r=0; r<2; r++)
          for(int s=0; s<2; s++)
            sub_cor[r][s] = phii[i+r][j+s];
      area += get_face(0, piece, sub_cor);
    }
  area = AMIN1(area, (Real)1.0);
  return area;
}
//-------------------------------------------------------------------------------------------------
//  Get volume fraction of a 3D cell from level set on cell corners
//-------------------------------------------------------------------------------------------------
Real GeometryCalculator::get_volume(Real phii[3][3][3])
{
  Real vof = 0.0;
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++)
      {
        Real sub_cor[2][2][2];
        for(int r=0; r<2; r++)
          for(int s=0; s<2; s++)
            for(int t=0; t<2; t++)
              sub_cor[r][s][t] = phii[i+r][j+s][k+t];
        vof += get_subvolume(0, 0.125, sub_cor);
      }
  vof = AMIN1(vof, (Real)1.0);
  return vof;
}
//-------------------------------------------------------------------------------------------------
//  Get surface area of a 3D cell from level set on cell corners
//-------------------------------------------------------------------------------------------------
Real GeometryCalculator::get_area_MQ(Real phii[3][3][3], Real dx, Real dy, Real dz)
{
  Real area = 0.0;
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++)
      {
        Polygonise_MQ::GRIDCELL sub_grid;
        std::vector<Polygonise_MQ::TRIANGLE> tries; tries.clear();

        //find all the vertices in current grid
        int vertlist[8][3] = {{0,0,0},{0,1,0},{1,1,0},{1,0,0},
                              {0,0,1},{0,1,1},{1,1,1},{1,0,1}};

        for (int ivert = 0; ivert < 8; ivert++){
          int      vx, vy, vz;
          vx = i + vertlist[ivert][0];
          vy = j + vertlist[ivert][1];
          vz = k + vertlist[ivert][2];
          sub_grid.p[ivert].x = vx*0.5*dx;
          sub_grid.p[ivert].y = vy*0.5*dy;
          sub_grid.p[ivert].z = vz*0.5*dz;

          sub_grid.val[ivert] = phii[vx][vy][vz];
        }

        int ntris = Polygonise_MQ::Polygonise(sub_grid, double(0.0), tries);

        if (ntris > 0){
          std::vector<Polygonise_MQ::TRIANGLE>::iterator ittri;
          for (ittri=tries.begin(); ittri!=tries.end(); ++ittri){
            Real area_ = get_area_tri ((*ittri).p[0], (*ittri).p[1], (*ittri).p[2]);
            area += area_;
          }
        }
      }
  return area;
}
//-------------------------------------------------------------------------------------------------
// mark if this is a full face with a 9-point stencil
//-------------------------------------------------------------------------------------------------
Real GeometryCalculator::get_area_tri(Polygonise_MQ::XYZ p1, Polygonise_MQ::XYZ p2, Polygonise_MQ::XYZ p3)
{
  my_real vector_1;
  vector_1.i = p2.x - p1.x;
  vector_1.j = p2.y - p1.y;
  vector_1.k = p2.z - p1.z;
  my_real vector_2;
  vector_2.i = p3.x - p2.x;
  vector_2.j = p3.y - p2.y;
  vector_2.k = p3.z - p2.z;
  my_real vector_3;
  vector_3.i = p3.x - p1.x;
  vector_3.j = p3.y - p1.y;
  vector_3.k = p3.z - p1.z;

  Real dist_1 = sqrt(vector_1.i*vector_1.i + vector_1.j*vector_1.j + vector_1.k*vector_1.k);
  Real dist_2 = sqrt(vector_2.i*vector_2.i + vector_2.j*vector_2.j + vector_2.k*vector_2.k);
  Real dist_3 = sqrt(vector_3.i*vector_3.i + vector_3.j*vector_3.j + vector_3.k*vector_3.k);

  Real area = 0.5*(dist_3*dist_1)*sin(acos((vector_1.i * vector_3.i + vector_1.j * vector_3.j + vector_1.k * vector_3.k)/(dist_3*dist_1  + 1.e-20)));

  if (std::isnan(area)) area = 0.;
  
  return area;
}
//-------------------------------------------------------------------------------------------------
// mark if this is a full face with a 9-point stencil
//-------------------------------------------------------------------------------------------------
int GeometryCalculator::mark_face(Real cor[3][3])
{
  int flag = 0;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      if(cor[i][j] > 0.0)
        flag += 1;
      else if(cor[i][j] < 0.0)
        flag -= 1;
      else;
    }

  if(flag == 9)
    return 1; //full face
  else if( flag == -9)
    return 0; //empty face
  else
    return -1;  //cut face
}
//-------------------------------------------------------------------------------------------------
//  mark if this is a full face with a 4-point stencil
//-------------------------------------------------------------------------------------------------
int GeometryCalculator::check_face(Real cor[2][2])
{
  int flag = 0;
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
    {
      if(cor[i][j] > 0.0)
        flag += 1;
      else if(cor[i][j] < 0.0)
        flag -= 1;
      else;
    }

  if(flag == 4)
    return 1;
  else if( flag == -4)
    return 0;
  else
    return -1;
}
//-------------------------------------------------------------------------------------------------
//  mark if this is a full cell with a 27-point stencil
//-------------------------------------------------------------------------------------------------
int GeometryCalculator::mark_volume(Real cor[3][3][3])
{
  int flag = 0;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
      {
        if(cor[i][j][k] > 0.0)
          flag += 1;
        else if(cor[i][j][k] < 0.0)
          flag -= 1;
        else;
      }

  if(flag == 27)
    return 1;
  else if( flag == -27)
    return 0;
  else
    return -1;
}
//-------------------------------------------------------------------------------------------------
//  mark if this is a full cell with a 8-point stencil
//-------------------------------------------------------------------------------------------------
int GeometryCalculator::check_volume(Real cor[2][2][2])
{
  int flag = 0;
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      for(int k=0; k<2; k++)
      {
        if(cor[i][j][k] > 0.0)
          flag += 1;
        else if(cor[i][j][k] < 0.0)
          flag -= 1;
        else;
      }

  if(flag == 8)
    return 1;
  else if( flag == -8)
    return 0;
  else
    return -1;
}
//-------------------------------------------------------------------------------------------------
//  Get area of a sub-face with adaptive refinement
//-------------------------------------------------------------------------------------------------
Real GeometryCalculator::get_face(int level_phi, Real multi, Real cor[2][2])
{
  int flag = check_face(cor);
  if (flag == 1)
    return multi;
  else if(flag == 0)
    return 0.0;
  else
  {
    Real vof = 0.0;
    if(level_phi<=5)
    {
      multi *=0.25;
      level_phi += 1;
      Real rf_cor[3][3];
      //step1
      for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
          rf_cor[2*i][2*j] = cor[i][j];
      //step2
      for(int i=0; i<3; i=i+2)
        rf_cor[i][1] = 0.5*(rf_cor[i][0] + rf_cor[i][2]);
      //step3
      for(int j=0; j<3; j++)
        rf_cor[1][j] = 0.5*(rf_cor[0][j] + rf_cor[2][j]);

      for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
        {
          Real sub_cor[2][2];
          for(int r=0; r<2; r++)
            for(int s=0; s<2; s++)
              sub_cor[r][s] = rf_cor[i+r][j+s];
          vof += get_face(level_phi, multi, sub_cor);
        }
    }
    else
    {
      Real center_value = 0.0;
      for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
          center_value += cor[i][j];

      if(center_value>0.0)
        vof = 0.75*multi;
      else if(center_value==0.0)
        vof = 0.5*multi;
      else
        vof = 0.25*multi;
    }
    return vof;
  }
}
//-------------------------------------------------------------------------------------------------
//  Get volume of a sub-cell with adaptive refinement
//-------------------------------------------------------------------------------------------------
Real GeometryCalculator::get_subvolume(int level_phi, Real multi, Real cor[2][2][2])
{
  int flag = check_volume(cor);
  if (flag == 1)
    return multi;
  else if(flag == 0)
    return 0.0;
  else
  {
    Real vof = 0.0;
    if(level_phi<=5)
    {
      multi *=0.125;
      level_phi += 1;
      Real rf_cor[3][3][3];
      //step1
      for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
          for(int k=0; k<2; k++)
            rf_cor[2*i][2*j][2*k] = cor[i][j][k];
      //step2
      for(int i=0; i<3; i=i+2)
        for(int j=0; j<3; j=j+2)
          rf_cor[i][j][1] = 0.5*(rf_cor[i][j][0] + rf_cor[i][j][2]);
      //step3
      for(int k=0; k<3; k++)
        for(int i=0; i<3; i=i+2)
          rf_cor[i][1][k] = 0.5*(rf_cor[i][0][k] + rf_cor[i][2][k]);

      //step4
      for(int j=0; j<3; j++)
        for(int k=0; k<3; k++)
          rf_cor[1][j][k] = 0.5*(rf_cor[0][j][k] + rf_cor[2][j][k]);

      for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
          for(int k=0; k<2; k++)
          {
            Real sub_cor[2][2][2];
            for(int r=0; r<2; r++)
              for(int s=0; s<2; s++)
                for(int t=0; t<2; t++)
                  sub_cor[r][s][t] = rf_cor[i+r][j+s][k+t];
            vof += get_subvolume(level_phi, multi, sub_cor);
          }
    }
    else
    {
      Real center_value = 0.0;
      for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
          for(int k=0; k<2; k++)
            center_value += cor[i][j][k];

      if(center_value>0.0)
        vof = 0.75*multi;
      else if(center_value==0.0)
        vof = 0.5*multi;
      else
        vof = 0.25*multi;
    }
    return vof;
  }
}

//-------------------------------------------------------------------------------------------------
//  Get the length of a segment intersection with an AABB
//-------------------------------------------------------------------------------------------------
Real GeometryCalculator::Get_segment_length(my_real min, my_real max, my_real v0, my_real v1)
{
  Real dist = 0.;

  double min_[3] = {min.i,min.j,min.k};
  double max_[3] = {max.i,max.j,max.k};
  double v0_ [3] = { v0.i, v0.j, v0.k};
  double v1_ [3] = { v1.i, v1.j, v1.k};
  double p0_ [3] = {   0.,   0.,   0.};
  double p1_ [3] = {   0.,   0.,   0.};

  if (AABB_calculator::Segment_AABB_Intersection(min_, max_, v0_, v1_, p0_, p1_)){
    double pp_ [3] = {p1_[0]-p0_[0], p1_[1]-p0_[1], p1_[2]-p0_[2]};
    dist = sqrt (pp_[0]*pp_[0] + pp_[1]*pp_[1] + pp_[2]*pp_[2]);
  }

  return dist;
}