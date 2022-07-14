#ifndef PARALLELFUNC_H
#define PARALLELFUNC_H

#include "glbcls.h"
#include "float.h"
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

using namespace std;
using namespace tbb;

#define MAXK 10

void ParallelGetTimestep(int n, Particle **Particle_temp, Real &local_current_timestep, Real &local_max_velo);
void ParallelAccumulation(int n, Particle **Particle_temp, my_real &local_total_force);
void Parallelmassenergy(int n, Particle **Particle_temp,  Real &local_total_mass,  Real &local_total_e);
void ParallelAccumulationEnergy(int n, Particle **Particle_temp, Real &local_total_energy);
void ParallelAccumulationSum(int n, my_real *force_accumulate, my_real &local_total_force);
void ParallelGetMassCenter(int n, Real total_mass, Particle **Particle_temp, my_real &mass_center);
void ParallelGetMeanVelocity(int n, Particle **Particle_temp, my_real &v_avg);
#ifdef _MPI_
void ParallelGetLocalBound(int n, Particle **Particle_temp, my_real &local_box_l, my_real &local_box_r);
void ParallelGetDegree(VertexIterator vit_f, VertexIterator vit_s, Graph *graph, int &maximal_degree);
#endif

// radix sort with keys
inline int get_value(int a,int d)
{
  int b=a;
  for (;d>0&&a>0;d--)
    b/=MAXK;
  return b%MAXK;
}
inline int count_sort (concurrent_vector <std::pair <int,p_Particle>> &key, int n,int d)
{
  int k[MAXK] = {0};  
  int * temp,*b;
  p_Particle *temp2;
  int dt =1;
    
   temp = new int[n];
  temp2 = new p_Particle[n];
      b = new int[n];

  if (NULL == temp)
    return 0 ;

  for (int i=0;i<n;i++){
    b[i] = get_value(key[i].first,d);
    k[b[i]]++;
  }

  for (int i=1;i<10;i++)
    k[i]+=k[i-1];

  for (int i=n-1;i>=0;i--){  
    temp[--k[b[i]]]=key[i].first; 
    temp2[k[b[i]]]=key[i].second;
  }

  for (int i = 0; i < n; i++){
    key[i].first = temp[i]; 
    key[i].second = temp2[i];
  }

  delete [] temp;
  delete [] temp2;
  delete [] b;

  return 1 ;  
}

inline void sort_by_key(concurrent_vector <std::pair <int,p_Particle>> &key, int total_num_subcell)
{
  int d = 1;
  int numCell = total_num_subcell;
  while(numCell/10){
    numCell/=10;
    ++d;
  }
  for ( int i=0;i<=d;i++)
    count_sort(key, key.size(), i);
}
#endif

