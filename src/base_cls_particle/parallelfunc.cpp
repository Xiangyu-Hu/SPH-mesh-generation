#include "parallelfunc.h"
#include "glbfunc.h"
//-------------------------------
// SPH:: get minmal timestep
//-------------------------------

class ApplyGetTimestep{
  Particle **my_particle;
public:
  Real min_dt;
  Real max_v;
  void operator()( const blocked_range<int>& r ) {
    Particle **Particle_temp = my_particle;

    for( int num=r.begin(); num!=r.end(); ++num) 
    {
      Real local_t = Particle_temp[num]->timestep;
#ifndef _TRANSPORT_VELOCITY_
      Real local_max_v = get_distance(Particle_temp[num]->v);
#endif
#ifdef _TRANSPORT_VELOCITY_
      Real local_max_v = get_distance(Particle_temp[num]->v_tv);
#endif
      min_dt = AMIN1(min_dt, local_t);
      max_v = AMAX1(max_v, local_max_v);
    }
  }

  ApplyGetTimestep(ApplyGetTimestep &x, split)
  {
    my_particle = x.my_particle;
    min_dt = LONG_MAX;
    max_v = 0.0;
  }

  void join( const ApplyGetTimestep& y)
  {
    min_dt = AMIN1(min_dt, y.min_dt);
    max_v = AMAX1(max_v, y.max_v);
  }

  ApplyGetTimestep(Particle **Particle_temp)
  {
    my_particle = Particle_temp;
    min_dt = LONG_MAX;
    max_v = 0.0;
  }
};

void ParallelGetTimestep(int n, Particle **Particle_temp, Real &local_current_timestep, Real &local_max_velo)
{
  ApplyGetTimestep GetTimestep(Particle_temp);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<int>(0, n), GetTimestep, ap);
  local_current_timestep = GetTimestep.min_dt;
  local_max_velo = GetTimestep.max_v;

}

//------------------------------------------
// SPH:: accumulate the force calculation
//------------------------------------------

class ApplyAccumulation{
  Particle **my_particle;
public:
  my_real sum;
  void operator()( const blocked_range<int>& r ) {
    Particle **Particle_temp = my_particle;

    for( int num=r.begin(); num!=r.end(); ++num) 
    {
      my_real local_force;
#ifdef _ALE_
      Real mass = 1.;
#else
      Real mass = Particle_temp[num]->mass;
#endif
      local_force.i = Particle_temp[num]->F1.i*mass;
      local_force.j = Particle_temp[num]->F1.j*mass;
      local_force.k = Particle_temp[num]->F1.k*mass;
      
      sum.i += local_force.i;
      sum.j += local_force.j;
      sum.k += local_force.k;
    }
  }

  ApplyAccumulation(ApplyAccumulation &x, split)
  {
    my_particle = x.my_particle;
    sum.i = 0.;
    sum.j = 0.;
    sum.k = 0.;
  }

  void join( const ApplyAccumulation& y)
  {
    sum.i += y.sum.i;
    sum.j += y.sum.j;
    sum.k += y.sum.k;
  }

  ApplyAccumulation(Particle **Particle_temp)
  {
    my_particle = Particle_temp;
    sum.i = 0.;
    sum.j = 0.;
    sum.k = 0.;
  }
};

void ParallelAccumulation(int n, Particle **Particle_temp, my_real &local_total_force)
{
  ApplyAccumulation Accumulation(Particle_temp);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<int>(0, n), Accumulation, ap);
  local_total_force.i = Accumulation.sum.i;
  local_total_force.j = Accumulation.sum.j;
  local_total_force.k = Accumulation.sum.k;
}
#ifdef _ALE_
//------------------------------------------
// SPH:: get glbl mass energy
//------------------------------------------

class ApplyParallelmassenergy{
  Particle **my_particle;
public:
  Real sum_mass;
  Real sum_e;
  void operator()( const blocked_range<int>& r ) {
    Particle **Particle_temp = my_particle;

    for( int num=r.begin(); num!=r.end(); ++num)
    {
      sum_mass += Particle_temp[num]->mass;
      sum_e += Particle_temp[num]->U3;
    }
  }

  ApplyParallelmassenergy(ApplyParallelmassenergy &x, split)
  {
    my_particle = x.my_particle;
    sum_mass = 0.;
    sum_e = 0.;
  }

  void join( const ApplyParallelmassenergy& y)
  {
    sum_mass += y.sum_mass;
    sum_e += y.sum_e;
  }

  ApplyParallelmassenergy(Particle **Particle_temp)
  {
    my_particle = Particle_temp;
    sum_mass = 0.;
    sum_e = 0.;
  }
};

void Parallelmassenergy(int n, Particle **Particle_temp,  Real &local_total_mass,  Real &local_total_e)
{
  ApplyParallelmassenergy Parallelmassenergy(Particle_temp);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<int>(0, n), Parallelmassenergy, ap);
  local_total_mass = Parallelmassenergy.sum_mass;
  local_total_e = Parallelmassenergy.sum_e;
}
#endif
//------------------------------------------
// SPH:: accumulate the kinetic energy
//------------------------------------------

class ApplyAccumulationEnergy{
  Particle **my_particle;
public:
  Real total;
  void operator()( const blocked_range<int>& r ) {
    Particle **Particle_temp = my_particle;

    for( int num=r.begin(); num!=r.end(); ++num) 
    {
      my_real local_velocity;
      local_velocity.i = Particle_temp[num]->v.i;
      local_velocity.j = Particle_temp[num]->v.j;
      local_velocity.k = Particle_temp[num]->v.k;
      
      total += 0.5*(local_velocity.i*local_velocity.i + local_velocity.j*local_velocity.j + local_velocity.k*local_velocity.k);
    }
  }

  ApplyAccumulationEnergy(ApplyAccumulationEnergy &x, split)
  {
    my_particle = x.my_particle;
    total = 0.;
  }

  void join( const ApplyAccumulationEnergy& y)
  {
    total += y.total;
  }

  ApplyAccumulationEnergy(Particle **Particle_temp)
  {
    my_particle = Particle_temp;
    total = 0.;
  }
};

void ParallelAccumulationEnergy(int n, Particle **Particle_temp, Real &local_total_energy)
{
  ApplyAccumulationEnergy AccumulationEnergy(Particle_temp);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<int>(0, n), AccumulationEnergy, ap);
  local_total_energy = AccumulationEnergy.total;
}

//------------------------------------------
// SPH:: get mean velocity
//------------------------------------------

class ApplyGetMeanVelocity{
  Particle **my_particle;
public:
  my_real v_avg;
  void operator()( const blocked_range<int>& r ) {
    Particle **Particle_temp = my_particle;

    for( int num=r.begin(); num!=r.end(); ++num) 
    {
      my_real local_velocity;
      local_velocity.i = Particle_temp[num]->v.i;
      local_velocity.j = Particle_temp[num]->v.j;
      local_velocity.k = Particle_temp[num]->v.k;
      
      v_avg.i += local_velocity.i;
      v_avg.j += local_velocity.j;
      v_avg.k += local_velocity.k;
    }
  }

  ApplyGetMeanVelocity(ApplyGetMeanVelocity &x, split)
  {
    my_particle = x.my_particle;
     v_avg.i = 0.;
     v_avg.j = 0.;
     v_avg.k = 0.;
  }

  void join( const ApplyGetMeanVelocity& y)
  {
    v_avg.i += y.v_avg.i;
    v_avg.j += y.v_avg.j;
    v_avg.k += y.v_avg.k;
  }

  ApplyGetMeanVelocity(Particle **Particle_temp)
  {
    my_particle = Particle_temp;
    v_avg.i = 0.;
    v_avg.j = 0.;
    v_avg.k = 0.;
  }
};

void ParallelGetMeanVelocity(int n, Particle **Particle_temp, my_real &v_avg)
{
  ApplyGetMeanVelocity GetMeanVelocity(Particle_temp);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<int>(0, n), GetMeanVelocity, ap);
  v_avg.i = GetMeanVelocity.v_avg.i/n;
  v_avg.j = GetMeanVelocity.v_avg.j/n;
  v_avg.k = GetMeanVelocity.v_avg.k/n;
}
//------------------------------------------
// SPH:: get mean velocity
//------------------------------------------

class ApplyGetMassCenter{
  Particle **my_particle;
public:
  my_real mass_center;
  void operator()( const blocked_range<int>& r ) {
    Particle **Particle_temp = my_particle;

    for( int num=r.begin(); num!=r.end(); ++num) 
    {
      my_real local_coord;
      local_coord.i = Particle_temp[num]->coord.i;
      local_coord.j = Particle_temp[num]->coord.j;
      local_coord.k = Particle_temp[num]->coord.k;
      Real mass = Particle_temp[num]->p_mass;
      
      mass_center.i += mass*local_coord.i;
      mass_center.j += mass*local_coord.j;
      mass_center.k += mass*local_coord.k;
    }
  }

  ApplyGetMassCenter(ApplyGetMassCenter &x, split)
  {
    my_particle = x.my_particle;
     mass_center.i = 0.;
     mass_center.j = 0.;
     mass_center.k = 0.;
  }

  void join( const ApplyGetMassCenter& y)
  {
    mass_center.i += y.mass_center.i;
    mass_center.j += y.mass_center.j;
    mass_center.k += y.mass_center.k;
  }

  ApplyGetMassCenter(Particle **Particle_temp)
  {
    my_particle = Particle_temp;
     mass_center.i = 0.;
     mass_center.j = 0.;
     mass_center.k = 0.;
  }
};

void ParallelGetMassCenter(int n, Real total_mass, Particle **Particle_temp, my_real &mass_center)
{
  ApplyGetMassCenter GetMassCenter(Particle_temp);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<int>(0, n), GetMassCenter, ap);
  mass_center.i = GetMassCenter.mass_center.i/total_mass;
  mass_center.j = GetMassCenter.mass_center.j/total_mass;
  mass_center.k = GetMassCenter.mass_center.k/total_mass;
}

//------------------------------------------
// SPH:: accumulation
//------------------------------------------

class ApplyAccumulationSum{
  my_real *my_force_accumulate;
public:
  my_real sum;
  void operator()( const blocked_range<int>& r ) {
    my_real *force_accumulate = my_force_accumulate;
    for( int num=r.begin(); num!=r.end(); ++num) 
    {
      sum.i += force_accumulate[num].i;
      sum.j += force_accumulate[num].j;
      sum.k += force_accumulate[num].k;
    }
  }

  ApplyAccumulationSum(ApplyAccumulationSum &x, split)
  {
    my_force_accumulate = x.my_force_accumulate;
    sum.i = 0.;
    sum.j = 0.;
    sum.k = 0.;
  }

  void join( const ApplyAccumulationSum& y)
  {
    sum.i += y.sum.i;
    sum.j += y.sum.j;
    sum.k += y.sum.k;
  }

  ApplyAccumulationSum(my_real *force_accumulate)
  {
    my_force_accumulate = force_accumulate;
    sum.i = 0.;
    sum.j = 0.;
    sum.k = 0.;
  }
};

void ParallelAccumulationSum(int n, my_real *force_accumulate, my_real &local_total_force)
{
  ApplyAccumulationSum AccumulationSum(force_accumulate);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<int>(0, n), AccumulationSum, ap);
  local_total_force.i += AccumulationSum.sum.i;
  local_total_force.j += AccumulationSum.sum.j;
  local_total_force.k += AccumulationSum.sum.k;
}
#ifdef _MPI_
//------------------------------------------
// SPH:: get local boundary to build
//  local tree
//------------------------------------------
class ApplyGetLocalBound{
  Particle **my_particle;
public:
  my_real    local_box_l;
  my_real    local_box_r;

  void operator()( const blocked_range<int>& r ) {
    Particle **Particle_temp = my_particle;
    for( int num=r.begin(); num!=r.end(); ++num){
      my_real coord;
      my_set_data (coord, Particle_temp[num]->coord);
      local_box_l.i = DIM_X==1 ? AMIN1(local_box_l.i, coord.i) : 0.;
      local_box_l.j = DIM_Y==1 ? AMIN1(local_box_l.j, coord.j) : 0.;
      local_box_l.k = DIM_Z==1 ? AMIN1(local_box_l.k, coord.k) : 0.;
      local_box_r.i = DIM_X==1 ? AMAX1(local_box_r.i, coord.i) : 0.;
      local_box_r.j = DIM_Y==1 ? AMAX1(local_box_r.j, coord.j) : 0.;
      local_box_r.k = DIM_Z==1 ? AMAX1(local_box_r.k, coord.k) : 0.;
    }
  }
  ApplyGetLocalBound (ApplyGetLocalBound &x, split)
  {
    my_particle = x.my_particle;
    local_box_l.i = LONG_MAX;
    local_box_l.j = LONG_MAX;
    local_box_l.k = LONG_MAX;
    local_box_r.i = LONG_MIN;
    local_box_r.j = LONG_MIN;
    local_box_r.k = LONG_MIN;
  }
  void join( const ApplyGetLocalBound &y)
  {
    local_box_l.i = AMIN1(local_box_l.i, y.local_box_l.i);
    local_box_l.j = AMIN1(local_box_l.j, y.local_box_l.j);
    local_box_l.k = AMIN1(local_box_l.k, y.local_box_l.k);
    local_box_r.i = AMAX1(local_box_r.i, y.local_box_r.i);
    local_box_r.j = AMAX1(local_box_r.j, y.local_box_r.j);
    local_box_r.k = AMAX1(local_box_r.k, y.local_box_r.k);
  }
  ApplyGetLocalBound (Particle **Particle_temp)
  {
    my_particle = Particle_temp;
    local_box_l.i = LONG_MAX;
    local_box_l.j = LONG_MAX;
    local_box_l.k = LONG_MAX;
    local_box_r.i = LONG_MIN;
    local_box_r.j = LONG_MIN;
    local_box_r.k = LONG_MIN;
  }
};

void ParallelGetLocalBound(int n, Particle **Particle_temp, my_real &local_box_l, my_real &local_box_r)
{
  ApplyGetLocalBound GetLocalBound(Particle_temp);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<int>(0, n), GetLocalBound, ap);
  local_box_l.i = GetLocalBound.local_box_l.i;
  local_box_l.j = GetLocalBound.local_box_l.j;
  local_box_l.k = GetLocalBound.local_box_l.k;
  local_box_r.i = GetLocalBound.local_box_r.i;
  local_box_r.j = GetLocalBound.local_box_r.j;
  local_box_r.k = GetLocalBound.local_box_r.k;
}

//------------------------------------------
// SPH:: get maximal degree
//------------------------------------------

class ApplyGetDegree{
   Graph *my_graph;
public:
  int maximal_degree;
    void operator()( const blocked_range<VertexIterator>& r ) {
  Graph *graph = my_graph;
    for (VertexIterator vit = r.begin(); vit != r.end(); ++vit) 
    {
      VertexDescriptor vd = *vit; 
      maximal_degree = AMAX1(maximal_degree,int(degree(vd,*graph)));
    }
  }

  ApplyGetDegree(ApplyGetDegree &x, split)
  {
    my_graph = x.my_graph;
    maximal_degree = 0;
  }

  void join( const ApplyGetDegree& y)
  {
    maximal_degree = AMAX1(maximal_degree, y.maximal_degree);
  }

  ApplyGetDegree(Graph *graph)
  {
    my_graph = graph;
    maximal_degree = 0;
  }
};

void ParallelGetDegree(VertexIterator vit_f, VertexIterator vit_s, Graph *graph, int &maximal_degree)
{
  ApplyGetDegree GetDegree(graph);
  static affinity_partitioner ap;
  parallel_reduce( blocked_range<VertexIterator>(vit_f, vit_s), GetDegree, ap);
  maximal_degree = GetDegree.maximal_degree;

}
#endif
