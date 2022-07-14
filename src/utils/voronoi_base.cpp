#include "voronoi_base.h"
#include "glbfunc.h"

//-------------------------------------------------------
// Initialization for Vo_particle
//-------------------------------------------------------
void Vo_particle::Initialize(int index)
{
  id     = index;
  color  = index;
  mass   = 0.;
  P      = 0.;
  rho    = 0.;
  vol    = 0.;
  h      = 0.;
  h_min  = 0.;

  my_set_const (coord, 0.);
  my_set_const (a, 0.);
  my_set_const (f, 0.);
  my_set_const (v, 0.);
  my_set_const (CVT_shift, 0.);
  my_set_const (mass_center, 0.);

}
//-------------------------------------------------------
// Reset VP information
//-------------------------------------------------------
void Vo_particle::Reset()
{
  mass    = 0.;
  energy  = 0.;
  my_set_const (mass_center, 0.);
}