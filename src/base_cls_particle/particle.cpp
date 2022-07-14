#include "glbfunc.h"
#include "sph.h"
#include "level_infor.h"
#include "particle.h"

/***************************************************/
/*                                                 */
/*      Functions defined in class "Particle"      */
/*                                                 */
/***************************************************/
//-----------------------------------------------------------------
// Get the velocity of particle for Taylor Green Vortex simulation
// t = 0s
//-----------------------------------------------------------------
void Particle_base::Get_velocity_for_TGV()
{
  my_set_const (v, 0.0);
#if DIM ==  2
  v.i = -1*cos(2*PI*coord.i)*sin(2*PI*coord.j);
  v.j = sin(2*PI*coord.i)*cos(2*PI*coord.j);
  v.k = 0.0;
#elif DIM == 3
  v.i = cos(2*PI*coord.i)*sin(2*PI*coord.j)*sin(2*PI*coord.k);
  v.j = -0.5*sin(2*PI*coord.i)*cos(2*PI*coord.j)*sin(2*PI*coord.k);
  v.k = -0.5*sin(2*PI*coord.i)*sin(2*PI*coord.j)*cos(2*PI*coord.k);
#endif
}
//-----------------------------------------------------
// get the value of kernel function
//-----------------------------------------------------
Real Particle_base::Kernel_function_1D(Real dist, Real h)
{
#ifdef _HYPERBOLIC_KERNEL_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 1.0/7.0/powern(h, 1);

  if (s < 1)
    val = C * ( 6.0 - 6.0*s + s*s*s );
  else if (s < 2)
    val = C * (2-s) * (2-s) * (2-s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _CUBIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 0.666666667 / powern(h, 1);

  if (s < 1)
    val = C * ( 1 - 1.5*s*s + 0.75*s*s*s );
  else if (s < 2)
    val = C/4.0 * (2-s) * (2-s) * (2-s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _WINDLAND_C4_
    cout<<"WendlandQuintic: Dim 1 not supported\n";

  return 0.;
#endif
#ifdef _GAUSSIAN_
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);

  Real   C = 0.;

  C = fac/h;

  if (s < 3)
    val = exp(-1.*s*s)*C;
  else
    val = 0.0;

  return val;
#endif
#ifdef _QUINTIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  C = 1.0 / 120. / powern(h, 1);

  if (s < 1)
    val = C * (powern ((3. - s), 5) - 6.*powern ((2. - s), 5) + 15.*powern ((1. - s), 5));
  else if (s < 2)
    val = C * (powern ((3. - s), 5) - 6.*powern ((2. - s), 5));
  else if (s < 3)
    val = C * (powern ((3. - s), 5));
  else
    val = 0.0;

  return val;
#endif
}
//-----------------------------------------------------
// get the value of kernel function
//-----------------------------------------------------
Real Particle_base::Kernel_function_2D(Real dist, Real h)
{
#ifdef _HYPERBOLIC_KERNEL_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 1.0/(3.0*PI) / powern(h, 2);

  if (s < 1)
    val = C * ( 6.0 - 6.0*s + s*s*s );
  else if (s < 2)
    val = C * (2-s) * (2-s) * (2-s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _CUBIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 0.454728408 / powern(h, 2);

  if (s < 1)
    val = C * ( 1 - 1.5*s*s + 0.75*s*s*s );
  else if (s < 2)
    val = C/4.0 * (2-s) * (2-s) * (2-s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _WINDLAND_C4_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  C = 1.0 / PI / powern(h, 2) * 7. / 4.;

  if (s < 2)
    val = C * powern((1. - 0.5*s), 4)*(1. + 2.*s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _GAUSSIAN_
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);
  fac *= 0.5*2/sqrt(PI);

  Real   C = 0.;

  C = fac / powern(h, 2);

  if (s < 3)
    val = exp(-1.*s*s)*C;
  else
    val = 0.0;

  return val;
#endif
#ifdef _QUINTIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  C = 1.0 / PI * 7./478. / powern(h, 2);

  if (s < 1)
    val = C * (powern ((3. - s), 5) - 6.*powern ((2. - s), 5) + 15.*powern ((1. - s), 5));
  else if (s < 2)
    val = C * (powern ((3. - s), 5) - 6.*powern ((2. - s), 5));
  else if (s < 3)
    val = C * (powern ((3. - s), 5));
  else
    val = 0.0;

  return val;
#endif
}
//-----------------------------------------------------
// get the value of kernel function
//-----------------------------------------------------
Real Particle_base::Kernel_function_3D(Real dist, Real h)
{
#ifdef _HYPERBOLIC_KERNEL_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  C = 15.0/(62.0*PI) / powern(h, 3);

  if (s < 1)
    val = C * ( 6.0 - 6.0*s + s*s*s );
  else if (s < 2)
    val = C * (2-s) * (2-s) * (2-s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _CUBIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  C = 0.318309886 / powern(h, 3);

  if (s < 1)
    val = C * ( 1 - 1.5*s*s + 0.75*s*s*s );
  else if (s < 2)
    val = C/4.0 * (2-s) * (2-s) * (2-s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _WINDLAND_C4_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  C = 1.0 / PI / powern(h, 3) * 21. / 16.;

  if (s < 2)
    val = C * powern((1. - 0.5*s), 4)*(1. + 2.*s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _GAUSSIAN_
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);
  fac *= 0.5*2/sqrt(PI);
  fac *= 0.5*2/sqrt(PI);

  Real   C = 0.;
  C = fac / powern(h, 3);

  if (s < 3)
    val = exp(-1.*s*s)*C;
  else
    val = 0.0;

  return val;
#endif
#ifdef _QUINTIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  C = 1.0 / PI * 3./359. / powern(h, 3);

  if (s < 1)
    val = C * (powern ((3. - s), 5) - 6.*powern ((2. - s), 5) + 15.*powern ((1. - s), 5));
  else if (s < 2)
    val = C * (powern ((3. - s), 5) - 6.*powern ((2. - s), 5));
  else if (s < 3)
    val = C * (powern ((3. - s), 5));
  else
    val = 0.0;

  return val;
#endif
}
//-----------------------------------------------------
// get the value of kernel function
//-----------------------------------------------------
Real Particle_base::Kernel_function(Real dist, Real h)
{
#ifdef _HYPERBOLIC_KERNEL_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 1.0/7.0/powern(h, 1);
  else if (DIM == 2)
    C = 1.0/(3.0*PI) / powern(h, 2);
  else if (DIM == 3)
    C = 15.0/(62.0*PI) / powern(h, 3);

  if (s < 1)
    val = C * ( 6.0 - 6.0*s + s*s*s );
  else if (s < 2)
    val = C * (2-s) * (2-s) * (2-s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _CUBIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 0.666666667 / powern(h, 1);
  else if (DIM == 2)
    C = 0.454728408 / powern(h, 2);
  else if (DIM == 3)
    C = 0.318309886 / powern(h, 3);

  if (s < 1)
    val = C * ( 1 - 1.5*s*s + 0.75*s*s*s );
  else if (s < 2)
    val = C/4.0 * (2-s) * (2-s) * (2-s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _WINDLAND_C4_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    cout<<"WendlandQuintic: Dim 1 not supported\n";
  else if (DIM == 2)
    C = 1.0 / PI / powern(h, 2) * 7. / 4.;
  else if (DIM == 3)
    C = 1.0 / PI / powern(h, 3) * 21. / 16.;

  if (s < 2)
    val = C * powern((1. - 0.5*s), 4)*(1. + 2.*s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _GAUSSIAN_
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);
  if (DIM > 1)
    fac *= 0.5*2/sqrt(PI);
  if (DIM > 2)
    fac *= 0.5*2/sqrt(PI);

  Real   C = 0.;

  if (DIM == 1)
    C = fac/h;
  else if (DIM == 2)
    C = fac / powern(h, 2);
  else if (DIM == 3)
    C = fac / powern(h, 3);

  if (s < 3)
    val = exp(-1.*s*s)*C;
  else
    val = 0.0;

  return val;
#endif
#ifdef _QUINTIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 1.0 / 120. / powern(h, 1);
  else if (DIM == 2)
    C = 1.0 / PI * 7./478. / powern(h, 2);
  else if (DIM == 3)
    C = 1.0 / PI * 3./359. / powern(h, 3);

  if (s < 1)
    val = C * (powern ((3. - s), 5) - 6.*powern ((2. - s), 5) + 15.*powern ((1. - s), 5));
  else if (s < 2)
    val = C * (powern ((3. - s), 5) - 6.*powern ((2. - s), 5));
  else if (s < 3)
    val = C * (powern ((3. - s), 5));
  else
    val = 0.0;

  return val;
#endif
}
//-----------------------------------------------------
// get the derivative value of kernel function
//-----------------------------------------------------
Real Particle_base::Derivative_kernel_function_1D(Real dist, my_real dr, Real h)
{
#ifdef _HYPERBOLIC_KERNEL_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 1.0/7.0/powern(h, 2);

  if (s < 1)
    val = C * (-6.0 + 3.*s*s );
  else if (s < 2)
    val = C * (-3.0) * (2-s) * (2-s);
  else
    val = 0.0;

  val /= (dist + 1.e-20);
  return val;
#endif
#ifdef _CUBIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 0.666666667 / powern(h, 2);

  if (s < 1)
    val = 3*C*(-s + 3*(s*s)/4 );
  else if (s < 2)
    val = -3/4.0*C*(2-s)*(2-s);
  else
    val = 0.0;

  val /= (dist + 1.e-20);
  return val;
#endif
#ifdef _WINDLAND_C4_
    cout<<"WendlandQuintic: Dim 1 not supported\n";
    return 0.;
#endif
#ifdef _GAUSSIAN_
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);

  Real   C = 0.;

    C = fac/h;

  if (s < 3)
    val = -2.*s*exp(-1.*s*s)*C/h/(dist+1.e-15);
  else
    val = 0.0;

  return val;
#endif
#ifdef _QUINTIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 1.0 / 120. / powern(h, 1);

  if (s < 1){
    val = C * (-5. * powern ((3. - s), 4) + 30.*powern ((2. - s), 4) - 75.*powern ((1. - s), 4));
    val = val/h/(dist+1.e-15);
  }
  else if (s < 2){
    val = C * (-5. * powern ((3. - s), 4) + 30.*powern ((2. - s), 4));
    val = val/h/(dist+1.e-15);
  }
  else if (s < 3){
    val = C * (-5. * powern ((3. - s), 4));
    val = val/h/(dist+1.e-15);
  }
  else
    val = 0.0;

  return val;
#endif
}
//-----------------------------------------------------
// get the derivative value of kernel function
//-----------------------------------------------------
Real Particle_base::Derivative_kernel_function_2D(Real dist, my_real dr, Real h)
{
#ifdef _HYPERBOLIC_KERNEL_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 1.0/(3.0*PI) / powern(h, 3);

  if (s < 1)
    val = C * (-6.0 + 3.*s*s );
  else if (s < 2)
    val = C * (-3.0) * (2-s) * (2-s);
  else
    val = 0.0;

  val /= (dist + 1.e-20);
  return val;
#endif
#ifdef _CUBIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 0.454728408 / powern(h, 3);

  if (s < 1)
    val = 3*C*(-s + 3*(s*s)/4 );
  else if (s < 2)
    val = -3/4.0*C*(2-s)*(2-s);
  else
    val = 0.0;

  val /= (dist + 1.e-20);
  return val;
#endif
#ifdef _WINDLAND_C4_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 1.0 / PI / powern(h, 2) * 7. / 4.;

  if (s < 2)
    val = -5.*C*powern((1. - 0.5*s), 3)*s/h/(dist+1.e-15);
  else
    val = 0.0;

  return val;
#endif
#ifdef _GAUSSIAN_
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);
  fac *= 0.5*2/sqrt(PI);

  Real   C = 0.;

    C = fac / powern(h, 2);

  if (s < 3)
    val = -2.*s*exp(-1.*s*s)*C/h/(dist+1.e-15);
  else
    val = 0.0;

  return val;
#endif
#ifdef _QUINTIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 1.0 / PI * 7./478. / powern(h, 2);

  if (s < 1){
    val = C * (-5. * powern ((3. - s), 4) + 30.*powern ((2. - s), 4) - 75.*powern ((1. - s), 4));
    val = val/h/(dist+1.e-15);
  }
  else if (s < 2){
    val = C * (-5. * powern ((3. - s), 4) + 30.*powern ((2. - s), 4));
    val = val/h/(dist+1.e-15);
  }
  else if (s < 3){
    val = C * (-5. * powern ((3. - s), 4));
    val = val/h/(dist+1.e-15);
  }
  else
    val = 0.0;

  return val;
#endif
}
//-----------------------------------------------------
// get the derivative value of kernel function
//-----------------------------------------------------
Real Particle_base::Derivative_kernel_function_3D(Real dist, my_real dr, Real h)
{
#ifdef _HYPERBOLIC_KERNEL_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 15.0/(62.0*PI) / powern(h, 4);

  if (s < 1)
    val = C * (-6.0 + 3.*s*s );
  else if (s < 2)
    val = C * (-3.0) * (2-s) * (2-s);
  else
    val = 0.0;

  val /= (dist + 1.e-20);
  return val;
#endif
#ifdef _CUBIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 0.318309886 / powern(h, 4);

  if (s < 1)
    val = 3*C*(-s + 3*(s*s)/4 );
  else if (s < 2)
    val = -3/4.0*C*(2-s)*(2-s);
  else
    val = 0.0;

  val /= (dist + 1.e-20);
  return val;
#endif
#ifdef _WINDLAND_C4_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 1.0 / PI / powern(h, 3) * 21. / 16.;

  if (s < 2)
    val = -5.*C*powern((1. - 0.5*s), 3)*s/h/(dist+1.e-15);
  else
    val = 0.0;

  return val;
#endif
#ifdef _GAUSSIAN_
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);
  fac *= 0.5*2/sqrt(PI);
  fac *= 0.5*2/sqrt(PI);

  Real   C = 0.;

    C = fac / powern(h, 3);

  if (s < 3)
    val = -2.*s*exp(-1.*s*s)*C/h/(dist+1.e-15);
  else
    val = 0.0;

  return val;
#endif
#ifdef _QUINTIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

    C = 1.0 / PI * 3./359. / powern(h, 3);

  if (s < 1){
    val = C * (-5. * powern ((3. - s), 4) + 30.*powern ((2. - s), 4) - 75.*powern ((1. - s), 4));
    val = val/h/(dist+1.e-15);
  }
  else if (s < 2){
    val = C * (-5. * powern ((3. - s), 4) + 30.*powern ((2. - s), 4));
    val = val/h/(dist+1.e-15);
  }
  else if (s < 3){
    val = C * (-5. * powern ((3. - s), 4));
    val = val/h/(dist+1.e-15);
  }
  else
    val = 0.0;

  return val;
#endif
}
//-----------------------------------------------------
// get the derivative value of kernel function
//-----------------------------------------------------
Real Particle_base::Derivative_kernel_function(Real dist, my_real dr, Real h)
{
#ifdef _HYPERBOLIC_KERNEL_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 1.0/7.0/powern(h, 2);
  else if (DIM == 2)
    C = 1.0/(3.0*PI) / powern(h, 3);
  else if (DIM == 3)
    C = 15.0/(62.0*PI) / powern(h, 4);

  if (s < 1)
    val = C * (-6.0 + 3.*s*s );
  else if (s < 2)
    val = C * (-3.0) * (2-s) * (2-s);
  else
    val = 0.0;

  val /= (dist + 1.e-20);
  return val;
#endif
#ifdef _CUBIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 0.666666667 / powern(h, 2);
  else if (DIM == 2)
    C = 0.454728408 / powern(h, 3);
  else if (DIM == 3)
    C = 0.318309886 / powern(h, 4);

  if (s < 1)
    val = 3*C*(-s + 3*(s*s)/4 );
  else if (s < 2)
    val = -3/4.0*C*(2-s)*(2-s);
  else
    val = 0.0;

  val /= (dist + 1.e-20);
  return val;
#endif
#ifdef _WINDLAND_C4_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    cout<<"WendlandQuintic: Dim 1 not supported\n";
  else if (DIM == 2)
    C = 1.0 / PI / powern(h, 2) * 7. / 4.;
  else if (DIM == 3)
    C = 1.0 / PI / powern(h, 3) * 21. / 16.;

  if (s < 2)
    val = -5.*C*powern((1. - 0.5*s), 3)*s/h/(dist+1.e-15);
  else
    val = 0.0;

  return val;
#endif
#ifdef _GAUSSIAN_
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);
  if (DIM > 1)
    fac *= 0.5*2/sqrt(PI);
  if (DIM > 2)
    fac *= 0.5*2/sqrt(PI);

  Real   C = 0.;

  if (DIM == 1)
    C = fac/h;
  else if (DIM == 2)
    C = fac / powern(h, 2);
  else if (DIM == 3)
    C = fac / powern(h, 3);

  if (s < 3)
    val = -2.*s*exp(-1.*s*s)*C/h/(dist+1.e-15);
  else
    val = 0.0;

  return val;
#endif
#ifdef _QUINTIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 1.0 / 120. / powern(h, 1);
  else if (DIM == 2)
    C = 1.0 / PI * 7./478. / powern(h, 2);
  else if (DIM == 3)
    C = 1.0 / PI * 3./359. / powern(h, 3);

  if (s < 1){
    val = C * (-5. * powern ((3. - s), 4) + 30.*powern ((2. - s), 4) - 75.*powern ((1. - s), 4));
    val = val/h/(dist+1.e-15);
  }
  else if (s < 2){
    val = C * (-5. * powern ((3. - s), 4) + 30.*powern ((2. - s), 4));
    val = val/h/(dist+1.e-15);
  }
  else if (s < 3){
    val = C * (-5. * powern ((3. - s), 4));
    val = val/h/(dist+1.e-15);
  }
  else
    val = 0.0;

  return val;
#endif
}
//-----------------------------------------------------
// get the derivative value of kernel function regarding
// to smooth length using a quintic spline function
//-----------------------------------------------------
Real Particle_base::Derivative_h_kernel_function
(Real dist, Real h)
{
#ifdef _HYPERBOLIC_KERNEL_
  cout<<"kernel type not supported\n";
  return 0.0;
#endif
#ifdef _CUBIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 0.666666667 / powern(h, 2);
  else if(DIM == 2)
    C = 0.454728408 / powern(h, 3);
  else if (DIM == 3)
    C = 0.318309886 / powern(h, 4);

  if (s < 1)
    val = C * ( -1.*DIM + 1.5*(DIM+2.)*s*s - 0.75*(DIM+3.)*s*s*s );
  else if (s < 2)
    val = C/4.0 * ( -1.*DIM*powern((2-s),3) +3.*s*powern((2.-s),2));
  else
    val = 0.0;

  return val;
#endif
#ifdef _WINDLAND_C4_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    cout<<"WendlandQuintic: Dim 1 not supported\n";
  else if (DIM == 2)
    C = 1.0 / PI / powern(h, 2) * 7. / 4.;
  else if (DIM == 3)
    C = 1.0 / PI / powern(h, 3) * 21. / 16.;

  Real w = 0.;
  Real dw = 0.;
  Real tmp = 1. - 0.5*s;
  if (s < 2){
    w   = tmp * tmp * tmp * tmp * (2.0*s + 1.0);
    dw  = -5.0 * s * tmp * tmp * tmp;
    val = -1. * C/h * ( dw*s + w*DIM );
  }
  else
    val = 0.0;

  return val;

#endif
#ifdef _GAUSSIAN_
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);
  if (DIM > 1)
    fac *= 0.5*2/sqrt(PI);
  if (DIM > 2)
    fac *= 0.5*2/sqrt(PI);

  Real   C = 0.;

  if (DIM == 1)
    C = fac/h;
  else if (DIM == 2)
    C = fac / powern(h, 2);
  else if (DIM == 3)
    C = fac / powern(h, 3);

  Real  w = 0.;
  Real dw = 0.;
  if (s < 3){
    w   = exp(-1.*s*s);
    dw  = -2.*s*w;
    val = -1.*C/h*(dw*s + w*DIM);
  }
  else
    val = 0.0;

  return val;
#endif
#ifdef _QUINTIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 1.0 / 120. / powern(h, 1);
  else if (DIM == 2)
    C = 1.0 / PI * 7./478. / powern(h, 2);
  else if (DIM == 3)
    C = 1.0 / PI * 3./359. / powern(h, 3);

  Real  w = 0.;
  Real dw = 0.;
  if (s < 1){
     w  = powern ((3. - s), 5);
     w -= 6. * powern ((2. - s), 5);
     w += 15. * powern ((1. - s), 5);

    dw  = -5. * powern ((3. - s), 4);
    dw += 30. * powern ((2. - s), 4);
    dw -= 75. * powern ((1. - s), 4);
  }
  else if (s < 2){
     w  = powern ((3. - s), 5);
     w -= 6. * powern ((2. - s), 5);

    dw  = -5. * powern ((3. - s), 4);
    dw += 30. * powern ((2. - s), 4);
  }
  else if (s < 3){
     w  = powern ((3. - s), 5);
    dw  = -5. * powern ((3. - s), 4);
  }
  else{
     w  = 0.;
    dw  = 0.;
  }
  val = -1. * C / h * (dw*s + w*DIM);

  return val;
#endif
}
//-----------------------------------------------------
// get particel energy
//-----------------------------------------------------
void Particle_base::Get_energy(Real &iEk, Real &iEp, Real &iEf)
{
  Real v_2 = get_distance (v);
       v_2 = powern(v_2, 2);
  Real x_2 = get_distance (coord);
       x_2 = powern(x_2, 2);
  iEk = 0.5*mass*v_2;
  iEp = vol*P;
  iEf = 0.5*mass*x_2;
}
//-----------------------------------------------------
// Set particle scale
//-----------------------------------------------------
void Particle_base::Set_particle_scale(SPH *sph, int &flag)
{
  Real max_scale = sph->max_scale;
  Real min_scale = sph->min_scale;
  int  stop = 0;
  while (stop == 0){
    if (h*CUT_OFF > scale + 1.e-15){
      scale *= SCALE_RATIO;
      if ( scale > max_scale + 1.e-15)
        flag = -1;
    }else if (h*CUT_OFF < scale/SCALE_RATIO - 1.e-15){
      scale /= SCALE_RATIO;
      if ( scale < min_scale - 1.e-15)
        flag = 1;
    }else{
      stop = 1;
    }
  }
}
#ifdef _MPI_
#ifdef _WEIGHTED_PARTITION_
//-----------------------------------------------------
// Get weighted mass for CVP patitioning
//-----------------------------------------------------
void Particle_base::Get_weighted_patitioning_mass(SPH *sph, int flag, int n_neighbor)
{
  Real weight = sph->weight;
  Real glbl_mass = 1./sph->glbl_total_neighbor_size;
  Real glbl_dist = 1./sph->glbl_total_dist_number;
  if (flag == 0)
    p_mass = weight*mass*glbl_mass + (1.-weight)*p_mass*glbl_dist;
  else if (flag == 1)
    p_mass = weight*n_neighbor*glbl_mass + (1.-weight)*p_mass*glbl_dist;
}
#endif
//-----------------------------------------------------
// Get mass for CVP patitioning
//-----------------------------------------------------
void Particle_base::Normalize_particle_mass(SPH *sph)
{
  p_mass /= sph->glbl_total_mass;
}
//-----------------------------------------------------
// Get mass for CVP patitioning
//-----------------------------------------------------
void Particle_base::Get_patitioning_mass(SPH *sph, int flag, int n_neighbor)
{
  p_mass = 0.;
#ifdef _MASS_
  if (flag == 0)
    p_mass = mass;
  else if (flag == 1)
    p_mass = n_neighbor;
#endif
#if defined(_DIST_NUMBER_) || defined(_WEIGHTED_PARTITION_)

  p_Level_info    current_level = sph->level_info[level-sph->Lmin];
  my_real coord_shift = my_minus_data (coord, sph->box_l);

#ifndef _SCLL_
  my_int  pos = get_cell_id (coord_shift, current_level->dcell, current_level->cell_start, current_level->cell_end);
#else
  my_int  pos        = get_cell_id (coord_shift, current_level->dsubcell, current_level->subcell_start, current_level->subcell_end);
#endif
  current_level->Get_patitioning_mass(sph, this, pos, flag);
#endif
}
//-----------------------------------------------------
// Set particle scale
//-----------------------------------------------------
void Particle_base::Reset_color(int current_color)
{
  color = current_color;
}
//-----------------------------------------------------
// shift particle position
//-----------------------------------------------------
void Particle_base::Shift_coordinate(my_real domain, my_real box_l, my_real box_r)
{
#if P_DIM_X == 1 && PERI_DIM_X == 1
      if (coord.i > box_r.i+1.e-20) coord.i -= domain.i;
      else if (coord.i < box_l.i-1.e-20) coord.i += domain.i;
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 1
      if (coord.j > box_r.j+1.e-20) coord.j -= domain.j;
      else if (coord.j < box_l.j-1.e-20) coord.j += domain.j;
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 1
      if (coord.k > box_r.k+1.e-20) coord.k -= domain.k;
      else if (coord.k < box_l.k-1.e-20) coord.k += domain.k;
#endif
}
//-----------------------------------------------------
// shift particle position according to direction
//-----------------------------------------------------
void Particle_base::Shift_coordinate_constrained(my_real domain, my_real box_l, my_real box_r, my_int dim)
{
#if P_DIM_X == 1 && PERI_DIM_X == 1
  if (dim.i != 0){
   if (coord.i > box_r.i+1.e-20) coord.i -= domain.i;
   else if (coord.i < box_l.i-1.e-20) coord.i += domain.i;
  } 
#endif
#if P_DIM_Y == 1 && PERI_DIM_Y == 1
  if (dim.j != 0){
    if (coord.j > box_r.j+1.e-20) coord.j -= domain.j;
    else if (coord.j < box_l.j-1.e-20) coord.j += domain.j;
  }
#endif
#if P_DIM_Z == 1 && PERI_DIM_Z == 1
  if (dim.k != 0){
    if (coord.k > box_r.k+1.e-20) coord.k -= domain.k;
    else if (coord.k < box_l.k-1.e-20) coord.k += domain.k;
  }
#endif
}
#endif
