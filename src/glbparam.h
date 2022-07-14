#ifndef GLBPARAM_H
#define GLBPARAM_H

#ifdef _INCPRS_
#define SOLVER SPH_incprs
#define Particle Particle_incprs
#define WATER 0
#define WALL 1
#endif
#ifdef _CPRS_
#define SOLVER SPH_cprs
#define Particle Particle_cprs
#endif
#ifdef _GSPH_
#define SOLVER SPH_gsph
#define Particle Particle_gsph
#endif
#ifdef _ALE_
#define SOLVER SPH_ale
#define Particle Particle_ale
#define CONDITION_NUMBER_DANGER 1000.
#endif
#ifdef _MESH_GENERATION_
#define SOLVER SPH_mesh_generation
#define Particle Particle_mesh_generation
#ifdef _DIM2_
  #define ICPX 4 //internal cells per package in X direction
  #define ICPY 4 //internal cells per package in Y direction
  #define ICPZ 1 //internal cells per package in Z direction
  #define BwidthX 4 //outer cells per package in X direction
  #define BwidthY 4 //outer cells per package in Y direction
  #define BwidthZ 0 //outer cells per package in Z direction
  #define TCPX 12 //total cells per package in X direction
  #define TCPY 12 //total cells per package in Y direction
  #define TCPZ 1 //total cells per package in Z direction
#endif
#ifdef _DIM3_
  #define ICPX 4 //internal cells per package in X direction
  #define ICPY 4 //internal cells per package in Y direction
  #define ICPZ 4 //internal cells per package in Z direction
  #define BwidthX 4 //outer cells per package in X direction
  #define BwidthY 4 //outer cells per package in Y direction
  #define BwidthZ 4 //outer cells per package in Z direction
  #define TCPX 12 //total cells per package in X direction
  #define TCPY 12 //total cells per package in Y direction
  #define TCPZ 12 //total cells per package in Z direction
#endif
#define NORMAL_CELL 0
#define CUT_CELL 1
#define SINGULARITY_CELL 1
#define SEGMENT_CELL 2
#define NARROW_BAND_SMALL 2
#define NARROW_BAND_LARGE 4
#define Emax 1
#define REAL_PARTICLE 0
#define SURFACE_PARTICLE 1
#define SINGULARITY_PARTICLE 2
#define SEGMENT_PARTICLE 3
#define GHOST_BC_PARTICLE 4
#define GHOST_BUFFER_PARTICLE 5
#define DUMMY_PARTICLE 6
#define NAGTIVE_PARTICLE -1
#define FILL 1
#define FULL 0
#define RELAX 2
#endif
#ifdef _EULARIAN_CAT_
#define SOLVER Eularian_cat
#define ICPX 4 //internal cells per package in X direction
#define ICPY 4 //internal cells per package in Y direction
#define ICPZ 1 //internal cells per package in Z direction
#define BwidthX 4 //outer cells per package in X direction
#define BwidthY 4 //outer cells per package in Y direction
#define BwidthZ 0 //outer cells per package in Z direction
#define TCPX 12 //total cells per package in X direction
#define TCPY 12 //total cells per package in Y direction
#define TCPZ 1 //total cells per package in Z direction
#define Emax 5
enum VARIABLES_cell_avg {MASS, MX, MY, MZ, ENERGY};
enum VARIABLES_primary {DEN, VX, VY, VZ, PRESSURE};
#endif

#ifdef _DIM1_
#define DIM 1
#define DIM_X 1
#define DIM_Y 0
#define DIM_Z 0
#endif
#ifdef _DIM2_
#define DIM 2
#define DIM_X 1
#define DIM_Y 1
#define DIM_Z 0
#endif
#ifdef _DIM3_
#define DIM 3
#define DIM_X 1
#define DIM_Y 1
#define DIM_Z 1
#endif

#ifdef _P_DIM1_
#define P_DIM 1
#define P_DIM_X 1
#define P_DIM_Y 0
#define P_DIM_Z 0
#endif
#ifdef _P_DIM2_
#define P_DIM 2
#define P_DIM_X 1
#define P_DIM_Y 1
#define P_DIM_Z 0
#endif
#ifdef _P_DIM3_
#define P_DIM 3
#define P_DIM_X 1
#define P_DIM_Y 1
#define P_DIM_Z 1
#endif

#ifdef _NON_PERI_
#define PERI_DIM 0
#define PERI_DIM_X 0
#define PERI_DIM_Y 0
#define PERI_DIM_Z 0
#endif
#ifdef _PERI_DIM1_
#define PERI_DIM 1
#define PERI_DIM_X 1
#define PERI_DIM_Y 0
#define PERI_DIM_Z 0
#endif
#ifdef _PERI_DIM2_
#define PERI_DIM 2
#define PERI_DIM_X 1
#define PERI_DIM_Y 1
#define PERI_DIM_Z 0
#endif
#ifdef _PERI_DIM3_
#define PERI_DIM 3
#define PERI_DIM_X 1
#define PERI_DIM_Y 1
#define PERI_DIM_Z 1
#endif

#ifdef _NON_SYM_
#define SYM_DIM 0
#define SYM_DIM_X 0
#define SYM_DIM_Y 0
#define SYM_DIM_Z 0
#endif
#ifdef _SYM_DIM1_
#define SYM_DIM 1
#define SYM_DIM_X 1
#define SYM_DIM_Y 0
#define SYM_DIM_Z 0
#endif
#ifdef _SYM_DIM2_
#define SYM_DIM 2
#define SYM_DIM_X 1
#define SYM_DIM_Y 1
#define SYM_DIM_Z 0
#endif
#ifdef _SYM_DIM3_
#define SYM_DIM 3
#define SYM_DIM_X 1
#define SYM_DIM_Y 1
#define SYM_DIM_Z 1
#endif

#if defined(_MESH_GENERATION_) && defined (_DIM3_)
  #define CELL_RATIO  1.45     // the ratio of cell size vs particle distance
  #define SCALE_RATIO 2      // the ratio between two consecutive level
  #define CUT_OFF     2.0
#elif defined(_MESH_GENERATION_) && defined (_DIM2_)
  #define CELL_RATIO  2.0     // the ratio of cell size vs particle distance
  #define SCALE_RATIO 2      // the ratio between two consecutive level
  #define CUT_OFF     2.0
#else
  #define CELL_RATIO  2.4     // the ratio of cell size vs particle distance
  #define SCALE_RATIO 2      // the ratio between two consecutive level
  #define CUT_OFF     2.0
#endif

#define PI 3.1415926
#define GRAVITY 1.

#endif
