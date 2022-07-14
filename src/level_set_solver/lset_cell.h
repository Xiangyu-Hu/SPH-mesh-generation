#ifndef _LEVEL_SET_CELL_
#define _LEVEL_SET_CELL_
#include "glbcls.h"
#include "boundary_condition.h"

class   Levelset_cell{

private:

public: 
  // variable
  Real phi, phi_1;        // for generating the unstructured mesh  
  Real n_x, n_y, n_z;     // normal direction
  Real curv, curv_1;      // measure the curvature
  Real psi, psi_1;        // for the distance function on the interface 
  Real scale;
  Real vol;
  Real area;
  Real LU[Emax];          // cell average values

  int tag_interface;
  int tag_characteristic;

  int idx_characteristic;

  Levelset_cell(){};
  // function
  void Initialize(Levelset_package *lset_pkg);
  void Set_phi (Real phi_);
  void Set_phi_1 ();
  void Get_normal ();
  void Update_extend_psi ();
  void Update_extend_curv ();
};

#endif