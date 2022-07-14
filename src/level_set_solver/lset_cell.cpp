#include "glbfunc.h"
#include "lset_cell.h"
#include "level_set.h"
#include "lset_package.h"

/***************************************************/
/*                                                 */
/*    Functions defined in class "Levelset_cell"   */
/*                                                 */
/***************************************************/
//------------------------------------------------------------
// initialze all the necessary parameters 
// for Levelset_cell
//------------------------------------------------------------
void Levelset_cell::Initialize(Levelset_package *lset_pkg)
{
  phi = 0.;
  phi_1 = 0.;
  n_x = 0.;
  n_y = 0.;
  n_z = 0.;
  curv = 0.;
  curv_1 = 0.;
  psi = 1.e-20;
  psi_1 = 0.;
  scale = 0.;
  vol = 0.;
  area = 0.;
  
  for (int i=0; i<Emax; i++){
    LU[i] = 0.;
  }
  tag_interface = 0;
  tag_characteristic = 0;

  idx_characteristic = -1;
}
//------------------------------------------------------------
// set phi
//------------------------------------------------------------
void Levelset_cell::Set_phi (Real phi_){
  phi = phi_;
}
//------------------------------------------------------------
// set phi_1
//------------------------------------------------------------
void Levelset_cell::Set_phi_1 (){
  phi_1 = phi;
}
//------------------------------------------------------------
// calculate normals
//------------------------------------------------------------
void Levelset_cell::Get_normal (){
  
}
//------------------------------------------------------------
// Update_Extend_psi
//------------------------------------------------------------
void Levelset_cell::Update_extend_psi (){
  psi += LU[0];
}
//------------------------------------------------------------
// calculate normals
//------------------------------------------------------------
void Levelset_cell::Update_extend_curv (){
  curv += LU[0];
}
