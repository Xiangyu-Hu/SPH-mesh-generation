#ifndef PARALLELFUNC_MESH_H
#define PARALLELFUNC_MESH_H

#include "glbcls.h"
#include "float.h"
#ifdef _MESH_GENERATION_
#include "level_set.h"
#include "lset_cell.h"
#include "lset_package.h"
#endif

using namespace std;
using namespace tbb;

void Parallel_get_global_scale(int n, std::vector<p_Levelset_package> &lset_pkg, Levelset *level_set);

void Parallel_get_psi_max(int n, std::vector<p_Levelset_package> &interface_pkg, Levelset *level_set);

void Parallel_get_total_volume_mass(int n, std::vector<p_Levelset_package> &lset_pkg, Levelset *level_set);

void Parallel_get_num_of_interface_cell(int n, std::vector<p_Levelset_package> &interface_pkg, Levelset *level_set, int &num_interface_cell);

#endif