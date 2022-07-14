#ifndef _BOUNDARY_SPECIFICATIONS_H
#define _BOUNDARY_SPECIFICATIONS_H

//All kinds of Boundary Conditions
enum BC_type {
	Interior,
	Inflow,
	Outflow,
	Symmetry,
	Periodic,
	Wall
};
#endif