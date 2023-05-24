#ifndef SCALAR_SOLVER_H
#define SCALAR_SOLVER_H
#include "common_solver.h"
Solver *scalar_solver_init(size_t sim_size, float dt, float diff, float visc,
                           float force, float source);
#endif
