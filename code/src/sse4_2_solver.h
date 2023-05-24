#ifndef SSE4_2_SOLVER_H
#define SSE4_2_SOLVER_H
#include "common_solver.h"
Solver *sse4_2_solver_init(size_t sim_size, float dt, float diff, float visc,
                           float force, float source);
#endif
