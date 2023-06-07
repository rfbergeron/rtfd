#ifndef SSE4_2_SOLVER_H
#define SSE4_2_SOLVER_H
#include "common_solver.h"
void sse4_2_solver_dens_step(Solver* solver);
void sse4_2_solver_vel_step(Solver* solver);
float* sse4_2_solver_ix(Solver* solver, size_t x, size_t y, MatrixType type);
void sse4_2_solver_clear(Solver* solver, MatrixType type);
Solver* sse4_2_solver_init(size_t sim_size, float dt, float diff, float visc,
                           float force, float source);
void sse4_2_solver_destroy(Solver* solver);
#endif
