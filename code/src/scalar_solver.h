#ifndef SCALAR_SOLVER_H
#define SCALAR_SOLVER_H
#include "common_solver.h"
void scalar_solver_dens_step(Solver* solver);
void scalar_solver_vel_step(Solver* solver);
float* scalar_solver_ix(Solver* solver, size_t x, size_t y, MatrixType type);
void scalar_solver_clear(Solver* solver, MatrixType type);
Solver* scalar_solver_init(size_t sim_size, float dt, float diff, float visc,
                           float force, float source);
void scalar_solver_destroy(Solver* solver);
#endif
