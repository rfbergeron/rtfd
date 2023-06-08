#ifndef CL_SOLVER_H
#define CL_SOLVER_H
#include "common_solver.h"
void cl_solver_dens_step(Solver* solver);
void cl_solver_vel_step(Solver* solver);
float* cl_solver_ix(Solver* solver, size_t x, size_t y, MatrixType type);
void cl_solver_clear(Solver* solver, MatrixType type);
Solver* cl_solver_init(size_t sim_size, float dt, float diff, float visc,
                       float force, float source);
void cl_solver_destroy(Solver* solver);
#endif
