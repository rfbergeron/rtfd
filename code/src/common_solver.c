#include "common_solver.h"

#include <stdlib.h>

void solver_dens_step(Solver* solver) { return solver->dens_fn(solver); }

void solver_vel_step(Solver* solver) { return solver->vel_fn(solver); }

float* solver_ix(Solver* solver, size_t x, size_t y, MatrixType type) {
  return solver->ix_fn(solver, x, y, type);
}

void solver_clear(Solver* solver, MatrixType type) {
  return solver->clear_fn(solver, type);
}

void solver_destroy(Solver* solver) {
  free(solver->d);
  free(solver->d_prev);
  free(solver->u);
  free(solver->u_prev);
  free(solver->v);
  free(solver->v_prev);
  free(solver);
}
