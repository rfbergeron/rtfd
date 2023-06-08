#include "common_solver.h"

#include <stdlib.h>

#include "cl_solver.h"
#include "scalar_solver.h"
#include "sse4_2_solver.h"

const size_t MAT_TYPE_IXS[] = {
    [SLV_MAT_D] = 0,  [SLV_MAT_U] = 1,  [SLV_MAT_V] = 2,
    [SLV_MAT_D0] = 3, [SLV_MAT_U0] = 4, [SLV_MAT_V0] = 5};

void solver_dens_step(Solver* solver) {
  switch (solver->type) {
    case SLV_SCALAR:
      return scalar_solver_dens_step(solver);
    case SLV_SSE4_2:
      return sse4_2_solver_dens_step(solver);
    case SLV_CL:
      return cl_solver_dens_step(solver);
    default:
      abort();
  }
}

void solver_vel_step(Solver* solver) {
  switch (solver->type) {
    case SLV_SCALAR:
      return scalar_solver_vel_step(solver);
    case SLV_SSE4_2:
      return sse4_2_solver_vel_step(solver);
    case SLV_CL:
      return cl_solver_vel_step(solver);
    default:
      abort();
  }
}

float* solver_ix(Solver* solver, size_t x, size_t y, MatrixType type) {
  switch (solver->type) {
    case SLV_SCALAR:
      return scalar_solver_ix(solver, x, y, type);
    case SLV_SSE4_2:
      return sse4_2_solver_ix(solver, x, y, type);
    case SLV_CL:
      return cl_solver_ix(solver, x, y, type);
    default:
      abort();
  }
}

void solver_clear(Solver* solver, MatrixType type) {
  switch (solver->type) {
    case SLV_SCALAR:
      return scalar_solver_clear(solver, type);
    case SLV_SSE4_2:
      return sse4_2_solver_clear(solver, type);
    case SLV_CL:
      return cl_solver_clear(solver, type);
    default:
      abort();
  }
}

Solver* solver_init(size_t sim_size, SolverType type, float dt, float diff,
                    float visc, float force, float source) {
  switch (type) {
    case SLV_SCALAR:
      return scalar_solver_init(sim_size, dt, diff, visc, force, source);
    case SLV_SSE4_2:
      return sse4_2_solver_init(sim_size, dt, diff, visc, force, source);
    case SLV_CL:
      return cl_solver_init(sim_size, dt, diff, visc, force, source);
    default:
      abort();
  }
}

void solver_destroy(Solver* solver) {
  switch (solver->type) {
    case SLV_SCALAR:
      return scalar_solver_destroy(solver);
    case SLV_SSE4_2:
      return sse4_2_solver_destroy(solver);
    case SLV_CL:
      return cl_solver_destroy(solver);
    default:
      abort();
  }
}
