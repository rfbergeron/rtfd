#include "cl_solver.h"

#include <stdio.h>
#include <stdlib.h>

#include "cl_common.h"
#include "cl_helper.h"

#define SWAP(x0, x)  \
  {                  \
    float* tmp = x0; \
    x0 = x;          \
    x = tmp;         \
  }

static void add_source(size_t sim_size, float* restrict x,
                       const float* restrict s, float dt, float multiplier) {
  for (size_t i = 0, size = ACTUAL_SIZE; i < size; i++)
    x[i] += dt * multiplier * s[i];
}

void cl_solver_dens_step(Solver* solver) {
  add_source(solver->sim_size, solver->h_buffers[0], solver->h_buffers[3],
             solver->dt, solver->source);
  SWAP(solver->h_buffers[3], solver->h_buffers[0]);
  const char* errmsg;
  cl_int status = cl_dens_step_full(solver->bundle, solver->sim_size,
                                    solver->h_buffers[0], solver->h_buffers[3],
                                    solver->h_buffers[1], solver->h_buffers[2],
                                    solver->diff, solver->dt, 20, &errmsg);
  if (status) abort();
}

void cl_solver_vel_step(Solver* solver) {
  add_source(solver->sim_size, solver->h_buffers[1], solver->h_buffers[4],
             solver->dt, solver->force);
  add_source(solver->sim_size, solver->h_buffers[2], solver->h_buffers[5],
             solver->dt, solver->force);
  SWAP(solver->h_buffers[4], solver->h_buffers[1]);
  SWAP(solver->h_buffers[5], solver->h_buffers[2]);
  const char* errmsg;
  cl_int status = cl_vel_step_full(solver->bundle, solver->sim_size,
                                   solver->h_buffers[1], solver->h_buffers[2],
                                   solver->h_buffers[4], solver->h_buffers[5],
                                   solver->visc, solver->dt, 20, &errmsg);
  if (status) abort();
}

float* cl_solver_ix(Solver* solver, size_t x, size_t y, MatrixType type) {
  const size_t sim_size = solver->sim_size,
               offset = IX(x + ROW_BEGIN, y + COL_BEGIN);
  switch (type) {
    case SLV_MAT_D:
      return solver->h_buffers[0] + offset;
    case SLV_MAT_U:
      return solver->h_buffers[1] + offset;
    case SLV_MAT_V:
      return solver->h_buffers[2] + offset;
    case SLV_MAT_D0:
      return solver->h_buffers[3] + offset;
    case SLV_MAT_U0:
      return solver->h_buffers[4] + offset;
    case SLV_MAT_V0:
      return solver->h_buffers[5] + offset;
    default:
      abort();
  }
}

void cl_solver_clear(Solver* solver, MatrixType type) {
  const size_t sim_size = solver->sim_size;
  for (size_t i = 0, size = ACTUAL_SIZE; i < size; i++) {
    if (type & SLV_MAT_D) solver->h_buffers[0][i] = 0.0f;
    if (type & SLV_MAT_U) solver->h_buffers[1][i] = 0.0f;
    if (type & SLV_MAT_V) solver->h_buffers[2][i] = 0.0f;
    if (type & SLV_MAT_D0) solver->h_buffers[3][i] = 0.0f;
    if (type & SLV_MAT_U0) solver->h_buffers[4][i] = 0.0f;
    if (type & SLV_MAT_V0) solver->h_buffers[5][i] = 0.0f;
  }
}

Solver* cl_solver_init(size_t sim_size, float dt, float diff, float visc,
                       float force, float source) {
  if (sim_size < 256 || sim_size % 256 != 0) abort();
  size_t size = ACTUAL_SIZE;

  Solver* solver = malloc(sizeof(Solver));
  if (!solver) return NULL;

  solver->type = SLV_CL;
  solver->sim_size = sim_size;
  solver->dt = dt;
  solver->diff = diff;
  solver->visc = visc;
  solver->force = force;
  solver->source = source;
  for (size_t i = 0; i < 6; ++i)
    solver->h_buffers[i] = malloc(size * sizeof(float));
  for (size_t i = 0; i < 6; ++i)
    if (!solver->h_buffers[i]) return solver_destroy(solver), NULL;

  cl_int status;
  const char* errmsg_out;
  solver->bundle = init_gpu_bundle(sim_size, &status, &errmsg_out);
  if (status) abort();

  return solver;
}

void cl_solver_destroy(Solver* solver) {
  for (size_t i = 0; i < 6; ++i) free(solver->h_buffers[i]);
  free_bundle(solver->bundle);
  free(solver);
}
