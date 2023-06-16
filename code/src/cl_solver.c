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

static void set_bnd(size_t sim_size, MatrixType type, float* x) {
  for (size_t i = 0; i < sim_size; i++) {
    x[IX(ROW_BEGIN - 1, COL_BEGIN + i)] =
        (type == SLV_MAT_U || type == SLV_MAT_U0)
            ? -x[IX(ROW_BEGIN, COL_BEGIN + i)]
            : x[IX(ROW_BEGIN, COL_BEGIN + i)];
    x[IX(ROW_END, COL_BEGIN + i)] = (type == SLV_MAT_U || type == SLV_MAT_U0)
                                        ? -x[IX(ROW_END - 1, COL_BEGIN + i)]
                                        : x[IX(ROW_END - 1, COL_BEGIN + i)];
    x[IX(ROW_BEGIN + i, COL_BEGIN - 1)] =
        (type == SLV_MAT_V || type == SLV_MAT_V0)
            ? -x[IX(ROW_BEGIN + i, COL_BEGIN)]
            : x[IX(ROW_BEGIN + i, COL_BEGIN)];
    x[IX(ROW_BEGIN + i, COL_END)] = (type == SLV_MAT_V || type == SLV_MAT_V0)
                                        ? -x[IX(ROW_BEGIN + i, COL_END - 1)]
                                        : x[IX(ROW_BEGIN + i, COL_END - 1)];
  }
}

static void set_corners(size_t sim_size, MatrixType type, float* x) {
  if (type == SLV_MAT_D || type == SLV_MAT_D0) {
    x[IX(ROW_BEGIN - 1, COL_BEGIN - 1)] = x[IX(ROW_BEGIN, COL_BEGIN)];
    x[IX(ROW_BEGIN - 1, COL_END)] = x[IX(ROW_BEGIN, COL_END - 1)];
    x[IX(ROW_END, COL_BEGIN - 1)] = x[IX(ROW_END - 1, COL_BEGIN)];
    x[IX(ROW_END, COL_END)] = x[IX(ROW_END - 1, COL_END - 1)];
  } else {
    x[IX(ROW_BEGIN - 1, COL_BEGIN - 1)] = x[IX(ROW_BEGIN - 1, COL_END)] =
        x[IX(ROW_END, COL_BEGIN - 1)] = x[IX(ROW_END, COL_END)] = 0.0f;
  }
}

static void solve(size_t sim_size, MatrixType type, float* restrict x,
                  const float* restrict x0, float* restrict x1, float a,
                  float c, cl_bundle bundle) {
  (void)x1;  // unused parameter
  const char* errmsg;
  int status = cl_solve_setup(bundle, sim_size, x, x0, &errmsg);
  if (status) abort();
  for (size_t k = 0; k < 20; ++k) {
    status =
        cl_solve_step(bundle, sim_size, a, c,
                      (bool[2]){type == SLV_MAT_V, type == SLV_MAT_U}, &errmsg);
    if (status) abort();
  }
  status = cl_solve_retrieve(bundle, sim_size, x, &errmsg);
  if (status) abort();
}

static void diffuse(size_t sim_size, MatrixType type, float* restrict x,
                    const float* restrict x0, float* restrict x1, float diff,
                    float dt, cl_bundle bundle) {
  float a = dt * diff * sim_size * sim_size;
  solve(sim_size, type, x, x0, x1, a, 1 + 4 * a, bundle);
}

static void advect(size_t sim_size, MatrixType type, float* restrict d,
                   const float* restrict d0, const float* restrict u,
                   const float* restrict v, float dt, cl_bundle bundle) {
  const char* errmsg = NULL;
  cl_int status = cl_advect_setup(bundle, sim_size, d0, u, v, &errmsg);
  if (status) abort();
  status = cl_advect(bundle, sim_size, dt,
                     (bool[2]){type == SLV_MAT_V, type == SLV_MAT_U}, &errmsg);
  if (status) abort();
  status = cl_advect_retrieve(bundle, sim_size, d, &errmsg);
  if (status) abort();
}

static void project(size_t sim_size, float* restrict u, float* restrict v,
                    float* restrict p, float* restrict div,
                    float* restrict scratch, cl_bundle bundle) {
  (void)p, (void)div, (void)scratch;  // unused parameters
  const char* errmsg = NULL;
  int status = cl_project_setup(bundle, sim_size, u, v, &errmsg);
  if (status) abort();
  status = cl_project_one(bundle, sim_size, &errmsg);
  if (status) abort();
  for (size_t k = 0; k < 20; ++k) {
    status =
        cl_solve_step(bundle, sim_size, 1, 4, (bool[2]){false, false}, &errmsg);
    if (status) abort();
  }

  // setup again; result of solver will be preserved
  status = cl_project_setup(bundle, sim_size, u, v, &errmsg);
  if (status) abort();
  status = cl_project_two(bundle, sim_size, &errmsg);
  if (status) abort();
  status = cl_project_retrieve(bundle, sim_size, u, v, &errmsg);
  if (status) abort();
}

void cl_solver_dens_step(Solver* solver) {
  add_source(solver->sim_size, solver->h_buffers[0], solver->h_buffers[3],
             solver->dt, solver->source);
  SWAP(solver->h_buffers[3], solver->h_buffers[0]);
  diffuse(solver->sim_size, SLV_MAT_D, solver->h_buffers[0],
          solver->h_buffers[3], solver->h_buffers[4], solver->diff, solver->dt,
          solver->bundle);
  SWAP(solver->h_buffers[3], solver->h_buffers[0]);
  set_corners(solver->sim_size, SLV_MAT_D, solver->h_buffers[3]);
  advect(solver->sim_size, SLV_MAT_D, solver->h_buffers[0],
         solver->h_buffers[3], solver->h_buffers[1], solver->h_buffers[2],
         solver->dt, solver->bundle);
}

void cl_solver_vel_step(Solver* solver) {
  add_source(solver->sim_size, solver->h_buffers[1], solver->h_buffers[4],
             solver->dt, solver->force);
  add_source(solver->sim_size, solver->h_buffers[2], solver->h_buffers[5],
             solver->dt, solver->force);
  SWAP(solver->h_buffers[4], solver->h_buffers[1]);
  diffuse(solver->sim_size, SLV_MAT_U, solver->h_buffers[1],
          solver->h_buffers[4], solver->h_buffers[3], solver->visc, solver->dt,
          solver->bundle);
  SWAP(solver->h_buffers[5], solver->h_buffers[2]);
  diffuse(solver->sim_size, SLV_MAT_V, solver->h_buffers[2],
          solver->h_buffers[5], solver->h_buffers[3], solver->visc, solver->dt,
          solver->bundle);
  project(solver->sim_size, solver->h_buffers[1], solver->h_buffers[2],
          solver->h_buffers[4], solver->h_buffers[5], solver->h_buffers[3],
          solver->bundle);
  SWAP(solver->h_buffers[4], solver->h_buffers[1]);
  SWAP(solver->h_buffers[5], solver->h_buffers[2]);
  set_corners(solver->sim_size, SLV_MAT_U, solver->h_buffers[4]);
  advect(solver->sim_size, SLV_MAT_U, solver->h_buffers[1],
         solver->h_buffers[4], solver->h_buffers[4], solver->h_buffers[5],
         solver->dt, solver->bundle);
  set_corners(solver->sim_size, SLV_MAT_V, solver->h_buffers[5]);
  advect(solver->sim_size, SLV_MAT_V, solver->h_buffers[2],
         solver->h_buffers[5], solver->h_buffers[4], solver->h_buffers[5],
         solver->dt, solver->bundle);
  project(solver->sim_size, solver->h_buffers[1], solver->h_buffers[2],
          solver->h_buffers[4], solver->h_buffers[5], solver->h_buffers[3],
          solver->bundle);
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
