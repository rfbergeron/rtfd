#include "scalar_solver.h"

#include <stdlib.h>

#define ROW_BORDER 1
#define COL_BORDER 1
#define IX(i, j) ((sim_size + 2 * COL_BORDER) * (i) + (j))
#define ACTUAL_SIZE ((sim_size + 2 * COL_BORDER) * (sim_size + 2 * ROW_BORDER))
#define COL_BEGIN COL_BORDER
#define COL_END (sim_size + COL_BORDER)
#define ROW_BEGIN ROW_BORDER
#define ROW_END (sim_size + ROW_BORDER)
#define SOLVER lin_solve
#define SOLVE(sim_size, type, x, x0, x1, a, c) \
  SOLVER(sim_size, type, x, x0, x1, a, c)
#define SWAP(x0, x)  \
  {                  \
    float* tmp = x0; \
    x0 = x;          \
    x = tmp;         \
  }

static void add_source(size_t sim_size, float* x, float* s, float dt,
                       float multiplier) {
  for (size_t i = 0, size = ACTUAL_SIZE; i < size; i++)
    x[i] += dt * multiplier * s[i];
}

// unfortunately, GCC is not smart enough to auto-vectorize this
/*
static void set_bnd_alt(size_t sim_size, MatrixType type, float* x) {
  switch (type) {
    case SLV_MAT_D:
    case SLV_MAT_D0:
      for (size_t j = COL_BEGIN; j < COL_END; ++j) {
        x[IX(ROW_BEGIN - 1, j)] = x[IX(ROW_BEGIN, j)];
        x[IX(ROW_END, j)] = x[IX(ROW_END - 1, j)];
      }
      for (size_t i = ROW_BEGIN; i < ROW_END; ++i) {
        x[IX(i, COL_BEGIN - 1)] = x[IX(i, COL_BEGIN)];
        x[IX(i, COL_END)] = x[IX(i, COL_END - 1)];
      }
      break;
    case SLV_MAT_U:
    case SLV_MAT_U0:
      for (size_t j = COL_BEGIN; j < COL_END; ++j) {
        x[IX(ROW_BEGIN - 1, j)] = -x[IX(ROW_BEGIN, j)];
        x[IX(ROW_END, j)] = -x[IX(ROW_END - 1, j)];
      }
      for (size_t i = ROW_BEGIN; i < ROW_END; ++i) {
        x[IX(i, COL_BEGIN - 1)] = x[IX(i, COL_BEGIN)];
        x[IX(i, COL_END)] = x[IX(i, COL_END - 1)];
      }
      break;
    case SLV_MAT_V:
    case SLV_MAT_V0:
      for (size_t j = COL_BEGIN; j < COL_END; ++j) {
        x[IX(ROW_BEGIN - 1, j)] = x[IX(ROW_BEGIN, j)];
        x[IX(ROW_END, j)] = x[IX(ROW_END - 1, j)];
      }
      for (size_t i = ROW_BEGIN; i < ROW_END; ++i) {
        x[IX(i, COL_BEGIN - 1)] = -x[IX(i, COL_BEGIN)];
        x[IX(i, COL_END)] = -x[IX(i, COL_END - 1)];
      }
      break;
    default:
      abort();
  }
  x[IX(ROW_BEGIN - 1, COL_BEGIN - 1)] =
      0.5f *
      (x[IX(ROW_BEGIN, COL_BEGIN - 1)] + x[IX(ROW_BEGIN - 1, COL_BEGIN)]);
  x[IX(ROW_BEGIN - 1, COL_END)] =
      0.5f * (x[IX(ROW_BEGIN, COL_END)] + x[IX(ROW_BEGIN - 1, COL_END - 1)]);
  x[IX(ROW_END, COL_BEGIN - 1)] =
      0.5f * (x[IX(ROW_END - 1, COL_BEGIN - 1)] + x[IX(ROW_END, COL_BEGIN)]);
  x[IX(ROW_END, COL_END)] =
      0.5f * (x[IX(ROW_END - 1, COL_END)] + x[IX(ROW_END, COL_END - 1)]);
}
*/

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
  x[IX(ROW_BEGIN - 1, COL_BEGIN - 1)] =
      0.5f *
      (x[IX(ROW_BEGIN, COL_BEGIN - 1)] + x[IX(ROW_BEGIN - 1, COL_BEGIN)]);
  x[IX(ROW_BEGIN - 1, COL_END)] =
      0.5f * (x[IX(ROW_BEGIN, COL_END)] + x[IX(ROW_BEGIN - 1, COL_END - 1)]);
  x[IX(ROW_END, COL_BEGIN - 1)] =
      0.5f * (x[IX(ROW_END - 1, COL_BEGIN - 1)] + x[IX(ROW_END, COL_BEGIN)]);
  x[IX(ROW_END, COL_END)] =
      0.5f * (x[IX(ROW_END - 1, COL_END)] + x[IX(ROW_END, COL_END - 1)]);
}

/*
static void jac_solve(size_t sim_size, MatrixType type, float* x, float* x0,
        float* x1, float a, float c) {
  for (size_t k = 0; k < 20; k++) {
    for (size_t i = ROW_BEGIN; i < ROW_END; i++) {
      for (size_t j = COL_BEGIN; j < COL_END; j++) {
        x1[IX(i, j)] =
            (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                 x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
            c;
      }
    }
    set_bnd(sim_size, type, x1);
    SWAP(x, x1);
  }
}
*/

static void lin_solve(size_t sim_size, MatrixType type, float* x, float* x0,
                      float* x1, float a, float c) {
  (void)x1;  // unused parameter
  for (size_t k = 0; k < 20; k++) {
    for (size_t i = ROW_BEGIN; i < ROW_END; i++) {
      for (size_t j = COL_BEGIN; j < COL_END; j++) {
        x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                           x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
                      c;
      }
    }
    set_bnd(sim_size, type, x);
  }
}

static void diffuse(size_t sim_size, MatrixType type, float* x, float* x0,
                    float* x1, float diff, float dt) {
  float a = dt * diff * sim_size * sim_size;
  SOLVE(sim_size, type, x, x0, x1, a, 1 + 4 * a);
}

// the compiler is smart enough to turn the if statements into calls to minss
// and maxss; there is little to no benefit to using them explicitly. you could
// use roundss and subss to get the fractional component of x and y, but all
// you'd be saving is one int-to-float conversion, which has roughly the same
// latency and throughput. we do however do that above because we should really
// be converting to 64-bit integers, which are too big to fit 4 of them into a
// single xmm register. there also isn't a way to convert packed floats to
// packed 64-bit integers until AVX512.
static void advect(size_t sim_size, MatrixType type, float* d, float* d0,
                   float* u, float* v, float dt) {
  float dt0 = dt * sim_size;
  for (size_t i = ROW_BEGIN; i < ROW_END; i++) {
    for (size_t j = COL_BEGIN; j < COL_END; j++) {
      float x = i - dt0 * u[IX(i, j)];
      float y = j - dt0 * v[IX(i, j)];
      if (x < ROW_BEGIN - 0.5f) x = ROW_BEGIN - 0.5f;
      if (x > ROW_END - 0.5f) x = ROW_END - 0.5f;
      size_t i0 = x;
      size_t i1 = i0 + 1;
      if (y < COL_BEGIN - 0.5f) y = COL_BEGIN - 0.5f;
      if (y > COL_END - 0.5f) y = COL_END - 0.5f;
      size_t j0 = y;
      size_t j1 = j0 + 1;
      float s1 = x - i0;
      float s0 = 1 - s1;
      float t1 = y - j0;
      float t0 = 1 - t1;
      d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                    s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
    }
  }
  set_bnd(sim_size, type, d);
}

static void project(size_t sim_size, float* u, float* v, float* p, float* div,
                    float* scratch) {
  for (size_t i = ROW_BEGIN; i < ROW_END; i++) {
    for (size_t j = COL_BEGIN; j < COL_END; j++) {
      div[IX(i, j)] = -0.5f *
                      (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] -
                       v[IX(i, j - 1)]) /
                      sim_size;
      p[IX(i, j)] = 0;
    }
  }
  set_bnd(sim_size, SLV_MAT_D, div);
  set_bnd(sim_size, SLV_MAT_D, p);

  SOLVE(sim_size, SLV_MAT_D, p, div, scratch, 1, 4);

  for (size_t i = ROW_BEGIN; i < ROW_END; i++) {
    for (size_t j = COL_BEGIN; j < COL_END; j++) {
      u[IX(i, j)] -= 0.5f * sim_size * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
      v[IX(i, j)] -= 0.5f * sim_size * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
    }
  }
  set_bnd(sim_size, SLV_MAT_U, u);
  set_bnd(sim_size, SLV_MAT_V, v);
}

void scalar_solver_dens_step(Solver* solver) {
  add_source(solver->sim_size, solver->d, solver->d_prev, solver->dt,
             solver->source);
  SWAP(solver->d_prev, solver->d);
  diffuse(solver->sim_size, SLV_MAT_D, solver->d, solver->d_prev,
          solver->u_prev, solver->diff, solver->dt);
  SWAP(solver->d_prev, solver->d);
  advect(solver->sim_size, SLV_MAT_D, solver->d, solver->d_prev, solver->u,
         solver->v, solver->dt);
}

void scalar_solver_vel_step(Solver* solver) {
  add_source(solver->sim_size, solver->u, solver->u_prev, solver->dt,
             solver->force);
  add_source(solver->sim_size, solver->v, solver->v_prev, solver->dt,
             solver->force);
  SWAP(solver->u_prev, solver->u);
  diffuse(solver->sim_size, SLV_MAT_U, solver->u, solver->u_prev,
          solver->d_prev, solver->visc, solver->dt);
  SWAP(solver->v_prev, solver->v);
  diffuse(solver->sim_size, SLV_MAT_V, solver->v, solver->v_prev,
          solver->d_prev, solver->visc, solver->dt);
  project(solver->sim_size, solver->u, solver->v, solver->u_prev,
          solver->v_prev, solver->d_prev);
  SWAP(solver->u_prev, solver->u);
  SWAP(solver->v_prev, solver->v);
  advect(solver->sim_size, SLV_MAT_U, solver->u, solver->u_prev, solver->u_prev,
         solver->v_prev, solver->dt);
  advect(solver->sim_size, SLV_MAT_V, solver->v, solver->v_prev, solver->u_prev,
         solver->v_prev, solver->dt);
  project(solver->sim_size, solver->u, solver->v, solver->u_prev,
          solver->v_prev, solver->d_prev);
}

float* scalar_solver_ix(Solver* solver, size_t x, size_t y, MatrixType type) {
  const size_t sim_size = solver->sim_size,
               offset = IX(x + ROW_BEGIN, y + COL_BEGIN);
  switch (type) {
    case SLV_MAT_D:
      return solver->d + offset;
    case SLV_MAT_D0:
      return solver->d_prev + offset;
    case SLV_MAT_U:
      return solver->u + offset;
    case SLV_MAT_U0:
      return solver->u_prev + offset;
    case SLV_MAT_V:
      return solver->v + offset;
    case SLV_MAT_V0:
      return solver->v_prev + offset;
    default:
      abort();
  }
}

void scalar_solver_clear(Solver* solver, MatrixType type) {
  const size_t sim_size = solver->sim_size;
  for (size_t i = 0, size = ACTUAL_SIZE; i < size; i++) {
    if (type & SLV_MAT_D) solver->d[i] = 0.0f;
    if (type & SLV_MAT_U) solver->u[i] = 0.0f;
    if (type & SLV_MAT_V) solver->v[i] = 0.0f;
    if (type & SLV_MAT_D0) solver->d_prev[i] = 0.0f;
    if (type & SLV_MAT_U0) solver->u_prev[i] = 0.0f;
    if (type & SLV_MAT_V0) solver->v_prev[i] = 0.0f;
  }
}

Solver* scalar_solver_init(size_t sim_size, float dt, float diff, float visc,
                           float force, float source) {
  size_t size = ACTUAL_SIZE;

  Solver* solver = malloc(sizeof(Solver));
  if (!solver) return NULL;

  solver->type = SLV_SCALAR;
  solver->sim_size = sim_size;
  solver->dt = dt;
  solver->diff = diff;
  solver->visc = visc;
  solver->force = force;
  solver->source = source;
  solver->u = malloc(size * sizeof(float));
  solver->v = malloc(size * sizeof(float));
  solver->u_prev = malloc(size * sizeof(float));
  solver->v_prev = malloc(size * sizeof(float));
  solver->d = malloc(size * sizeof(float));
  solver->d_prev = malloc(size * sizeof(float));

  if (!solver->u || !solver->v || !solver->u_prev || !solver->v_prev ||
      !solver->d || !solver->d_prev) {
    solver_destroy(solver);
    return NULL;
  }

  return solver;
}

void scalar_solver_destroy(Solver* solver) {
  free(solver->d);
  free(solver->d_prev);
  free(solver->u);
  free(solver->u_prev);
  free(solver->v);
  free(solver->v_prev);
  free(solver);
}
