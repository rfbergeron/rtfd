
#ifndef COMMON_SOLVER_H
#define COMMON_SOLVER_H
#include <stddef.h>

typedef enum matrix_type {
  SLV_MAT_D = 1 << 0,
  SLV_MAT_U = 1 << 1,
  SLV_MAT_V = 1 << 2,
  SLV_MAT_D0 = 1 << 3,
  SLV_MAT_U0 = 1 << 4,
  SLV_MAT_V0 = 1 << 5,
  SLV_MAT_CURR = SLV_MAT_D | SLV_MAT_U | SLV_MAT_V,
  SLV_MAT_PREV = SLV_MAT_D0 | SLV_MAT_U0 | SLV_MAT_V0,
  SLV_MAT_ALL = SLV_MAT_CURR | SLV_MAT_PREV
} MatrixType;

typedef struct solver {
  size_t sim_size;
  float *u, *v, *u_prev, *v_prev;
  float *d, *d_prev;
  void (*vel_fn)(struct solver *), (*dens_fn)(struct solver *),
      (*clear_fn)(struct solver *, MatrixType);
  float *(*ix_fn)(struct solver *, size_t, size_t, MatrixType);
  float dt, diff, visc;
  float force, source;
} Solver;

void solver_destroy(Solver *solver);
void solver_dens_step(Solver *solver);
void solver_vel_step(Solver *solver);
float *solver_ix(Solver *solver, size_t x, size_t y, MatrixType type);
void solver_clear(Solver *solver, MatrixType type);
#endif
