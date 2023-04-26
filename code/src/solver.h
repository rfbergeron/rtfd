#ifndef SOLVER_H
#define SOLVER_H
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
  size_t sim_size, row_border, col_border;
  float *u, *v, *u_prev, *v_prev;
  float *d, *d_prev;
  float dt, diff, visc;
  float force, source;
} Solver;

Solver *solver_init(size_t sim_size, float dt, float diff, float visc,
                    float force, float source);
void solver_destroy(Solver *solver);
void solver_dens_step(Solver *solver);
void solver_vel_step(Solver *solver);
float *solver_at(Solver *solver, size_t x, size_t y, MatrixType type);
void solver_clear(Solver *solver, MatrixType type);
#endif
