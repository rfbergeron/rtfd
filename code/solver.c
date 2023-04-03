#include "solver.h"

#include <stdlib.h>

#define SOLVER jac_solve
#define SOLVE(N, type, x, x0, a, c) SOLVER(N, type, x, x0, a, c)

typedef enum matrix_type { MAT_FLUID, MAT_U_VEL, MAT_V_VEL } MatrixType;

#define SWAP(x0, x)  \
  {                  \
    float* tmp = x0; \
    x0 = x;          \
    x = tmp;         \
  }

void add_source(size_t N, float* x, float* s, float dt) {
  for (size_t i = 0, size = ACTUALSIZE; i < size; i++) x[i] += dt * s[i];
}

void set_bnd(size_t N, MatrixType type, float* x) {
  for (size_t i = 1; i <= N; i++) {
    x[IX(0, i)] = type == MAT_U_VEL ? -x[IX(1, i)] : x[IX(1, i)];
    x[IX(N + 1, i)] = type == MAT_U_VEL ? -x[IX(N, i)] : x[IX(N, i)];
    x[IX(i, 0)] = type == MAT_V_VEL ? -x[IX(i, 1)] : x[IX(i, 1)];
    x[IX(i, N + 1)] = type == MAT_V_VEL ? -x[IX(i, N)] : x[IX(i, N)];
  }
  x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
  x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
  x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
  x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

void jac_solve(size_t N, MatrixType type, float* x, float* x0, float a,
               float c) {
  float* x1 = malloc(ACTUALSIZE * sizeof(float));

  for (size_t k = 0; k < 20; k++) {
    for (size_t i = 1; i <= N; i++) {
      for (size_t j = 1; j <= N; j++) {
        x1[IX(i, j)] =
            (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                 x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
            c;
      }
    }
    set_bnd(N, type, x1);
    SWAP(x, x1);
  }

  free(x1);
}

void lin_solve(size_t N, MatrixType type, float* x, float* x0, float a,
               float c) {
  for (size_t k = 0; k < 20; k++) {
    for (size_t i = 1; i <= N; i++) {
      for (size_t j = 1; j <= N; j++) {
        x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                           x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
                      c;
      }
    }
    set_bnd(N, type, x);
  }
}

void diffuse(size_t N, MatrixType type, float* x, float* x0, float diff,
             float dt) {
  float a = dt * diff * N * N;
  SOLVE(N, type, x, x0, a, 1 + 4 * a);
}

void advect(size_t N, MatrixType type, float* d, float* d0, float* u, float* v,
            float dt) {
  float dt0 = dt * N;
  for (size_t i = 1; i <= N; i++) {
    for (size_t j = 1; j <= N; j++) {
      float x = i - dt0 * u[IX(i, j)];
      float y = j - dt0 * v[IX(i, j)];
      if (x < 0.5f) x = 0.5f;
      if (x > N + 0.5f) x = N + 0.5f;
      size_t i0 = x;
      size_t i1 = i0 + 1;
      if (y < 0.5f) y = 0.5f;
      if (y > N + 0.5f) y = N + 0.5f;
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
  set_bnd(N, type, d);
}

void project(size_t N, float* u, float* v, float* p, float* div) {
  for (size_t i = 1; i <= N; i++) {
    for (size_t j = 1; j <= N; j++) {
      div[IX(i, j)] = -0.5f *
                      (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] -
                       v[IX(i, j - 1)]) /
                      N;
      p[IX(i, j)] = 0;
    }
  }
  set_bnd(N, MAT_FLUID, div);
  set_bnd(N, MAT_FLUID, p);

  SOLVE(N, MAT_FLUID, p, div, 1, 4);

  for (size_t i = 1; i <= N; i++) {
    for (size_t j = 1; j <= N; j++) {
      u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
      v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
    }
  }
  set_bnd(N, MAT_U_VEL, u);
  set_bnd(N, MAT_V_VEL, v);
}

void dens_step(size_t N, float* x, float* x0, float* u, float* v, float diff,
               float dt) {
  add_source(N, x, x0, dt);
  SWAP(x0, x);
  diffuse(N, MAT_FLUID, x, x0, diff, dt);
  SWAP(x0, x);
  advect(N, MAT_FLUID, x, x0, u, v, dt);
}

void vel_step(size_t N, float* u, float* v, float* u0, float* v0, float visc,
              float dt) {
  add_source(N, u, u0, dt);
  add_source(N, v, v0, dt);
  SWAP(u0, u);
  diffuse(N, MAT_U_VEL, u, u0, visc, dt);
  SWAP(v0, v);
  diffuse(N, MAT_V_VEL, v, v0, visc, dt);
  project(N, u, v, u0, v0);
  SWAP(u0, u);
  SWAP(v0, v);
  advect(N, MAT_U_VEL, u, u0, u0, v0, dt);
  advect(N, MAT_V_VEL, v, v0, u0, v0, dt);
  project(N, u, v, u0, v0);
}
