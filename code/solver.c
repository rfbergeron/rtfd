#include "solver.h"

#include <immintrin.h>
#include <stdlib.h>
#include <string.h>

#define SOLVER sse2_solve
#define SOLVE(N, type, x, x0, a, c) SOLVER(N, type, x, x0, a, c)
#define PROJECTOR sse2_project
#define PROJECT(N, u, v, p, div) PROJECTOR(N, u, v, p, div)
#define ADDER sse2_add_source
#define ADD_SOURCE(N, x, s, dt) ADDER(N, x, s, dt)
#define BOUNDER sse2_set_bnd
#define SET_BND(N, type, x) BOUNDER(N, type, x)

typedef enum matrix_type { MAT_FLUID, MAT_U_VEL, MAT_V_VEL } MatrixType;

#define SWAP(x0, x)  \
  {                  \
    float* tmp = x0; \
    x0 = x;          \
    x = tmp;         \
  }

void sse2_add_source(size_t N, float* x, float* s, float dt) {
  __m128 dt_vec = _mm_set1_ps(dt);
  for (size_t i = 0, size = ACTUALSIZE; i < size; i += 4) {
    __m128 x_current = _mm_load_ps(x + i), s_current = _mm_load_ps(s + i);
    s_current = _mm_mul_ps(s_current, dt_vec);
    x_current = _mm_add_ps(x_current, s_current);
    _mm_store_ps(x + i, x_current);
  }
}

void add_source(size_t N, float* x, float* s, float dt) {
  for (size_t i = 0, size = ACTUALSIZE; i < size; i++) x[i] += dt * s[i];
}

void sse2_set_bnd(size_t N, MatrixType type, float* x) {
  switch (type) {
    case MAT_FLUID:
      for (size_t j = 0; j < N; j += 4) {
        _mm_store_ps(x + IX(ROWBEGIN - 1, COLBEGIN + j),
                     _mm_load_ps(x + IX(ROWBEGIN, COLBEGIN + j)));
        _mm_store_ps(x + IX(ROWEND, COLBEGIN + j),
                     _mm_load_ps(x + IX(ROWEND - 1, COLBEGIN + j)));
        for (size_t i = j; i < j + 4; ++i) {
          x[IX(ROWBEGIN + i, COLBEGIN - 1)] = x[IX(ROWBEGIN + i, COLBEGIN)];
          x[IX(ROWBEGIN + i, COLEND)] = x[IX(ROWBEGIN + i, COLEND - 1)];
        }
      }
      break;
    case MAT_U_VEL:
      for (size_t j = 0; j < N; j += 4) {
        _mm_store_ps(x + IX(ROWBEGIN - 1, COLBEGIN + j),
                     _mm_sub_ps(_mm_setzero_ps(),
                                _mm_load_ps(x + IX(ROWBEGIN, COLBEGIN + j))));
        _mm_store_ps(x + IX(ROWEND, COLBEGIN + j),
                     _mm_sub_ps(_mm_setzero_ps(),
                                _mm_load_ps(x + IX(ROWEND - 1, COLBEGIN + j))));
        for (size_t i = j; i < j + 4; ++i) {
          x[IX(ROWBEGIN + i, COLBEGIN - 1)] = x[IX(ROWBEGIN + i, COLBEGIN)];
          x[IX(ROWBEGIN + i, COLEND)] = x[IX(ROWBEGIN + i, COLEND - 1)];
        }
      }
      break;
    case MAT_V_VEL:
      for (size_t j = 0; j < N; j += 4) {
        _mm_store_ps(x + IX(ROWBEGIN - 1, COLBEGIN + j),
                     _mm_load_ps(x + IX(ROWBEGIN, COLBEGIN + j)));
        _mm_store_ps(x + IX(ROWEND, COLBEGIN + j),
                     _mm_load_ps(x + IX(ROWEND - 1, COLBEGIN + j)));
        for (size_t i = j; i < j + 4; ++i) {
          x[IX(ROWBEGIN + i, COLBEGIN - 1)] = -x[IX(ROWBEGIN + i, COLBEGIN)];
          x[IX(ROWBEGIN + i, COLEND)] = -x[IX(ROWBEGIN + i, COLEND - 1)];
        }
      }
      break;
    default:
      abort();
  }
  x[IX(ROWBEGIN - 1, COLBEGIN - 1)] =
      0.5f * (x[IX(ROWBEGIN, COLBEGIN - 1)] + x[IX(ROWBEGIN - 1, COLBEGIN)]);
  x[IX(ROWBEGIN - 1, COLEND)] =
      0.5f * (x[IX(ROWBEGIN, COLEND)] + x[IX(ROWBEGIN - 1, COLEND - 1)]);
  x[IX(ROWEND, COLBEGIN - 1)] =
      0.5f * (x[IX(ROWEND - 1, COLBEGIN - 1)] + x[IX(ROWEND, COLBEGIN)]);
  x[IX(ROWEND, COLEND)] =
      0.5f * (x[IX(ROWEND - 1, COLEND)] + x[IX(ROWEND, COLEND - 1)]);
}

void set_bnd(size_t N, MatrixType type, float* x) {
  for (size_t i = 0; i < N; i++) {
    x[IX(ROWBEGIN - 1, COLBEGIN + i)] = type == MAT_U_VEL
                                            ? -x[IX(ROWBEGIN, COLBEGIN + i)]
                                            : x[IX(ROWBEGIN, COLBEGIN + i)];
    x[IX(ROWEND, COLBEGIN + i)] = type == MAT_U_VEL
                                      ? -x[IX(ROWEND - 1, COLBEGIN + i)]
                                      : x[IX(ROWEND - 1, COLBEGIN + i)];
    x[IX(ROWBEGIN + i, COLBEGIN - 1)] = type == MAT_V_VEL
                                            ? -x[IX(ROWBEGIN + i, COLBEGIN)]
                                            : x[IX(ROWBEGIN + i, COLBEGIN)];
    x[IX(ROWBEGIN + i, COLEND)] = type == MAT_V_VEL
                                      ? -x[IX(ROWBEGIN + i, COLEND - 1)]
                                      : x[IX(ROWBEGIN + i, COLEND - 1)];
  }
  x[IX(ROWBEGIN - 1, COLBEGIN - 1)] =
      0.5f * (x[IX(ROWBEGIN, COLBEGIN - 1)] + x[IX(ROWBEGIN - 1, COLBEGIN)]);
  x[IX(ROWBEGIN - 1, COLEND)] =
      0.5f * (x[IX(ROWBEGIN, COLEND)] + x[IX(ROWBEGIN - 1, COLEND - 1)]);
  x[IX(ROWEND, COLBEGIN - 1)] =
      0.5f * (x[IX(ROWEND - 1, COLBEGIN - 1)] + x[IX(ROWEND, COLBEGIN)]);
  x[IX(ROWEND, COLEND)] =
      0.5f * (x[IX(ROWEND - 1, COLEND)] + x[IX(ROWEND, COLEND - 1)]);
}

void sse2_solve(size_t N, MatrixType type, float* x, float* x0, float a,
                float c) {
  float* x1 = aligned_alloc(64, ACTUALSIZE * sizeof(float));
  __m128 c_inv_vec = _mm_set1_ps(1.0f / c), a_vec = _mm_set1_ps(a);

  for (size_t k = 0; k < 20; k++) {
    for (size_t i = ROWBEGIN; i < ROWEND; ++i) {
      for (size_t j = COLBEGIN; j < COLEND; j += 4) {
        __m128 above = _mm_load_ps(x + IX(i + 1, j));
        __m128 current = _mm_load_ps(x + IX(i, j));
        __m128 below = _mm_load_ps(x + IX(i - 1, j));
        __m128 x0_vec = _mm_load_ps(x0 + IX(i, j));

        /* row-by-row addition */
        __m128 dest = _mm_add_ps(above, below);

        /* column-by-column addition */
        __m128 shiftr =
            _mm_shuffle_ps(current, current, _MM_SHUFFLE(0, 0, 1, 2));
        __m128 tempr = _mm_load_ss(x + IX(i, j - 1));
        shiftr = _mm_move_ss(shiftr, tempr);
        dest = _mm_add_ps(dest, shiftr);

        /* shuffle after since we don't need to save slot zero but we do need
         * the value we just loaded to be in the high slot
         */
        __m128 templ = _mm_load_ss(x + IX(i, j + 4));
        __m128 shiftl = _mm_move_ss(current, templ);
        shiftl = _mm_shuffle_ps(shiftl, shiftl, _MM_SHUFFLE(1, 2, 3, 0));
        dest = _mm_add_ps(dest, shiftl);

        /* multiply by a; add x0; divide by c */
        dest = _mm_mul_ps(dest, a_vec);
        dest = _mm_add_ps(dest, x0_vec);
        dest = _mm_mul_ps(dest, c_inv_vec);

        /* send it back */
        _mm_store_ps(x1 + IX(i, j), dest);
      }
    }
    SET_BND(N, type, x1);
    SWAP(x, x1);
  }
  free(x1);
}

void jac_solve(size_t N, MatrixType type, float* x, float* x0, float a,
               float c) {
  float* x1 = aligned_alloc(64, ACTUALSIZE * sizeof(float));

  for (size_t k = 0; k < 20; k++) {
    for (size_t i = ROWBEGIN; i < ROWEND; i++) {
      for (size_t j = COLBEGIN; j < COLEND; j++) {
        x1[IX(i, j)] =
            (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                 x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
            c;
      }
    }
    SET_BND(N, type, x1);
    SWAP(x, x1);
  }

  free(x1);
}

void lin_solve(size_t N, MatrixType type, float* x, float* x0, float a,
               float c) {
  for (size_t k = 0; k < 20; k++) {
    for (size_t i = ROWBEGIN; i < ROWEND; i++) {
      for (size_t j = COLBEGIN; j < COLEND; j++) {
        x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                           x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
                      c;
      }
    }
    SET_BND(N, type, x);
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
  for (size_t i = ROWBEGIN; i < ROWEND; i++) {
    for (size_t j = COLBEGIN; j < COLEND; j++) {
      float x = i - dt0 * u[IX(i, j)];
      float y = j - dt0 * v[IX(i, j)];
      if (x < ROWBEGIN - 0.5f) x = ROWBEGIN - 0.5f;
      if (x > ROWEND - 0.5f) x = ROWEND - 0.5f;
      size_t i0 = x;
      size_t i1 = i0 + 1;
      if (y < COLBEGIN - 0.5f) y = COLBEGIN - 0.5f;
      if (y > COLEND - 0.5f) y = COLEND - 0.5f;
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
  SET_BND(N, type, d);
}

void sse2_project(size_t N, float* u, float* v, float* p, float* div) {
  __m128 multiplier = _mm_set1_ps(-0.5f / N);
  memset(p, 0, ACTUALSIZE * sizeof(float));
  for (size_t i = ROWBEGIN; i < ROWEND; i++) {
    for (size_t j = COLBEGIN; j < COLEND; j += 4) {
      __m128 u_above = _mm_load_ps(u + IX(i + 1, j));
      __m128 u_below = _mm_load_ps(u + IX(i - 1, j));
      __m128 u_div = _mm_sub_ps(u_above, u_below);

      __m128 v_right = _mm_load_ps(v + IX(i, j));
      __m128 v_right_temp = _mm_load_ss(v + IX(i, j + 4));
      v_right = _mm_move_ss(v_right, v_right_temp);
      v_right = _mm_shuffle_ps(v_right, v_right, _MM_SHUFFLE(0, 3, 2, 1));

      __m128 v_left = _mm_load_ps(v + IX(i, j));
      v_left = _mm_shuffle_ps(v_left, v_left, _MM_SHUFFLE(2, 1, 0, 0));
      __m128 v_left_temp = _mm_load_ss(v + IX(i, j - 1));
      v_left = _mm_move_ss(v_left, v_left_temp);

      __m128 v_div = _mm_sub_ps(v_right, v_left);

      __m128 result = _mm_add_ps(u_div, v_div);
      result = _mm_mul_ps(result, multiplier);
      _mm_store_ps(div + IX(i, j), result);
    }
  }
  SET_BND(N, MAT_FLUID, div);

  SOLVE(N, MAT_FLUID, p, div, 1, 4);

  multiplier = _mm_set1_ps(0.5f * N);
  for (size_t i = ROWBEGIN; i < ROWEND; i++) {
    for (size_t j = COLBEGIN; j < COLEND; j += 4) {
      __m128 p_above = _mm_load_ps(p + IX(i + 1, j));
      __m128 p_below = _mm_load_ps(p + IX(i - 1, j));
      __m128 p_vert_diff = _mm_sub_ps(p_above, p_below);
      p_vert_diff = _mm_mul_ps(p_vert_diff, multiplier);

      __m128 u_current = _mm_load_ps(u + IX(i, j));
      u_current = _mm_sub_ps(u_current, p_vert_diff);
      _mm_store_ps(u + IX(i, j), u_current);

      __m128 p_right = _mm_load_ps(p + IX(i, j));
      __m128 p_right_temp = _mm_load_ss(p + IX(i, j + 4));
      p_right = _mm_move_ss(p_right, p_right_temp);
      p_right = _mm_shuffle_ps(p_right, p_right, _MM_SHUFFLE(0, 3, 2, 1));

      __m128 p_left = _mm_load_ps(p + IX(i, j));
      p_left = _mm_shuffle_ps(p_left, p_left, _MM_SHUFFLE(2, 1, 0, 0));
      __m128 p_left_temp = _mm_load_ss(p + IX(i, j - 1));
      p_left = _mm_move_ss(p_left, p_left_temp);

      __m128 p_horz_diff = _mm_sub_ps(p_right, p_left);
      p_horz_diff = _mm_mul_ps(p_horz_diff, multiplier);

      __m128 v_current = _mm_load_ps(v + IX(i, j));
      v_current = _mm_sub_ps(v_current, p_horz_diff);
      _mm_store_ps(v + IX(i, j), v_current);
    }
  }
  SET_BND(N, MAT_U_VEL, u);
  SET_BND(N, MAT_V_VEL, v);
}

void project(size_t N, float* u, float* v, float* p, float* div) {
  for (size_t i = ROWBEGIN; i < ROWEND; i++) {
    for (size_t j = COLBEGIN; j < COLEND; j++) {
      div[IX(i, j)] = -0.5f *
                      (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] -
                       v[IX(i, j - 1)]) /
                      N;
      p[IX(i, j)] = 0;
    }
  }
  SET_BND(N, MAT_FLUID, div);
  SET_BND(N, MAT_FLUID, p);

  SOLVE(N, MAT_FLUID, p, div, 1, 4);

  for (size_t i = ROWBEGIN; i < ROWEND; i++) {
    for (size_t j = COLBEGIN; j < COLEND; j++) {
      u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
      v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
    }
  }
  SET_BND(N, MAT_U_VEL, u);
  SET_BND(N, MAT_V_VEL, v);
}

void dens_step(size_t N, float* x, float* x0, float* u, float* v, float diff,
               float dt) {
  ADD_SOURCE(N, x, x0, dt);
  SWAP(x0, x);
  diffuse(N, MAT_FLUID, x, x0, diff, dt);
  SWAP(x0, x);
  advect(N, MAT_FLUID, x, x0, u, v, dt);
}

void vel_step(size_t N, float* u, float* v, float* u0, float* v0, float visc,
              float dt) {
  ADD_SOURCE(N, u, u0, dt);
  ADD_SOURCE(N, v, v0, dt);
  SWAP(u0, u);
  diffuse(N, MAT_U_VEL, u, u0, visc, dt);
  SWAP(v0, v);
  diffuse(N, MAT_V_VEL, v, v0, visc, dt);
  PROJECT(N, u, v, u0, v0);
  SWAP(u0, u);
  SWAP(v0, v);
  advect(N, MAT_U_VEL, u, u0, u0, v0, dt);
  advect(N, MAT_V_VEL, v, v0, u0, v0, dt);
  PROJECT(N, u, v, u0, v0);
}
