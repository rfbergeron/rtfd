#include "solver.h"

#include <immintrin.h>
#include <stdalign.h>
#include <stdlib.h>
#include <string.h>

#define SOLVER sse2_solve
#define SOLVE(N, type, x, x0, a, c) SOLVER(N, type, x, x0, a, c)
#define PROJECTOR sse2_project
#define PROJECT(N, u, v, p, div) PROJECTOR(N, u, v, p, div)
#define ADDER add_source
#define ADD_SOURCE(N, x, s, dt) ADDER(N, x, s, dt)
#define BOUNDER sse2_set_bnd
#define SET_BND(N, type, x) BOUNDER(N, type, x)
#define ADVECTOR sse4_2_advect
#define ADVECT(N, type, d, d0, u, v, dt) ADVECTOR(N, type, d, d0, u, v, dt)

typedef enum matrix_type { MAT_FLUID, MAT_U_VEL, MAT_V_VEL } MatrixType;

#define SWAP(x0, x)  \
  {                  \
    float* tmp = x0; \
    x0 = x;          \
    x = tmp;         \
  }

// NOTE: vectorizing this one by hand is not necessary; it is simple enough that
// the compiler can figure it out on its own
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

// GCC does not hoist the conditions of the ternary operator out of the loops,
// nor does it turn the subps instructions into the magic xorps instruction used
// with scalar operands. It also does not hoist the condition out if the ternary
// operators are replaced with if statements.
void sse2_set_bnd_alt(size_t N, MatrixType type, float* x) {
  for (size_t j = COLBEGIN; j < COLEND; j += 4) {
    _mm_store_ps(
        x + IX(ROWBEGIN - 1, j),
        type == MAT_U_VEL
            ? _mm_sub_ps(_mm_setzero_ps(), _mm_load_ps(x + IX(ROWBEGIN, j)))
            : _mm_load_ps(x + IX(ROWBEGIN, j)));
    _mm_store_ps(
        x + IX(ROWEND, j),
        type == MAT_U_VEL
            ? _mm_sub_ps(_mm_setzero_ps(), _mm_load_ps(x + IX(ROWEND - 1, j)))
            : _mm_load_ps(x + IX(ROWEND - 1, j)));
  }
  for (size_t i = ROWBEGIN; i < ROWEND; ++i) {
    x[IX(i, COLBEGIN - 1)] =
        type == MAT_V_VEL ? -x[IX(i, COLBEGIN)] : x[IX(i, COLBEGIN)];
    x[IX(i, COLEND)] =
        type == MAT_V_VEL ? -x[IX(i, COLEND - 1)] : x[IX(i, COLEND - 1)];
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

void sse2_set_bnd(size_t N, MatrixType type, float* x) {
  switch (type) {
    case MAT_FLUID:
      for (size_t j = COLBEGIN; j < COLEND; j += 4) {
        _mm_store_ps(x + IX(ROWBEGIN - 1, j), _mm_load_ps(x + IX(ROWBEGIN, j)));
        _mm_store_ps(x + IX(ROWEND, j), _mm_load_ps(x + IX(ROWEND - 1, j)));
      }
      for (size_t i = ROWBEGIN; i < ROWEND; ++i) {
        x[IX(i, COLBEGIN - 1)] = x[IX(i, COLBEGIN)];
        x[IX(i, COLEND)] = x[IX(i, COLEND - 1)];
      }
      break;
    case MAT_U_VEL:
      for (size_t j = COLBEGIN; j < COLEND; j += 4) {
        _mm_store_ps(
            x + IX(ROWBEGIN - 1, j),
            _mm_mul_ps(_mm_set1_ps(-1.0f), _mm_load_ps(x + IX(ROWBEGIN, j))));
        _mm_store_ps(
            x + IX(ROWEND, j),
            _mm_mul_ps(_mm_set1_ps(-1.0f), _mm_load_ps(x + IX(ROWEND - 1, j))));
      }
      for (size_t i = ROWBEGIN; i < ROWEND; ++i) {
        x[IX(i, COLBEGIN - 1)] = x[IX(i, COLBEGIN)];
        x[IX(i, COLEND)] = x[IX(i, COLEND - 1)];
      }
      break;
    case MAT_V_VEL:
      for (size_t j = COLBEGIN; j < COLEND; j += 4) {
        _mm_store_ps(x + IX(ROWBEGIN - 1, j), _mm_load_ps(x + IX(ROWBEGIN, j)));
        _mm_store_ps(x + IX(ROWEND, j), _mm_load_ps(x + IX(ROWEND - 1, j)));
      }
      for (size_t i = ROWBEGIN; i < ROWEND; ++i) {
        x[IX(i, COLBEGIN - 1)] = -x[IX(i, COLBEGIN)];
        x[IX(i, COLEND)] = -x[IX(i, COLEND - 1)];
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

// unfortunately, GCC is not smart enough to auto-vectorize this
void set_bnd_alt(size_t N, MatrixType type, float* x) {
  switch (type) {
    case MAT_FLUID:
      for (size_t j = COLBEGIN; j < COLEND; ++j) {
        x[IX(ROWBEGIN - 1, j)] = x[IX(ROWBEGIN, j)];
        x[IX(ROWEND, j)] = x[IX(ROWEND - 1, j)];
      }
      for (size_t i = ROWBEGIN; i < ROWEND; ++i) {
        x[IX(i, COLBEGIN - 1)] = x[IX(i, COLBEGIN)];
        x[IX(i, COLEND)] = x[IX(i, COLEND - 1)];
      }
      break;
    case MAT_U_VEL:
      for (size_t j = COLBEGIN; j < COLEND; ++j) {
        x[IX(ROWBEGIN - 1, j)] = -x[IX(ROWBEGIN, j)];
        x[IX(ROWEND, j)] = -x[IX(ROWEND - 1, j)];
      }
      for (size_t i = ROWBEGIN; i < ROWEND; ++i) {
        x[IX(i, COLBEGIN - 1)] = x[IX(i, COLBEGIN)];
        x[IX(i, COLEND)] = x[IX(i, COLEND - 1)];
      }
      break;
    case MAT_V_VEL:
      for (size_t j = COLBEGIN; j < COLEND; ++j) {
        x[IX(ROWBEGIN - 1, j)] = x[IX(ROWBEGIN, j)];
        x[IX(ROWEND, j)] = x[IX(ROWEND - 1, j)];
      }
      for (size_t i = ROWBEGIN; i < ROWEND; ++i) {
        x[IX(i, COLBEGIN - 1)] = -x[IX(i, COLBEGIN)];
        x[IX(i, COLEND)] = -x[IX(i, COLEND - 1)];
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
  size_t size = ACTUALSIZE;
  if (size % 16 != 0) size += 16 - size % 16;
  float* x1 = aligned_alloc(64, size * sizeof(float));
  __m128 c_inv_vec = _mm_set1_ps(1.0f / c), a_vec = _mm_set1_ps(a);

  for (size_t k = 0; k < 20; k++) {
    for (size_t i = ROWBEGIN; i < ROWEND; ++i) {
      for (size_t j = COLBEGIN; j < COLEND; j += 4) {
        __m128 above = _mm_load_ps(x + IX(i + 1, j));
        __m128 below = _mm_load_ps(x + IX(i - 1, j));

        /* row-by-row addition */
        __m128 dest = _mm_add_ps(above, below);

        /* column-by-column addition */
        __m128 left = _mm_load_ps(x + IX(i, j));
        left = _mm_shuffle_ps(left, left, _MM_SHUFFLE(2, 1, 0, 0));
        left = _mm_move_ss(left, _mm_load_ss(x + IX(i, j - 1)));
        dest = _mm_add_ps(dest, left);

        /* shuffle after since we don't need to save slot zero but we do need
         * the value we just loaded to be in the high slot
         */
        __m128 right = _mm_load_ps(x + IX(i, j));
        right = _mm_move_ss(right, _mm_load_ss(x + IX(i, j + 4)));
        right = _mm_shuffle_ps(right, right, _MM_SHUFFLE(0, 3, 2, 1));
        dest = _mm_add_ps(dest, right);

        /* multiply by a; add x0; divide by c */
        dest = _mm_mul_ps(dest, a_vec);
        dest = _mm_add_ps(dest, _mm_load_ps(x0 + IX(i, j)));
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
  size_t size = ACTUALSIZE;
  if (size % 16 != 0) size += 16 - size % 16;
  float* x1 = aligned_alloc(64, size * sizeof(float));

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

// SSE2 is not powerful enough to meaningfully vectorize this function; the
// oldest extension set that makes sense to use to implement this is SSE4.2
// (ish). SSE4.1 might be enough but really what I'm targeting is GCC's
// x86-64-v3 feature level, which lumps together a bunch of extensions, the
// most recent of which is SSE4.2.
void sse4_2_advect(size_t N, MatrixType type, float* d, float* d0, float* u,
                   float* v, float dt) {
  alignas(16) static const float j_init[4] = {COLBEGIN, COLBEGIN + 1.0f,
                                              COLBEGIN + 2.0f, COLBEGIN + 3.0f},
                                 x_min[4] = {ROWBEGIN - 0.5f, ROWBEGIN - 0.5f,
                                             ROWBEGIN - 0.5f, ROWBEGIN - 0.5f},
                                 y_min[4] = {COLBEGIN - 0.5f, COLBEGIN - 0.5f,
                                             COLBEGIN - 0.5f, COLBEGIN - 0.5f};
  __m128 dt0_vec = _mm_set1_ps(dt * N);
  // NOTE: if the way ROW/COL/END/BEGIN is calculated changes, this will break
  __m128 x_min_vec = _mm_load_ps(x_min);
  __m128 x_max_vec = _mm_add_ps(x_min_vec, _mm_set1_ps(N));
  __m128 y_min_vec = _mm_load_ps(y_min);
  __m128 y_max_vec = _mm_add_ps(y_min_vec, _mm_set1_ps(N));
  for (size_t i = ROWBEGIN; i < ROWEND; i++) {
    __m128 i_vec = _mm_set1_ps(i);
    __m128 j_vec = _mm_load_ps(j_init);
    for (size_t j = COLBEGIN; j < COLEND;
         j += 4, j_vec = _mm_add_ps(j_vec, _mm_set1_ps(4.0f))) {
      __m128 u_vec = _mm_load_ps(u + IX(i, j));
      __m128 x_vec = _mm_sub_ps(i_vec, _mm_mul_ps(dt0_vec, u_vec));
      x_vec = _mm_min_ps(_mm_max_ps(x_vec, x_min_vec), x_max_vec);
      __m128 s1_vec = _mm_sub_ps(x_vec, _mm_floor_ps(x_vec));
      __m128 s0_vec = _mm_sub_ps(_mm_set1_ps(1.0f), s1_vec);
      float i0s[4];
      _mm_store_ps(i0s, x_vec);

      __m128 v_vec = _mm_load_ps(v + IX(i, j));
      __m128 y_vec = _mm_sub_ps(j_vec, _mm_mul_ps(dt0_vec, v_vec));
      y_vec = _mm_min_ps(_mm_max_ps(y_vec, y_min_vec), y_max_vec);
      __m128 t1_vec = _mm_sub_ps(y_vec, _mm_floor_ps(y_vec));
      __m128 t0_vec = _mm_sub_ps(_mm_set1_ps(1.0f), t1_vec);
      float j0s[4];
      _mm_store_ps(j0s, y_vec);

      // this sucks
      __m128 d0_i0_j0 = _mm_set_ps(d0[IX((size_t)i0s[3], (size_t)j0s[3])],
                                   d0[IX((size_t)i0s[2], (size_t)j0s[2])],
                                   d0[IX((size_t)i0s[1], (size_t)j0s[1])],
                                   d0[IX((size_t)i0s[0], (size_t)j0s[0])]);
      __m128 d0_i1_j0 = _mm_set_ps(d0[IX((size_t)i0s[3] + 1, (size_t)j0s[3])],
                                   d0[IX((size_t)i0s[2] + 1, (size_t)j0s[2])],
                                   d0[IX((size_t)i0s[1] + 1, (size_t)j0s[1])],
                                   d0[IX((size_t)i0s[0] + 1, (size_t)j0s[0])]);
      __m128 d0_i0_j1 = _mm_set_ps(d0[IX((size_t)i0s[3], (size_t)j0s[3] + 1)],
                                   d0[IX((size_t)i0s[2], (size_t)j0s[2] + 1)],
                                   d0[IX((size_t)i0s[1], (size_t)j0s[1] + 1)],
                                   d0[IX((size_t)i0s[0], (size_t)j0s[0] + 1)]);
      __m128 d0_i1_j1 =
          _mm_set_ps(d0[IX((size_t)i0s[3] + 1, (size_t)j0s[3] + 1)],
                     d0[IX((size_t)i0s[2] + 1, (size_t)j0s[2] + 1)],
                     d0[IX((size_t)i0s[1] + 1, (size_t)j0s[1] + 1)],
                     d0[IX((size_t)i0s[0] + 1, (size_t)j0s[0] + 1)]);

      _mm_store_ps(
          d + IX(i, j),
          _mm_add_ps(
              _mm_mul_ps(s0_vec, _mm_add_ps(_mm_mul_ps(t0_vec, d0_i0_j0),
                                            _mm_mul_ps(t1_vec, d0_i0_j1))),
              _mm_mul_ps(s1_vec, _mm_add_ps(_mm_mul_ps(t0_vec, d0_i1_j0),
                                            _mm_mul_ps(t1_vec, d0_i1_j1)))));
    }
  }
  SET_BND(N, type, d);
}

// the compiler is smart enough to turn the if statements into calls to minss
// and maxss; there is little to no benefit to using them explicitly. you could
// use roundss and subss to get the fractional component of x and y, but all
// you'd be saving is one int-to-float conversion, which has roughly the same
// latency and throughput. we do however do that above because we should really
// be converting to 64-bit integers, which are too big to fit 4 of them into a
// single xmm register. there also isn't a way to convert packed floats to
// packed 64-bit integers until AVX512.
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
      v_right = _mm_move_ss(v_right, _mm_load_ss(v + IX(i, j + 4)));
      v_right = _mm_shuffle_ps(v_right, v_right, _MM_SHUFFLE(0, 3, 2, 1));

      __m128 v_left = _mm_load_ps(v + IX(i, j));
      v_left = _mm_shuffle_ps(v_left, v_left, _MM_SHUFFLE(2, 1, 0, 0));
      v_left = _mm_move_ss(v_left, _mm_load_ss(v + IX(i, j - 1)));

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
  ADVECT(N, MAT_FLUID, x, x0, u, v, dt);
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
  ADVECT(N, MAT_U_VEL, u, u0, u0, v0, dt);
  ADVECT(N, MAT_V_VEL, v, v0, u0, v0, dt);
  PROJECT(N, u, v, u0, v0);
}
