#include "solver.h"

#include <immintrin.h>
#include <stdalign.h>
#include <stdlib.h>

#define IX(i, j) ((sim_size + 2 * col_border) * (i) + (j))
#define ACTUAL_SIZE ((sim_size + 2 * col_border) * (sim_size + 2 * row_border))
#define COL_BEGIN col_border
#define COL_END (sim_size + col_border)
#define ROW_BEGIN row_border
#define ROW_END (sim_size + row_border)

#define SOLVER sse2_solve
#define SOLVE(sim_size, row_border, col_border, type, x, x0, x1, a, c) \
  SOLVER(sim_size, row_border, col_border, type, x, x0, x1, a, c)
#define PROJECTOR sse2_project
#define PROJECT(sim_size, row_border, col_border, u, v, p, div, scratch) \
  PROJECTOR(sim_size, row_border, col_border, u, v, p, div, scratch)
#define ADDER add_source
#define ADD_SOURCE(sim_size, row_border, col_border, x, s, dt, multiplier) \
  ADDER(sim_size, row_border, col_border, x, s, dt, multiplier)
#define BOUNDER sse2_set_bnd
#define SET_BND(sim_size, row_border, col_border, type, x) \
  BOUNDER(sim_size, row_border, col_border, type, x)
#define ADVECTOR sse4_2_advect
#define ADVECT(sim_size, row_border, col_border, type, d, d0, u, v, dt) \
  ADVECTOR(sim_size, row_border, col_border, type, d, d0, u, v, dt)
#define SWAP(x0, x)  \
  {                  \
    float* tmp = x0; \
    x0 = x;          \
    x = tmp;         \
  }

// NOTE: vectorizing this one by hand is not necessary; it is simple enough that
// the compiler can figure it out on its own
static void sse2_add_source(size_t sim_size, size_t row_border,
                            size_t col_border, float* x, float* s, float dt,
                            float multiplier) {
  __m128 dt_vec = _mm_set1_ps(dt * multiplier);
  for (size_t i = 0, size = ACTUAL_SIZE; i < size; i += 4) {
    __m128 x_current = _mm_load_ps(x + i), s_current = _mm_load_ps(s + i);
    s_current = _mm_mul_ps(s_current, dt_vec);
    x_current = _mm_add_ps(x_current, s_current);
    _mm_store_ps(x + i, x_current);
  }
}

static void add_source(size_t sim_size, size_t row_border, size_t col_border,
                       float* x, float* s, float dt, float multiplier) {
  for (size_t i = 0, size = ACTUAL_SIZE; i < size; i++)
    x[i] += dt * multiplier * s[i];
}

// GCC does not hoist the conditions of the ternary operator out of the loops,
// nor does it turn the subps instructions into the magic xorps instruction used
// with scalar operands. It also does not hoist the condition out if the ternary
// operators are replaced with if statements.
static void sse2_set_bnd_alt(size_t sim_size, size_t row_border,
                             size_t col_border, MatrixType type, float* x) {
  for (size_t j = COL_BEGIN; j < COL_END; j += 4) {
    _mm_store_ps(
        x + IX(ROW_BEGIN - 1, j),
        (type == SLV_MAT_U || type == SLV_MAT_U0)
            ? _mm_sub_ps(_mm_setzero_ps(), _mm_load_ps(x + IX(ROW_BEGIN, j)))
            : _mm_load_ps(x + IX(ROW_BEGIN, j)));
    _mm_store_ps(
        x + IX(ROW_END, j),
        (type == SLV_MAT_U || type == SLV_MAT_U0)
            ? _mm_sub_ps(_mm_setzero_ps(), _mm_load_ps(x + IX(ROW_END - 1, j)))
            : _mm_load_ps(x + IX(ROW_END - 1, j)));
  }
  for (size_t i = ROW_BEGIN; i < ROW_END; ++i) {
    x[IX(i, COL_BEGIN - 1)] = (type == SLV_MAT_V || type == SLV_MAT_V0)
                                  ? -x[IX(i, COL_BEGIN)]
                                  : x[IX(i, COL_BEGIN)];
    x[IX(i, COL_END)] = (type == SLV_MAT_V || type == SLV_MAT_V0)
                            ? -x[IX(i, COL_END - 1)]
                            : x[IX(i, COL_END - 1)];
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

static void sse2_set_bnd(size_t sim_size, size_t row_border, size_t col_border,
                         MatrixType type, float* x) {
  switch (type) {
    case SLV_MAT_D:
    case SLV_MAT_D0:
      for (size_t j = COL_BEGIN; j < COL_END; j += 4) {
        _mm_store_ps(x + IX(ROW_BEGIN - 1, j),
                     _mm_load_ps(x + IX(ROW_BEGIN, j)));
        _mm_store_ps(x + IX(ROW_END, j), _mm_load_ps(x + IX(ROW_END - 1, j)));
      }
      for (size_t i = ROW_BEGIN; i < ROW_END; ++i) {
        x[IX(i, COL_BEGIN - 1)] = x[IX(i, COL_BEGIN)];
        x[IX(i, COL_END)] = x[IX(i, COL_END - 1)];
      }
      break;
    case SLV_MAT_U:
    case SLV_MAT_U0:
      for (size_t j = COL_BEGIN; j < COL_END; j += 4) {
        _mm_store_ps(
            x + IX(ROW_BEGIN - 1, j),
            _mm_mul_ps(_mm_set1_ps(-1.0f), _mm_load_ps(x + IX(ROW_BEGIN, j))));
        _mm_store_ps(x + IX(ROW_END, j),
                     _mm_mul_ps(_mm_set1_ps(-1.0f),
                                _mm_load_ps(x + IX(ROW_END - 1, j))));
      }
      for (size_t i = ROW_BEGIN; i < ROW_END; ++i) {
        x[IX(i, COL_BEGIN - 1)] = x[IX(i, COL_BEGIN)];
        x[IX(i, COL_END)] = x[IX(i, COL_END - 1)];
      }
      break;
    case SLV_MAT_V:
    case SLV_MAT_V0:
      for (size_t j = COL_BEGIN; j < COL_END; j += 4) {
        _mm_store_ps(x + IX(ROW_BEGIN - 1, j),
                     _mm_load_ps(x + IX(ROW_BEGIN, j)));
        _mm_store_ps(x + IX(ROW_END, j), _mm_load_ps(x + IX(ROW_END - 1, j)));
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

// unfortunately, GCC is not smart enough to auto-vectorize this
static void set_bnd_alt(size_t sim_size, size_t row_border, size_t col_border,
                        MatrixType type, float* x) {
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

static void set_bnd(size_t sim_size, size_t row_border, size_t col_border,
                    MatrixType type, float* x) {
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

void sse2_solve(size_t sim_size, size_t row_border, size_t col_border,
                MatrixType type, float* x, float* x0, float* x1, float a,
                float c) {
  __m128 c_inv_vec = _mm_set1_ps(1.0f / c), a_vec = _mm_set1_ps(a);

  for (size_t k = 0; k < 20; k++) {
    for (size_t i = ROW_BEGIN; i < ROW_END; ++i) {
      for (size_t j = COL_BEGIN; j < COL_END; j += 4) {
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
    SET_BND(sim_size, row_border, col_border, type, x1);
    SWAP(x, x1);
  }
}

static void jac_solve(size_t sim_size, size_t row_border, size_t col_border,
                      MatrixType type, float* x, float* x0, float* x1, float a,
                      float c) {
  for (size_t k = 0; k < 20; k++) {
    for (size_t i = ROW_BEGIN; i < ROW_END; i++) {
      for (size_t j = COL_BEGIN; j < COL_END; j++) {
        x1[IX(i, j)] =
            (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                 x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
            c;
      }
    }
    SET_BND(sim_size, row_border, col_border, type, x1);
    SWAP(x, x1);
  }
}

static void lin_solve(size_t sim_size, size_t row_border, size_t col_border,
                      MatrixType type, float* x, float* x0, float* x1, float a,
                      float c) {
  (void)x1;  // unused parameter
  for (size_t k = 0; k < 20; k++) {
    for (size_t i = ROW_BEGIN; i < ROW_END; i++) {
      for (size_t j = COL_BEGIN; j < COL_END; j++) {
        x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                           x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
                      c;
      }
    }
    SET_BND(sim_size, row_border, col_border, type, x);
  }
}

static void diffuse(size_t sim_size, size_t row_border, size_t col_border,
                    MatrixType type, float* x, float* x0, float* x1, float diff,
                    float dt) {
  float a = dt * diff * sim_size * sim_size;
  SOLVE(sim_size, row_border, col_border, type, x, x0, x1, a, 1 + 4 * a);
}

// SSE2 is not powerful enough to meaningfully vectorize this function; the
// oldest extension set that makes sense to use to implement this is SSE4.2
// (ish). SSE4.1 might be enough but really what I'm targeting is GCC's
// x86-64-v3 feature level, which lumps together a bunch of extensions, the
// most recent of which is SSE4.2.
static void sse4_2_advect(size_t sim_size, size_t row_border, size_t col_border,
                          MatrixType type, float* d, float* d0, float* u,
                          float* v, float dt) {
  __m128 dt0_vec = _mm_set1_ps(dt * sim_size);
  __m128 x_min_vec = _mm_set1_ps(ROW_BEGIN - 0.5f);
  __m128 x_max_vec = _mm_set1_ps(ROW_END - 0.5f);
  __m128 y_min_vec = _mm_set1_ps(COL_BEGIN - 0.5f);
  __m128 y_max_vec = _mm_set1_ps(COL_END - 0.5f);
  for (size_t i = ROW_BEGIN; i < ROW_END; i++) {
    __m128 i_vec = _mm_set1_ps(i);
    __m128 j_vec = _mm_set_ps(COL_BEGIN + 3.0f, COL_BEGIN + 2.0f,
                              COL_BEGIN + 1.0f, COL_BEGIN);
    for (size_t j = COL_BEGIN; j < COL_END;
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
  SET_BND(sim_size, row_border, col_border, type, d);
}

// the compiler is smart enough to turn the if statements into calls to minss
// and maxss; there is little to no benefit to using them explicitly. you could
// use roundss and subss to get the fractional component of x and y, but all
// you'd be saving is one int-to-float conversion, which has roughly the same
// latency and throughput. we do however do that above because we should really
// be converting to 64-bit integers, which are too big to fit 4 of them into a
// single xmm register. there also isn't a way to convert packed floats to
// packed 64-bit integers until AVX512.
static void advect(size_t sim_size, size_t row_border, size_t col_border,
                   MatrixType type, float* d, float* d0, float* u, float* v,
                   float dt) {
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
  SET_BND(sim_size, row_border, col_border, type, d);
}

static void sse2_project(size_t sim_size, size_t row_border, size_t col_border,
                         float* u, float* v, float* p, float* div,
                         float* scratch) {
  __m128 multiplier = _mm_set1_ps(-0.5f / sim_size);
  for (size_t i = 0, size = ACTUAL_SIZE; i < size; i++) p[i] = 0.0f;
  for (size_t i = ROW_BEGIN; i < ROW_END; i++) {
    for (size_t j = COL_BEGIN; j < COL_END; j += 4) {
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
  SET_BND(sim_size, row_border, col_border, SLV_MAT_D, div);

  SOLVE(sim_size, row_border, col_border, SLV_MAT_D, p, div, scratch, 1, 4);

  multiplier = _mm_set1_ps(0.5f * sim_size);
  for (size_t i = ROW_BEGIN; i < ROW_END; i++) {
    for (size_t j = COL_BEGIN; j < COL_END; j += 4) {
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
  SET_BND(sim_size, row_border, col_border, SLV_MAT_U, u);
  SET_BND(sim_size, row_border, col_border, SLV_MAT_V, v);
}

static void project(size_t sim_size, size_t row_border, size_t col_border,
                    float* u, float* v, float* p, float* div, float* scratch) {
  for (size_t i = ROW_BEGIN; i < ROW_END; i++) {
    for (size_t j = COL_BEGIN; j < COL_END; j++) {
      div[IX(i, j)] = -0.5f *
                      (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] -
                       v[IX(i, j - 1)]) /
                      sim_size;
      p[IX(i, j)] = 0;
    }
  }
  SET_BND(sim_size, row_border, col_border, SLV_MAT_D, div);
  SET_BND(sim_size, row_border, col_border, SLV_MAT_D, p);

  SOLVE(sim_size, row_border, col_border, SLV_MAT_D, p, div, scratch, 1, 4);

  for (size_t i = ROW_BEGIN; i < ROW_END; i++) {
    for (size_t j = COL_BEGIN; j < COL_END; j++) {
      u[IX(i, j)] -= 0.5f * sim_size * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
      v[IX(i, j)] -= 0.5f * sim_size * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
    }
  }
  SET_BND(sim_size, row_border, col_border, SLV_MAT_U, u);
  SET_BND(sim_size, row_border, col_border, SLV_MAT_V, v);
}

void solver_dens_step(Solver* solver) {
  ADD_SOURCE(solver->sim_size, solver->row_border, solver->col_border,
             solver->d, solver->d_prev, solver->dt, solver->source);
  SWAP(solver->d_prev, solver->d);
  diffuse(solver->sim_size, solver->row_border, solver->col_border, SLV_MAT_D,
          solver->d, solver->d_prev, solver->u_prev, solver->diff, solver->dt);
  SWAP(solver->d_prev, solver->d);
  ADVECT(solver->sim_size, solver->row_border, solver->col_border, SLV_MAT_D,
         solver->d, solver->d_prev, solver->u, solver->v, solver->dt);
}

void solver_vel_step(Solver* solver) {
  ADD_SOURCE(solver->sim_size, solver->row_border, solver->col_border,
             solver->u, solver->u_prev, solver->dt, solver->force);
  ADD_SOURCE(solver->sim_size, solver->row_border, solver->col_border,
             solver->v, solver->v_prev, solver->dt, solver->force);
  SWAP(solver->u_prev, solver->u);
  diffuse(solver->sim_size, solver->row_border, solver->col_border, SLV_MAT_U,
          solver->u, solver->u_prev, solver->d_prev, solver->visc, solver->dt);
  SWAP(solver->v_prev, solver->v);
  diffuse(solver->sim_size, solver->row_border, solver->col_border, SLV_MAT_V,
          solver->v, solver->v_prev, solver->d_prev, solver->visc, solver->dt);
  PROJECT(solver->sim_size, solver->row_border, solver->col_border, solver->u,
          solver->v, solver->u_prev, solver->v_prev, solver->d_prev);
  SWAP(solver->u_prev, solver->u);
  SWAP(solver->v_prev, solver->v);
  ADVECT(solver->sim_size, solver->row_border, solver->col_border, SLV_MAT_U,
         solver->u, solver->u_prev, solver->u_prev, solver->v_prev, solver->dt);
  ADVECT(solver->sim_size, solver->row_border, solver->col_border, SLV_MAT_V,
         solver->v, solver->v_prev, solver->u_prev, solver->v_prev, solver->dt);
  PROJECT(solver->sim_size, solver->row_border, solver->col_border, solver->u,
          solver->v, solver->u_prev, solver->v_prev, solver->d_prev);
}

float* solver_at(Solver* solver, size_t x, size_t y, MatrixType type) {
  const size_t sim_size = solver->sim_size, row_border = solver->row_border,
               col_border = solver->col_border,
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

void solver_clear(Solver* solver, MatrixType type) {
  const size_t sim_size = solver->sim_size, row_border = solver->row_border,
               col_border = solver->col_border;
  for (size_t i = 0, size = ACTUAL_SIZE; i < size; i++) {
    if (type & SLV_MAT_D) solver->d[i] = 0.0f;
    if (type & SLV_MAT_U) solver->u[i] = 0.0f;
    if (type & SLV_MAT_V) solver->v[i] = 0.0f;
    if (type & SLV_MAT_D0) solver->d_prev[i] = 0.0f;
    if (type & SLV_MAT_U0) solver->u_prev[i] = 0.0f;
    if (type & SLV_MAT_V0) solver->v_prev[i] = 0.0f;
  }
}

Solver* solver_init(size_t sim_size, float dt, float diff, float visc,
                    float force, float source) {
  const size_t row_border = 1, col_border = 4;
  size_t size = ACTUAL_SIZE;

  Solver* solver = malloc(sizeof(Solver));
  if (!solver) return NULL;

  solver->sim_size = sim_size;
  solver->row_border = row_border;
  solver->col_border = col_border;
  solver->dt = dt;
  solver->diff = diff;
  solver->visc = visc;
  solver->force = force;
  solver->source = source;
  solver->u = aligned_alloc(64, size * sizeof(float));
  solver->v = aligned_alloc(64, size * sizeof(float));
  solver->u_prev = aligned_alloc(64, size * sizeof(float));
  solver->v_prev = aligned_alloc(64, size * sizeof(float));
  solver->d = aligned_alloc(64, size * sizeof(float));
  solver->d_prev = aligned_alloc(64, size * sizeof(float));

  if (!solver->u || !solver->v || !solver->u_prev || !solver->v_prev ||
      !solver->d || !solver->d_prev) {
    solver_destroy(solver);
    return NULL;
  }

  return solver;
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
