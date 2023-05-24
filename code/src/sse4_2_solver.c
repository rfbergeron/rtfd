#include "sse4_2_solver.h"

#include <immintrin.h>
#include <stdlib.h>

#define ROW_BORDER 1
#define COL_BORDER 4
#define IX(i, j) ((sim_size + 2 * COL_BORDER) * (i) + (j))
#define ACTUAL_SIZE ((sim_size + 2 * COL_BORDER) * (sim_size + 2 * ROW_BORDER))
#define COL_BEGIN COL_BORDER
#define COL_END (sim_size + COL_BORDER)
#define ROW_BEGIN ROW_BORDER
#define ROW_END (sim_size + ROW_BORDER)
#define SWAP(x0, x)  \
  {                  \
    float* tmp = x0; \
    x0 = x;          \
    x = tmp;         \
  }

// NOTE: vectorizing this one by hand is not necessary; it is simple enough that
// the compiler can figure it out on its own
/*
static void sse2_add_source(size_t sim_size, float* x, float* s, float dt,
                            float multiplier) {
  __m128 dt_vec = _mm_set1_ps(dt * multiplier);
  for (size_t i = 0, size = ACTUAL_SIZE; i < size; i += 4) {
    __m128 x_current = _mm_load_ps(x + i), s_current = _mm_load_ps(s + i);
    s_current = _mm_mul_ps(s_current, dt_vec);
    x_current = _mm_add_ps(x_current, s_current);
    _mm_store_ps(x + i, x_current);
  }
}
*/

static void add_source(size_t sim_size, float* x, float* s, float dt,
                       float multiplier) {
  for (size_t i = 0, size = ACTUAL_SIZE; i < size; i++)
    x[i] += dt * multiplier * s[i];
}

// GCC does not hoist the conditions of the ternary operator out of the loops,
// nor does it turn the subps instructions into the magic xorps instruction used
// with scalar operands. It also does not hoist the condition out if the ternary
// operators are replaced with if statements.
/*
static void sse2_set_bnd_alt(size_t sim_size, MatrixType type, float* x) {
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
*/

static void sse2_set_bnd(size_t sim_size, MatrixType type, float* x) {
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

void sse2_solve(size_t sim_size, MatrixType type, float* x, float* x0,
                float* x1, float a, float c) {
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
    sse2_set_bnd(sim_size, type, x1);
    SWAP(x, x1);
  }
}

static void diffuse(size_t sim_size, MatrixType type, float* x, float* x0,
                    float* x1, float diff, float dt) {
  float a = dt * diff * sim_size * sim_size;
  sse2_solve(sim_size, type, x, x0, x1, a, 1 + 4 * a);
}

// SSE2 is not powerful enough to meaningfully vectorize this function; the
// oldest extension set that makes sense to use to implement this is SSE4.2
// (ish). SSE4.1 might be enough but really what I'm targeting is GCC's
// x86-64-v3 feature level, which lumps together a bunch of extensions, the
// most recent of which is SSE4.2.
static void sse4_2_advect(size_t sim_size, MatrixType type, float* d, float* d0,
                          float* u, float* v, float dt) {
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
  sse2_set_bnd(sim_size, type, d);
}

static void sse2_project(size_t sim_size, float* u, float* v, float* p,
                         float* div, float* scratch) {
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
  sse2_set_bnd(sim_size, SLV_MAT_D, div);

  sse2_solve(sim_size, SLV_MAT_D, p, div, scratch, 1, 4);

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
  sse2_set_bnd(sim_size, SLV_MAT_U, u);
  sse2_set_bnd(sim_size, SLV_MAT_V, v);
}

static void dens_step(Solver* solver) {
  add_source(solver->sim_size, solver->d, solver->d_prev, solver->dt,
             solver->source);
  SWAP(solver->d_prev, solver->d);
  diffuse(solver->sim_size, SLV_MAT_D, solver->d, solver->d_prev,
          solver->u_prev, solver->diff, solver->dt);
  SWAP(solver->d_prev, solver->d);
  sse4_2_advect(solver->sim_size, SLV_MAT_D, solver->d, solver->d_prev,
                solver->u, solver->v, solver->dt);
}

static void vel_step(Solver* solver) {
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
  sse2_project(solver->sim_size, solver->u, solver->v, solver->u_prev,
               solver->v_prev, solver->d_prev);
  SWAP(solver->u_prev, solver->u);
  SWAP(solver->v_prev, solver->v);
  sse4_2_advect(solver->sim_size, SLV_MAT_U, solver->u, solver->u_prev,
                solver->u_prev, solver->v_prev, solver->dt);
  sse4_2_advect(solver->sim_size, SLV_MAT_V, solver->v, solver->v_prev,
                solver->u_prev, solver->v_prev, solver->dt);
  sse2_project(solver->sim_size, solver->u, solver->v, solver->u_prev,
               solver->v_prev, solver->d_prev);
}

static float* ix(Solver* solver, size_t x, size_t y, MatrixType type) {
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

static void clear(Solver* solver, MatrixType type) {
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

Solver* sse4_2_solver_init(size_t sim_size, float dt, float diff, float visc,
                           float force, float source) {
  size_t size = ACTUAL_SIZE;

  Solver* solver = malloc(sizeof(Solver));
  if (!solver) return NULL;

  solver->sim_size = sim_size;
  solver->dt = dt;
  solver->diff = diff;
  solver->visc = visc;
  solver->force = force;
  solver->source = source;
  solver->vel_fn = vel_step;
  solver->dens_fn = dens_step;
  solver->clear_fn = clear;
  solver->ix_fn = ix;
  // minimum alignment for xmm is 16; align to 64 because cache
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
