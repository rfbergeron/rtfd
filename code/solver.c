#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxlib.h"
#define BORDER 8
#define IX(i, j) ((i) + (N + BORDER) * (j))
#define ACTUALSIZE ((N + BORDER) * (N + BORDER))
#define BSTART ((BORDER / 2) - 1)
#define BEND (BSTART + N + 1)
#define SWAP(x0, x)  \
  {                  \
    float *tmp = x0; \
    x0 = x;          \
    x = tmp;         \
  }

// TODO(Robert): change all occurrences of 'b' to 'flag' for clarity
// TODO(Robert): change storage scheme for fluid velocity vectors to make
// simd implementation easier and more sensible
// TODO(Robert): change types and declaration locations to be more appropriate
// for 'modern' C

void add_source(size_t N, int flag, float *x, float *s, float dt) {
  size_t size = ACTUALSIZE;
  if (flag) {
   // x,s are uv arrays
    for (size_t i = 0; i < 2 * size; ++i) x[i] += dt * s[i];
  } else {
    // x,s are density arrays
    for (size_t i = 0; i < size; ++i) x[i] += dt * s[i];
  }
}

void set_bnd(size_t N, int flag, float *x) {
  if (flag) {
    // setting vector bounds; x is twice as wide as it is tall; even indices are
    // horizontal(u) components while odd indices are vertical(v) components
    for (size_t i = BSTART + 1; i < BEND; i+=2) {
      // u components
      x[IX(BSTART, i)] = x[IX(BSTART + 1, i)];
      x[IX(BEND, i)] = x[IX(BEND - 1, i)];
      x[IX(i, BSTART)] = -x[IX(i, BSTART + 1)];
      x[IX(i, BEND)] = -x[IX(i, BEND - 1)];
      // horizontal offset for v components
      x[IX(BSTART, i+1)] = -x[IX(BSTART + 1, i)];
      x[IX(BEND, i+1)] = -x[IX(BEND - 1, i)];
      x[IX(i, BSTART+1)] = x[IX(i, BSTART + 1)];
      x[IX(i, BEND+1)] = x[IX(i, BEND - 1)];
    }

    // TODO(Robert): figure out how to correctly do the math to account for the
    // width of the uv matrix
    x[IX(BSTART, BSTART)] =
        0.5f * (x[IX(BSTART + 1, BSTART)] + x[IX(BSTART, BSTART + 1)]);
    x[IX(BSTART, BEND)] =
        0.5f * (x[IX(BSTART + 1, BEND)] + x[IX(BSTART, BEND - 1)]);
    x[IX(2 * BEND, BSTART)] =
        0.5f * (x[IX(BEND - 1, BSTART)] + x[IX(BEND, BSTART + 1)]);
    x[IX(2 * BEND, BEND)] = 0.5f * (x[IX(BEND - 1, BEND)] + x[IX(BEND, BEND - 1)]);
  } else {
    for (size_t i = BSTART + 1; i < BEND; ++i) {
      // iterating over column; v components negated
      x[IX(BSTART, i)] = flag == 1 ? -x[IX(BSTART + 1, i)] : x[IX(BSTART + 1, i)];
      x[IX(BEND, i)] = flag == 1 ? -x[IX(BEND - 1, i)] : x[IX(BEND - 1, i)];
      // iterating over row; u components negated
      x[IX(i, BSTART)] = flag == 2 ? -x[IX(i, BSTART + 1)] : x[IX(i, BSTART + 1)];
      x[IX(i, BEND)] = flag == 2 ? -x[IX(i, BEND - 1)] : x[IX(i, BEND - 1)];
    }

    x[IX(BSTART, BSTART)] =
        0.5f * (x[IX(BSTART + 1, BSTART)] + x[IX(BSTART, BSTART + 1)]);
    x[IX(BSTART, BEND)] =
        0.5f * (x[IX(BSTART + 1, BEND)] + x[IX(BSTART, BEND - 1)]);
    x[IX(BEND, BSTART)] =
        0.5f * (x[IX(BEND - 1, BSTART)] + x[IX(BEND, BSTART + 1)]);
    x[IX(BEND, BEND)] = 0.5f * (x[IX(BEND - 1, BEND)] + x[IX(BEND, BEND - 1)]);
  }
}

void lin_solve(int N, int b, float *x, float *x0, float a, float c) {
  int i, j, k;

  for (k = 0; k < 20; k++) {
    for (i = 1; i <= N; i++) {
      for (j = 1; j <= N; j++) {
        x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                           x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
                      c;
      }
    }
    set_bnd(N, b, x);
  }
}

void jac_solve(int N, int b, float *x, float *x0, float a, float c) {
  int i, j, k;
  size_t size = ACTUALSIZE;
  float *x1 = malloc(size * sizeof(float));

  for (k = 0; k < 20; k++) {
    for (i = BSTART + 1; i < BEND; i++) {
      for (j = BSTART + 1; j < BEND; j++) {
        x1[IX(i, j)] =
            (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                 x[IX(i, j - 1)] + x[IX(i, j + 1)])) /
            c;
      }
    }
    set_bnd(N, b, x1);
    SWAP(x, x1);
  }
  free(x1);
}

void simd_solve(size_t N, int b, float *x, float *x0, float a, float c) {
  size_t i, j, k;
  size_t size = ACTUALSIZE;
  /* assume cache size is 64 or a multiple of 64 and align to it
   */
  float *x1 = aligned_alloc(64, size * sizeof(float));
  float c_inv = 1.0f / c;
  /* no broadcast available until AVX; if we're targeting SSE2 this is
   * something like what _mm_set1_ps translates to
   */
  __m128 c_inv_vec = _mm_load_ss(&c_inv);
  c_inv_vec = _mm_shuffle_ps(c_inv_vec, c_inv_vec, _MM_SHUFFLE(0, 0, 0, 0));
  __m128 a_vec = _mm_load_ss(&a);
  a_vec = _mm_shuffle_ps(a_vec, a_vec, _MM_SHUFFLE(0, 0, 0, 0));

  /* increment i(columns) by four and work with matrix in 4x1 chunks
   * had to switch i and j so that it goes across rows first
   */
  for (k = 0; k < 20; k++) {
    for (j = BSTART + 1; j < BEND; ++j) {
      for (i = BSTART + 1; i < BEND; i += 4) {
        __m128 above = _mm_load_ps(x + IX(i, j - 1));
        __m128 current = _mm_load_ps(x + IX(i, j));
        __m128 below = _mm_load_ps(x + IX(i, j + 1));
        __m128 nought = _mm_load_ps(x0 + IX(i, j));

        /* row-by-row addition */
        __m128 dest = _mm_add_ps(above, below);

        /* column-by-column addition */
        __m128 shiftr =
            _mm_shuffle_ps(current, current, _MM_SHUFFLE(0, 0, 1, 2));
        __m128 tempr = _mm_load_ss(x + IX(i - 1, j));
        shiftr = _mm_move_ss(shiftr, tempr);
        dest = _mm_add_ps(dest, shiftr);

        /* shuffle after since we don't need to save slot zero but we do need
         * the value we just loaded to be in the high slot
         */
        __m128 templ = _mm_load_ss(x + IX(i + 4, j));
        __m128 shiftl = _mm_move_ss(current, templ);
        shiftl = _mm_shuffle_ps(shiftl, shiftl, _MM_SHUFFLE(1, 2, 3, 0));
        dest = _mm_add_ps(dest, shiftl);

        /* multiply by a; add x0; divide by c */
        dest = _mm_mul_ps(dest, a_vec);
        dest = _mm_add_ps(dest, nought);
        dest = _mm_mul_ps(dest, c_inv_vec);

        /* send it back */
        _mm_store_ps(x + IX(i, j), dest);
      }
    }
    set_bnd(N, b, x1);
    SWAP(x, x1);
  }
  free(x1);
}

void diffuse(int N, int b, float *x, float *x0, float diff, float dt) {
  float a = dt * diff * N * N;
  simd_solve(N, b, x, x0, a, 1 + 4 * a);
}

void advect(int N, int b, float *d, float *d0, float *u, float *v, float dt) {
  int i, j, i0, j0, i1, j1;
  float x, y, s0, t0, s1, t1, dt0;

  dt0 = dt * N;
  for (i = BSTART + 1; i < BEND; ++i) {
    for (j = BSTART + 1; j < BEND; ++j) {
      x = i - dt0 * u[IX(i, j)];
      y = j - dt0 * v[IX(i, j)];
      if (x < BSTART + 0.5f) x = BSTART + 0.5f;
      if (x > BEND - 0.5f) x = BEND - 0.5f;
      i0 = (int)x;
      i1 = i0 + 1;
      if (y < BSTART + 0.5f) y = BSTART + 0.5f;
      if (y > BEND - 0.5f) y = BEND - 0.5f;
      j0 = (int)y;
      j1 = j0 + 1;
      s1 = x - i0;
      s0 = 1 - s1;
      t1 = y - j0;
      t0 = 1 - t1;
      d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                    s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
    }
  }
  set_bnd(N, b, d);
}

void project(int N, float *u, float *v, float *p, float *div) {
  int i, j;

  for (i = BSTART + 1; i < BEND; ++i) {
    for (j = BSTART + 1; j < BEND; ++j) {
      div[IX(i, j)] = -0.5f *
                      (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] -
                       v[IX(i, j - 1)]) /
                      N;
      p[IX(i, j)] = 0;
    }
  }
  set_bnd(N, 0, div);
  set_bnd(N, 0, p);

  simd_solve(N, 0, p, div, 1, 4);

  for (i = BSTART + 1; i < BEND; ++i) {
    for (j = BSTART + 1; j < BEND; ++j) {
      u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
      v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
    }
  }
  set_bnd(N, 1, u);
  set_bnd(N, 2, v);
}

void dens_step(int N, float *x, float *x0, float *u, float *v, float diff,
               float dt) {
  DEBUGS('s', "beginning density step");
  add_source(N, x, x0, dt);
  SWAP(x0, x);
  diffuse(N, 0, x, x0, diff, dt);
  SWAP(x0, x);
  advect(N, 0, x, x0, u, v, dt);
}

void vel_step(int N, float *u, float *v, float *u0, float *v0, float visc,
              float dt) {
  DEBUGS('s', "beginning velocity step");
  add_source(N, u, u0, dt);
  add_source(N, v, v0, dt);
  SWAP(u0, u);
  diffuse(N, 1, u, u0, visc, dt);
  SWAP(v0, v);
  diffuse(N, 2, v, v0, visc, dt);
  project(N, u, v, u0, v0);
  SWAP(u0, u);
  SWAP(v0, v);
  advect(N, 1, u, u0, u0, v0, dt);
  advect(N, 2, v, v0, u0, v0, dt);
  project(N, u, v, u0, v0);
}

void man_step(int N, float *x, float *x0, float *u, float *v, float diff,
              float dt) {
  add_source(N, x, x0, dt);
}
