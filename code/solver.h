#ifndef SOLVER_H
#define SOLVER_H
#include <stddef.h>
#define COLBEGIN 4
#define COLEND (N + COLBEGIN)
#define ROWBEGIN 1
#define ROWEND (N + ROWBEGIN)
#define ROWSIZE (N + ROWBEGIN * 2)
#define COLSIZE (N + COLBEGIN * 2)
#define IX(i, j) ((N + 2 * COLBEGIN) * (i) + (j))
#define ACTUALSIZE ((N + 2 * COLBEGIN) * (N + 2 * ROWBEGIN))

void dens_step(size_t N, float* x, float* x0, float* u, float* v, float diff,
               float dt);
void vel_step(size_t N, float* u, float* v, float* u0, float* v0, float visc,
              float dt);
#endif
