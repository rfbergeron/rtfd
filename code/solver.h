#ifndef SOLVER_H
#define SOLVER_H
#include <stddef.h>
#define IX(i, j) ((j) + (N + 2) * (i))
#define ACTUALSIZE ((N + 2) * (N + 2))

void dens_step(size_t N, float* x, float* x0, float* u, float* v, float diff,
               float dt);
void vel_step(size_t N, float* u, float* v, float* u0, float* v0, float visc,
              float dt);
#endif
