#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "renderer.h"
#include "solver.h"

static size_t N = 64;
static float dt = 0.1f, diff, visc;
static float force = 5.0f, source = 100.0f;
static float *u, *v, *u_prev, *v_prev;
static float *dens, *dens_prev;

static const int WIN_WIDTH = 512, WIN_HEIGHT = 512;
static const char *WIN_TITLE = "Alias | wavefront";

static void free_data(void) {
  free(u);
  free(v);
  free(u_prev);
  free(v_prev);
  free(dens);
  free(dens_prev);
}

static void clear_data(void) {
  for (size_t i = 0, size = ACTUALSIZE; i < size; i++) {
    u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
  }
}

static int allocate_data(void) {
  size_t size = ACTUALSIZE;

  u = aligned_alloc(64, size * sizeof(float));
  v = aligned_alloc(64, size * sizeof(float));
  u_prev = aligned_alloc(64, size * sizeof(float));
  v_prev = aligned_alloc(64, size * sizeof(float));
  dens = aligned_alloc(64, size * sizeof(float));
  dens_prev = aligned_alloc(64, size * sizeof(float));

  if (!u || !v || !u_prev || !v_prev || !dens || !dens_prev) {
    fprintf(stderr, "cannot allocate data\n");
    return 1;
  }

  return 0;
}

int main(int argc, char **argv) {
  if (argc != 1 && argc != 6) {
    fprintf(stderr, "usage : %s N dt diff visc force source\n", argv[0]);
    fprintf(stderr, "where:\n");
    fprintf(stderr, "\t N      : grid resolution\n");
    fprintf(stderr, "\t dt     : time step\n");
    fprintf(stderr, "\t diff   : diffusion rate of the density\n");
    fprintf(stderr, "\t visc   : viscosity of the fluid\n");
    fprintf(stderr,
            "\t force  : scales the mouse movement that generate a force\n");
    fprintf(stderr, "\t source : amount of density that will be deposited\n");
    exit(EXIT_FAILURE);
  }

  if (argc == 1) {
    fprintf(
        stderr,
        "Using defaults : N=%zu dt=%g diff=%g visc=%g force = %g source=%g\n",
        N, dt, diff, visc, force, source);
  } else {
    N = atoi(argv[1]);
    dt = atof(argv[2]);
    diff = atof(argv[3]);
    visc = atof(argv[4]);
    force = atof(argv[5]);
    source = atof(argv[6]);
  }

  printf("\n\nHow to use this demo:\n\n");
  printf("\t Add densities with the right mouse button\n");
  printf(
      "\t Add velocities with the left mouse button and dragging the mouse\n");
  printf("\t Toggle density/velocity display with the 'v' key\n");
  printf("\t Clear the simulation by pressing the 'c' key\n");
  printf("\t Quit by pressing the 'q' key\n");

  if (allocate_data()) exit(EXIT_FAILURE);
  clear_data();

  Renderer *renderer = renderer_init(N, WIN_WIDTH, WIN_HEIGHT, WIN_TITLE);
  if (renderer == NULL) exit(EXIT_FAILURE);

  while (!renderer_should_close(renderer)) {
    renderer_get_input(renderer, dens_prev, u_prev, v_prev);
    if (renderer_should_clear(renderer)) clear_data();
    vel_step(N, u, v, u_prev, v_prev, visc, dt);
    dens_step(N, dens, dens_prev, u, v, diff, dt);
    renderer_update(renderer, dens, u, v);
    renderer_draw(renderer);
  }

  renderer_destroy(renderer);
  free_data();
  exit(EXIT_SUCCESS);
}
