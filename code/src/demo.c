#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "renderer.h"
#include "solver.h"

static size_t sim_size = 512;
static float dt = 0.1f, diff = 0.0f, visc = 0.0f;
static float force = 5.0f, source = 100.0f;

static const int WIN_WIDTH = 512, WIN_HEIGHT = 512;
static const char *WIN_TITLE = "Alias | wavefront";

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
        sim_size, dt, diff, visc, force, source);
  } else {
    sim_size = atoi(argv[1]);
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

  Solver *solver = solver_init(sim_size, dt, diff, visc, force, source);
  if (!solver) exit(EXIT_FAILURE);
  solver_clear(solver, SLV_MAT_ALL);

  Renderer *renderer =
      renderer_init(sim_size, WIN_WIDTH, WIN_HEIGHT, WIN_TITLE);
  if (renderer == NULL) {
    solver_destroy(solver);
    exit(EXIT_FAILURE);
  }

  while (!renderer_should_close(renderer)) {
    if (renderer_should_clear(renderer)) solver_clear(solver, SLV_MAT_CURR);
    solver_clear(solver, SLV_MAT_U0 | SLV_MAT_V0);
    renderer_get_input(renderer, solver);
    solver_vel_step(solver);
    // d_prev is used as scratch space for jacobi solvers; clear it here
    solver_clear(solver, SLV_MAT_D0);
    renderer_get_input(renderer, solver);
    solver_dens_step(solver);
    renderer_update(renderer, solver);
    renderer_draw(renderer);
  }

  renderer_destroy(renderer);
  solver_destroy(solver);
  exit(EXIT_SUCCESS);
}
