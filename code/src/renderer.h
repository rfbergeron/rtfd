#ifndef RENDERER_H
#define RENDERER_H
#include <GLFW/glfw3.h>
#include <stdbool.h>

#include "common_solver.h"
#include "glad/gl.h"

typedef struct renderer {
  GLFWwindow *window;
  float *vel_vertices, *dens_vertices;
  unsigned int *dens_triangles;
  size_t sim_size;
  double xpos, ypos, old_xpos, old_ypos;
  unsigned int vel_shader, dens_shader;
  unsigned int dens_ebo;
  unsigned int vel_vao, dens_vao;
  unsigned int vel_vbo, dens_vbo;
  bool draw_vel, should_clear, add_dens, add_vel;
} Renderer;

Renderer *renderer_init(size_t sim_size, int width, int height,
                        const char *title);
void renderer_destroy(Renderer *renderer);
void renderer_draw(Renderer *renderer);
void renderer_update(Renderer *renderer, Solver *solver);
void renderer_get_input(Renderer *renderer, Solver *solver);
bool renderer_should_close(Renderer *renderer);
bool renderer_should_clear(Renderer *renderer);
void renderer_dump_dens(Renderer *renderer);
void renderer_dump_vel(Renderer *renderer);
#endif
