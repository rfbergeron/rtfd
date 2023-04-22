#ifndef RENDERER_H
#define RENDERER_H
#include <GLFW/glfw3.h>
#include <stdbool.h>

#include "glad/gl.h"

typedef struct renderer {
  GLFWwindow *window;
  float *vel_vertices, *dens_vertices;
  unsigned int *dens_triangles;
  unsigned int vel_shader, dens_shader;
  size_t N;
  unsigned int dens_ebo;
  unsigned int vel_vao, dens_vao;
  unsigned int vel_vbo, dens_vbo;
  double xpos, ypos, old_xpos, old_ypos;
  bool draw_vel, should_clear, add_dens, add_vel;
} Renderer;

Renderer *renderer_init(size_t N, int width, int height, const char *title);
void renderer_destroy(Renderer *renderer);
void renderer_draw(Renderer *renderer);
void renderer_update(Renderer *renderer, float *d, float *u, float *v);
void renderer_get_input(Renderer *renderer, float *d, float *u, float *v);
bool renderer_should_close(Renderer *renderer);
bool renderer_should_clear(Renderer *renderer);
void renderer_dump_dens(Renderer *renderer);
void renderer_dump_vel(Renderer *renderer);
#endif
