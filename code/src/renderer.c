#include "renderer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// density vertex contents do not have borders, but there is an extra row/col
#define DENS_IX(i, j) (3 * (renderer->sim_size + 1) * (i) + 3 * (j))
// triangle array contents do not have borders; iterated in pairs; the first
// triangle in the pair corresponds to the lower-left corner of the cell, and
// the second to the upper-right corner of the cell
#define TRIANGLE_IX(i, j) (6 * renderer->sim_size * (i) + 6 * (j))
// velocity array contents do not have borders; iterated in pairs; the first
// vertex should always be centered on the cell and the second should be moved
// based on the velocity
#define VEL_IX(i, j) (4 * renderer->sim_size * (i) + 4 * (j))
// extra row and column of vertices for top and left borders
#define DENS_VERTEX_COUNT ((renderer->sim_size + 1) * (renderer->sim_size + 1))
// 2 triangles per non-border cell
#define DENS_TRIANGLE_COUNT (2 * renderer->sim_size * renderer->sim_size)
// 2 vertices per line, one line per non-border cell
#define VEL_VERTEX_COUNT (2 * renderer->sim_size * renderer->sim_size)

#define LOG_BUFFER_SIZE 512
#define SHADER_BUFFER_SIZE 4096
#define NO_PROGRAM 0

static void generate_dens_triangles(Renderer *renderer) {
  // generate vertices; the order of the loops is reversed since j corresponds
  // to the y axis and we want the vertices to be laid out such that the first
  // row of the matrix corresponds to the first (bottom) row of vertices on the
  // screen. so, the layout of the vertices in this array should be such that
  // a logical row in the array corresponds to a row of vertices on the screen
  // from left to right. rows should start from the bottom of the screen and
  // move up to the top.
  for (size_t j = 0; j <= renderer->sim_size; ++j) {
    for (size_t i = 0; i <= renderer->sim_size; ++i) {
      renderer->dens_vertices[DENS_IX(i, j)] =
          2.0f * i / (float)renderer->sim_size - 1.0f;
      renderer->dens_vertices[DENS_IX(i, j) + 1] =
          2.0f * j / (float)renderer->sim_size - 1.0f;
      renderer->dens_vertices[DENS_IX(i, j) + 2] = 0;
    }
  }
  // generate triangle indices; the first triangle in the pair should be in the
  // bottom left corner, while the second should be in the top right
  for (size_t i = 0; i < renderer->sim_size; ++i) {
    for (size_t j = 0; j < renderer->sim_size; ++j) {
      renderer->dens_triangles[TRIANGLE_IX(i, j) + 1] =
          renderer->dens_triangles[TRIANGLE_IX(i, j) + 4] =
              DENS_IX(i + 1, j) / 3;
      renderer->dens_triangles[TRIANGLE_IX(i, j) + 2] =
          renderer->dens_triangles[TRIANGLE_IX(i, j) + 5] =
              DENS_IX(i, j + 1) / 3;
      renderer->dens_triangles[TRIANGLE_IX(i, j)] = DENS_IX(i, j) / 3;
      renderer->dens_triangles[TRIANGLE_IX(i, j) + 3] =
          DENS_IX(i + 1, j + 1) / 3;
    }
  }
}

static void generate_vel_lines(Renderer *renderer) {
  const float half_cell_width = 1.0f / renderer->sim_size;
  for (size_t i = 0; i < renderer->sim_size; ++i) {
    for (size_t j = 0; j < renderer->sim_size; ++j) {
      renderer->vel_vertices[VEL_IX(i, j)] =
          renderer->vel_vertices[VEL_IX(i, j) + 2] =
              2.0f * i / (float)renderer->sim_size - 1.0f + half_cell_width;
      renderer->vel_vertices[VEL_IX(i, j) + 1] =
          renderer->vel_vertices[VEL_IX(i, j) + 3] =
              2.0f * j / (float)renderer->sim_size - 1.0f + half_cell_width;
    }
  }
}

static void setup_dens_objects(Renderer *renderer) {
  static const size_t POSITION_VBO_OFFSET = 0, POSITION_ELEM_COUNT = 2;
  static const size_t COLOR_VBO_OFFSET = 2, COLOR_ELEM_COUNT = 1;
  static const size_t BUFFER_STRIDE = 3;
  static const unsigned int POSITION_ATTR = 0, COLOR_ATTR = 1;

  unsigned int buffers[2], vertex_array;
  glGenBuffers(2, buffers);
  glGenVertexArrays(1, &vertex_array);
  glBindVertexArray(vertex_array);
  glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
  glBufferData(GL_ARRAY_BUFFER, 3 * DENS_VERTEX_COUNT * sizeof(float),
               renderer->dens_vertices, GL_DYNAMIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[1]);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
               3 * DENS_TRIANGLE_COUNT * sizeof(unsigned int),
               renderer->dens_triangles, GL_STATIC_DRAW);
  glVertexAttribPointer(POSITION_ATTR, POSITION_ELEM_COUNT, GL_FLOAT, GL_FALSE,
                        BUFFER_STRIDE * sizeof(float),
                        (void *)(POSITION_VBO_OFFSET * sizeof(float)));
  glEnableVertexAttribArray(POSITION_ATTR);
  glVertexAttribPointer(COLOR_ATTR, COLOR_ELEM_COUNT, GL_FLOAT, GL_FALSE,
                        BUFFER_STRIDE * sizeof(float),
                        (void *)(COLOR_VBO_OFFSET * sizeof(float)));
  glEnableVertexAttribArray(COLOR_ATTR);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  renderer->dens_vbo = buffers[0];
  renderer->dens_ebo = buffers[1];
  renderer->dens_vao = vertex_array;
}

static void setup_vel_objects(Renderer *renderer) {
  static const size_t POSITION_VBO_OFFSET = 0, POSITION_ELEM_COUNT = 2;
  static const size_t BUFFER_STRIDE = 2;
  static const size_t POSITION_ATTR = 0;

  unsigned int vertex_buffer, vertex_array;
  glGenBuffers(1, &vertex_buffer);
  glGenVertexArrays(1, &vertex_array);
  glBindVertexArray(vertex_array);
  glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
  glBufferData(GL_ARRAY_BUFFER, 2 * VEL_VERTEX_COUNT * sizeof(float),
               renderer->vel_vertices, GL_DYNAMIC_DRAW);
  glVertexAttribPointer(POSITION_ATTR, POSITION_ELEM_COUNT, GL_FLOAT, GL_FALSE,
                        BUFFER_STRIDE * sizeof(float),
                        (void *)(POSITION_VBO_OFFSET * sizeof(float)));
  glEnableVertexAttribArray(POSITION_ATTR);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  glUseProgram(renderer->vel_shader);
  int colorLocation = glGetUniformLocation(renderer->vel_shader, "color");
  glUniform4f(colorLocation, 1.0f, 1.0f, 1.0f, 1.0f);
  glUseProgram(NO_PROGRAM);

  renderer->vel_vao = vertex_array;
  renderer->vel_vbo = vertex_buffer;
}

// TODO(Robert): error handling with gotos
unsigned int shader_init(const char *vertexPath, const char *fragmentPath) {
  int success;
  char infoLog[LOG_BUFFER_SIZE];
  char shaderBuffer[SHADER_BUFFER_SIZE];

  FILE *vertexFile = fopen(vertexPath, "r");
  if (vertexFile == NULL) {
    perror("Failed to open vertex shader source");
    return NO_PROGRAM;
  }

  shaderBuffer[0] = '\0';
  while (SHADER_BUFFER_SIZE - strlen(shaderBuffer) - 1 > 0 &&
         fgets(shaderBuffer + strlen(shaderBuffer),
               SHADER_BUFFER_SIZE - strlen(shaderBuffer), vertexFile))
    ;
  if (ferror(vertexFile)) {
    fprintf(stderr, "Error reading vertex shader source\n");
    fclose(vertexFile);
    return NO_PROGRAM;
  } else if (!feof(vertexFile)) {
    fprintf(stderr, "Buffer too small to hold vertex shader source\n");
    fclose(vertexFile);
    return NO_PROGRAM;
  } else {
    fclose(vertexFile);
  }

  unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
  const char *vertexShaderSource = shaderBuffer;
  glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
  glCompileShader(vertexShader);
  glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
  if (!success) {
    glGetShaderInfoLog(vertexShader, LOG_BUFFER_SIZE, NULL, infoLog);
    glDeleteShader(vertexShader);
    fprintf(stderr, "Failed to compile vertex shader\n%s\n", infoLog);
    return NO_PROGRAM;
  }

  FILE *fragmentFile = fopen(fragmentPath, "r");
  if (fragmentFile == NULL) {
    glDeleteShader(vertexShader);
    perror("Failed to open fragment shader source");
    return NO_PROGRAM;
  }

  shaderBuffer[0] = '\0';
  while (SHADER_BUFFER_SIZE - strlen(shaderBuffer) - 1 > 0 &&
         fgets(shaderBuffer + strlen(shaderBuffer),
               SHADER_BUFFER_SIZE - strlen(shaderBuffer), fragmentFile))
    ;
  if (ferror(fragmentFile)) {
    glDeleteShader(vertexShader);
    fprintf(stderr, "Error reading fragment shader source\n");
    fclose(fragmentFile);
    return NO_PROGRAM;
  } else if (!feof(fragmentFile)) {
    glDeleteShader(vertexShader);
    fprintf(stderr, "Buffer too small to hold fragment shader source\n");
    fclose(fragmentFile);
    return NO_PROGRAM;
  } else {
    fclose(fragmentFile);
  }

  unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
  const char *fragmentShaderSource = shaderBuffer;
  glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
  glCompileShader(fragmentShader);
  glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
  if (!success) {
    glDeleteShader(vertexShader);
    glGetShaderInfoLog(fragmentShader, LOG_BUFFER_SIZE, NULL, infoLog);
    glDeleteShader(fragmentShader);
    fprintf(stderr, "Failed to compile fragment shader\n%s\n", infoLog);
    return NO_PROGRAM;
  }

  unsigned int shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);
  glLinkProgram(shaderProgram);
  glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
  if (!success) {
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    glGetProgramInfoLog(shaderProgram, LOG_BUFFER_SIZE, NULL, infoLog);
    glDeleteProgram(shaderProgram);
    fprintf(stderr, "Failed to link shader program\n%s\n", infoLog);
    return NO_PROGRAM;
  }

  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);
  return shaderProgram;
}

static void framebuffer_size_callback(GLFWwindow *window, int width,
                                      int height) {
  (void)window;  // unused parameters
  glViewport(0, 0, width, height);
}

/*
static void key_callback(GLFWwindow *window, int key, int scancode, int action,
                         int mods) {
  (void)scancode, (void)mods;  // unused parameters
  Renderer *renderer = glfwGetWindowUserPointer(window);
  if (action != GLFW_PRESS) return;
  switch (key) {
    case GLFW_KEY_Q:
      glfwSetWindowShouldClose(window, GLFW_TRUE);
      break;
    case GLFW_KEY_C:
      renderer->should_clear = true;
      break;
    case GLFW_KEY_V:
      renderer->draw_vel = !renderer->draw_vel;
      break;
  }
}
*/

static void character_callback(GLFWwindow *window, unsigned int codepoint) {
  Renderer *renderer = glfwGetWindowUserPointer(window);
  switch (codepoint) {
    case 'q':
    case 'Q':
      glfwSetWindowShouldClose(window, GLFW_TRUE);
      break;
    case 'c':
    case 'C':
      renderer->should_clear = true;
      break;
    case 'v':
    case 'V':
      renderer->draw_vel = !renderer->draw_vel;
      break;
  }
}

static void cursor_position_callback(GLFWwindow *window, double xpos,
                                     double ypos) {
  Renderer *renderer = glfwGetWindowUserPointer(window);
  renderer->old_xpos = renderer->xpos, renderer->old_ypos = renderer->ypos;
  renderer->xpos = xpos, renderer->ypos = ypos;
}

static void mouse_button_callback(GLFWwindow *window, int button, int action,
                                  int mods) {
  (void)mods;  // unused parameters;
  Renderer *renderer = glfwGetWindowUserPointer(window);
  if (button == GLFW_MOUSE_BUTTON_LEFT) {
    if (action == GLFW_PRESS)
      renderer->add_vel = true;
    else if (action == GLFW_RELEASE)
      renderer->add_vel = false;
  } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
    if (action == GLFW_PRESS)
      renderer->add_dens = true;
    else if (action == GLFW_RELEASE)
      renderer->add_dens = false;
  }
}

Renderer *renderer_init(size_t sim_size, int width, int height,
                        const char *title) {
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
#endif

  Renderer *renderer = malloc(sizeof(Renderer));
  if (renderer == NULL) {
    fprintf(stderr, "Failed to allocate memory\n");
    glfwTerminate();
    return NULL;
  }

  renderer->window = glfwCreateWindow(width, height, title, NULL, NULL);
  if (renderer->window == NULL) {
    fprintf(stderr, "Failed to create GLFW window\n");
    free(renderer);
    glfwTerminate();
    return NULL;
  }
  glfwMakeContextCurrent(renderer->window);
  glfwSetFramebufferSizeCallback(renderer->window, framebuffer_size_callback);
  glfwSetCursorPosCallback(renderer->window, cursor_position_callback);
  glfwSetMouseButtonCallback(renderer->window, mouse_button_callback);
  // glfwSetKeyCallback(renderer->window, key_callback);
  glfwSetCharCallback(renderer->window, character_callback);
  glfwSetWindowUserPointer(renderer->window, renderer);

  if (!gladLoadGL((GLADloadfunc)glfwGetProcAddress)) {
    fprintf(stderr, "Failed to initialize GLAD\n");
    free(renderer);
    glfwTerminate();
    return NULL;
  }

  renderer->sim_size = sim_size;
  renderer->draw_vel = renderer->should_clear = false;
  renderer->add_dens = renderer->add_vel = false;

  renderer->vel_shader =
      shader_init("src/vel_shader.vert", "src/vel_shader.frag");
  if (renderer->vel_shader == NO_PROGRAM) {
    free(renderer);
    glfwTerminate();
    return NULL;
  }
  renderer->dens_shader =
      shader_init("src/dens_shader.vert", "src/dens_shader.frag");
  if (renderer->dens_shader == NO_PROGRAM) {
    glDeleteProgram(renderer->vel_shader);
    free(renderer);
    glfwTerminate();
    return NULL;
  }

  // 2 elements per vertex: x, y
  renderer->vel_vertices = malloc(sizeof(float) * 2 * VEL_VERTEX_COUNT);
  // 3 elements per vertex: x, y, and color
  renderer->dens_vertices = malloc(sizeof(float) * 3 * DENS_VERTEX_COUNT);
  // 3 indices per triangle
  renderer->dens_triangles =
      malloc(sizeof(unsigned int) * 3 * DENS_TRIANGLE_COUNT);
  if (renderer == NULL || renderer->vel_vertices == NULL ||
      renderer->dens_vertices == NULL || renderer->dens_triangles == NULL) {
    fprintf(stderr, "Failed to allocate memory\n");
    glDeleteProgram(renderer->vel_shader);
    glDeleteProgram(renderer->dens_shader);
    free(renderer->vel_vertices);
    free(renderer->dens_vertices);
    free(renderer->dens_triangles);
    free(renderer);
    glfwTerminate();
    return NULL;
  }

  generate_vel_lines(renderer);
  setup_vel_objects(renderer);
  generate_dens_triangles(renderer);
  setup_dens_objects(renderer);
  return renderer;
}

void renderer_destroy(Renderer *renderer) {
  glDeleteVertexArrays(1, &renderer->dens_vao);
  glDeleteVertexArrays(1, &renderer->vel_vao);
  glDeleteBuffers(1, &renderer->vel_vbo);
  glDeleteBuffers(1, &renderer->dens_vbo);
  glDeleteBuffers(1, &renderer->dens_ebo);
  glDeleteProgram(renderer->vel_shader);
  glDeleteProgram(renderer->dens_shader);
  free(renderer->vel_vertices);
  free(renderer->dens_vertices);
  free(renderer->dens_triangles);
  free(renderer);
  glfwTerminate();
}

void renderer_update(Renderer *renderer, Solver *solver) {
  // iterate density information separately to account for the extra row/col
  for (size_t i = 0; i <= renderer->sim_size; ++i) {
    for (size_t j = 0; j <= renderer->sim_size; ++j) {
      renderer->dens_vertices[DENS_IX(i, j) + 2] =
          *solver_at(solver, i, j, SLV_MAT_D);
    }
  }

  for (size_t i = 0; i < renderer->sim_size; ++i) {
    for (size_t j = 0; j < renderer->sim_size; ++j) {
      renderer->vel_vertices[VEL_IX(i, j) + 2] =
          renderer->vel_vertices[VEL_IX(i, j)] +
          *solver_at(solver, i, j, SLV_MAT_U);
      renderer->vel_vertices[VEL_IX(i, j) + 3] =
          renderer->vel_vertices[VEL_IX(i, j) + 1] +
          *solver_at(solver, i, j, SLV_MAT_V);
    }
  }

  glBindVertexArray(renderer->dens_vao);
  glBindBuffer(GL_ARRAY_BUFFER, renderer->dens_vbo);
  glBufferSubData(GL_ARRAY_BUFFER, 0, 3 * DENS_VERTEX_COUNT * sizeof(float),
                  renderer->dens_vertices);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glBindVertexArray(renderer->vel_vao);
  glBindBuffer(GL_ARRAY_BUFFER, renderer->vel_vbo);
  glBufferSubData(GL_ARRAY_BUFFER, 0, 2 * VEL_VERTEX_COUNT * sizeof(float),
                  renderer->vel_vertices);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

void renderer_draw(Renderer *renderer) {
  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
  if (renderer->draw_vel) {
    glUseProgram(renderer->vel_shader);
    glBindVertexArray(renderer->vel_vao);
    glDrawArrays(GL_LINES, 0, 2 * VEL_VERTEX_COUNT);
  } else {
    glUseProgram(renderer->dens_shader);
    glBindVertexArray(renderer->dens_vao);
    glDrawElements(GL_TRIANGLES, DENS_TRIANGLE_COUNT * 3, GL_UNSIGNED_INT, 0);
  }
  glUseProgram(NO_PROGRAM);
  glBindVertexArray(0);
  glfwSwapBuffers(renderer->window);
  glfwPollEvents();
}

void renderer_get_input(Renderer *renderer, Solver *solver) {
  int width, height;
  // get window size in screen coordinates; cursor position has been recorded
  // in screen coordinates relative to the upper-left corner of the window
  glfwGetWindowSize(renderer->window, &width, &height);

  if (renderer->xpos < 0 || renderer->xpos >= width || renderer->ypos < 0 ||
      renderer->ypos >= height)
    return;
  size_t x = renderer->xpos / width * renderer->sim_size;
  size_t y = (height - renderer->ypos) / height * renderer->sim_size;

  if (renderer->add_dens) {
    *solver_at(solver, x, y, SLV_MAT_D0) = 1.0f;
  }

  if (renderer->add_vel) {
    float u_force = (renderer->xpos - renderer->old_xpos);
    float v_force = (renderer->old_ypos - renderer->ypos);
    *solver_at(solver, x, y, SLV_MAT_U0) = u_force;
    *solver_at(solver, x, y, SLV_MAT_V0) = v_force;
  }
}

bool renderer_should_close(Renderer *renderer) {
  return glfwWindowShouldClose(renderer->window);
}

bool renderer_should_clear(Renderer *renderer) {
  return renderer->should_clear ? (renderer->should_clear = false, true)
                                : false;
}
