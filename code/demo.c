/*
 * ======================================================================
 * demo.c --- protoype to show off the simple solver
 * ----------------------------------------------------------------------
 * Author : Jos Stam (jstam@aw.sgi.com)
 * Creation Date : Jan 9 2003
 *
 * Description:
 *
 *      This code is a simple prototype that demonstrates how to use the
 *      code provided in my GDC2003 paper entitles "Real-Time Fluid Dynamics
 *      for Games". This code uses OpenGL and GLUT for graphics and interface
 *
 * =======================================================================
 */

#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include "auxlib.h"

// macros
#define BORDER 8
#define IX(i, j) ((i) + (N + BORDER) * (j))
#define ACTUALSIZE ((N + BORDER) * (N + BORDER))
#define BSTART (BORDER / 2) - 1
#define BEND (BSTART + N + 1)

// external definitions (from solver.c)
extern void dens_step(int N, float *x, float *x0, float *u, float *v,
                      float diff, float dt);
extern void man_step(int N, float *x, float *x0, float *u, float *v,
                      float diff, float dt);
extern void vel_step(int N, float *u, float *v, float *u0, float *v0,
                     float visc, float dt);

// global variables
static int N, dvel, manual_step;
static float dt, diff, visc, force, source;

static float *u, *v, *u_prev, *v_prev;
static float *dens, *dens_prev;

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int omx, omy, mx, my;

/*
 * ----------------------------------------------------------------------
 * free/clear/allocate simulation data
 * ----------------------------------------------------------------------
 */

static void free_data(void) {
  if (u) free(u);
  if (v) free(v);
  if (u_prev) free(u_prev);
  if (v_prev) free(v_prev);
  if (dens) free(dens);
  if (dens_prev) free(dens_prev);
}

static void clear_data(void) {
  int i, size = ACTUALSIZE;

  for (i = 0; i < size; i++) {
    u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
  }
}

static int allocate_data(void) {
  // allocate a full 4 extra elements in each direction
  int size = ACTUALSIZE;

  u = aligned_alloc(64, size * sizeof(float));
  v = aligned_alloc(64, size * sizeof(float));
  u_prev = aligned_alloc(64, size * sizeof(float));
  v_prev = aligned_alloc(64, size * sizeof(float));
  dens = aligned_alloc(64, size * sizeof(float));
  dens_prev = aligned_alloc(64, size * sizeof(float));

  if (!u || !v || !u_prev || !v_prev || !dens || !dens_prev) {
    fprintf(stderr, "cannot allocate data\n");
    return (0);
  }

  return (1);
}

/*
 * ----------------------------------------------------------------------
 * OpenGL specific drawing routines
 * ----------------------------------------------------------------------
 */

static void pre_display(void) {
  glViewport(0, 0, win_x, win_y);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0, 1.0, 0.0, 1.0);
  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
}

static void post_display(void) { glutSwapBuffers(); }

static void draw_velocity(void) {
  int i, j;
  float x, y, h;

  h = 1.0f / N;

  glColor3f(1.0f, 1.0f, 1.0f);
  glLineWidth(1.0f);

  glBegin(GL_LINES);

  for (i = BSTART+1; i < BEND; i++) {
    x = (i - 0.5f) * h;
    for (j = BSTART+1; j < BEND; j++) {
      y = (j - 0.5f) * h;

      glVertex2f(x, y);
      glVertex2f(x + u[IX(i, j)], y + v[IX(i, j)]);
    }
  }

  glEnd();
}

static void draw_density(void) {
  int i, j;
  float x, y, h, d00, d01, d10, d11;

  h = 1.0f / N;

  glBegin(GL_QUADS);

  for (i = BSTART; i < BEND; i++) {
    x = (i - BSTART + 1) * h;
    for (j = BSTART; j < BEND; j++) {
      y = (j - BSTART + 1) * h;

      d00 = dens[IX(i, j)];
      d01 = dens[IX(i, j + 1)];
      d10 = dens[IX(i + 1, j)];
      d11 = dens[IX(i + 1, j + 1)];

      glColor3f(d00, d00, d00);
      glVertex2f(x, y);
      glColor3f(d10, d10, d10);
      glVertex2f(x + h, y);
      glColor3f(d11, d11, d11);
      glVertex2f(x + h, y + h);
      glColor3f(d01, d01, d01);
      glVertex2f(x, y + h);
    }
  }

  glEnd();
}

/*
 * ----------------------------------------------------------------------
 * relates mouse movements to forces sources
 * ----------------------------------------------------------------------
 */

static void get_from_UI(float *d, float *u, float *v) {
  int i, j, size = ACTUALSIZE;

  for (i = 0; i < size; i++) {
    u[i] = v[i] = d[i] = 0.0f;
  }

  if (!mouse_down[0] && !mouse_down[2]) return;

  // why N + 1?
  i = (int)((mx / (float)win_x) * N + 1);
  j = (int)(((win_y - my) / (float)win_y) * N + 1);

  if (i < 1 || i > N || j < 1 || j > N) return;

  if (mouse_down[0]) {
    u[IX(i, j)] = force * (mx - omx);
    v[IX(i, j)] = force * (omy - my);
  }

  if (mouse_down[2]) {
    d[IX(i, j)] = source;
  }

  omx = mx;
  omy = my;

  return;
}

/*
 * ----------------------------------------------------------------------
 * GLUT callback routines
 * ----------------------------------------------------------------------
 */

static void key_func(unsigned char key, int x, int y) {
  switch (key) {
    case 'c':
    case 'C':
      clear_data();
      break;

    case 'q':
    case 'Q':
      free_data();
      exit(0);
      break;

    case 'v':
    case 'V':
      dvel = !dvel;
      break;

    case 'm':
    case 'M':
      manual_step = !manual_step;
      break;

    case 0x20:
      if(manual_step) {
        vel_step(N, u, v, u_prev, v_prev, visc, dt);
        dens_step(N, dens, dens_prev, u, v, diff, dt);
      }
      break;
  }
}

static void mouse_func(int button, int state, int x, int y) {
  omx = mx = x;
  omx = my = y;

  mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func(int x, int y) {
  mx = x;
  my = y;
}

static void reshape_func(int width, int height) {
  glutSetWindow(win_id);
  glutReshapeWindow(width, height);

  win_x = width;
  win_y = height;
}

static void idle_func(void) {
  get_from_UI(dens_prev, u_prev, v_prev);

  if(manual_step) {
    man_step(N, dens, dens_prev, u, v, diff, dt);
  } else {
    vel_step(N, u, v, u_prev, v_prev, visc, dt);
    dens_step(N, dens, dens_prev, u, v, diff, dt);
  }

  glutSetWindow(win_id);
  glutPostRedisplay();
}

static void display_func(void) {
  pre_display();

  if (dvel)
    draw_velocity();
  else
    draw_density();

  post_display();
}

/*
 * ----------------------------------------------------------------------
 * open_glut_window --- open a glut compatible window and set callbacks
 * ----------------------------------------------------------------------
 */

static void open_glut_window(void) {
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

  glutInitWindowPosition(0, 0);
  glutInitWindowSize(win_x, win_y);
  win_id = glutCreateWindow("Alias | wavefront");

  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
  glutSwapBuffers();
  glClear(GL_COLOR_BUFFER_BIT);
  glutSwapBuffers();

  pre_display();

  glutKeyboardFunc(key_func);
  glutMouseFunc(mouse_func);
  glutMotionFunc(motion_func);
  glutReshapeFunc(reshape_func);
  glutIdleFunc(idle_func);
  glutDisplayFunc(display_func);
}

void scan_options(int argc, char **argv) {
  opterr = 0;
  for (;;) {
    int option = getopt(argc, argv, "@:v:d:f:s:N:t:");

    if (option == EOF) break;
    switch (option) {
      case '@':
        debug_set_flags(optarg);
        break;
      case 'd':
        diff = atof(optarg);
        break;
      case 'v':
        visc = atof(optarg);
        break;
      case 'f':
        force = atof(optarg);
        break;
      case 's':
        source = atof(optarg);
        break;
      case 't':
        dt = atof(optarg);
        break;
      case 'N':
        N = atoi(optarg);
        break;
      case 'm':
      case 'M':
        manual_step = 1;
        break;
      default:
        fprintf(stderr, "-%c: invalid option\n", (char)optopt);
        break;
    }
  }
}

/*
 * ----------------------------------------------------------------------
 * main --- main routine
 * ----------------------------------------------------------------------
 */

int main(int argc, char **argv) {
  glutInit(&argc, argv);

  /* NOTE: because we are using aligned_alloc, we must have that (N+2)*(N+2)
   * is a multiple of 16
   * NOTE: trying out 64-byte alignment to see if this makes more efficient
   * use of the cache. It may me okay to only enforce that N is even, since
   * float is 4 bytes wide. This would mean that (N+2)^2 * sizeof(float) is
   * already a multiple of 16
   *
   * Actually, change of plans: N must be a multiple of 4; this should be
   * sufficient to achieve a total byte count that is a multiple of 64 with
   * 32-bit floats. Allocate 4 extra bytes on each size of the actual matrix
   * so that the vector scheme mirrors the scalar one.
   */
  N = 64;
  dvel = 0;
  manual_step = 0;
  dt = 0.1f;
  diff = 0.0f;
  visc = 0.0f;
  force = 5.0f;
  source = 100.0f;
  scan_options(argc, argv);
  printf("Using values : N=%d dt=%g diff=%g visc=%g force = %g source=%g\n", N,
         dt, diff, visc, force, source);

  printf("\n\nHow to use this demo:\n\n");
  printf("\t Add densities with the right mouse button\n");
  printf(
      "\t Add velocities with the left mouse button and dragging the mouse\n");
  printf("\t Toggle density/velocity display with the 'v' key\n");
  printf("\t Clear the simulation by pressing the 'c' key\n");
  printf("\t Quit by pressing the 'q' key\n");

  if (!allocate_data()) exit(1);
  DEBUGL('m', "data allocated");
  clear_data();
  DEBUGL('m', "data cleared");

  win_x = 512;
  win_y = 512;
  open_glut_window();

  DEBUGL('m', "entering main loop");
  glutMainLoop();

  exit(0);
}
