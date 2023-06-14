#ifndef CL_HELPER_H
#define CL_HELPER_H
#include <CL/cl.h>
#include <stdbool.h>

typedef struct _cl_bundle {
  cl_platform_id plat;
  cl_device_id dev;
  cl_context ctx;
  cl_command_queue h_cq;
  cl_command_queue d_cq;
  cl_kernel kernels[5];
  cl_mem d_buffers[4];
} *cl_bundle;

cl_bundle init_gpu_bundle(size_t sim_size, cl_int *status,
                          const char **errmsg_out);
cl_int free_bundle(cl_bundle bundle);
cl_int cl_solve_setup(cl_bundle bundle, size_t sim_size,
                      const float *restrict h_x, const float *restrict h_x0,
                      const char **errmsg_out);
cl_int cl_solve_retrieve(cl_bundle bundle, size_t sim_size, float *h_x,
                         const char **errmsg_out);
int cl_solve_step(cl_bundle bundle, unsigned int sim_size, float a, float c,
                  bool negate_axes[2], const char **errmsg_out);
cl_int cl_project_setup(cl_bundle bundle, size_t sim_size,
                        const float *restrict h_u, const float *restrict h_v,
                        const char **errmsg_out);
cl_int cl_project_retrieve(cl_bundle bundle, size_t sim_size,
                           float *restrict h_u, float *restrict h_v,
                           const char **errmsg_out);
cl_int cl_project_one(cl_bundle bundle, unsigned int sim_size,
                      const char **errmsg_out);
cl_int cl_project_two(cl_bundle bundle, unsigned int sim_size,
                      const char **errmsg_out);
cl_int cl_advect_setup(cl_bundle bundle, size_t sim_size,
                       const float *restrict h_x0, const float *restrict h_u,
                       const float *restrict h_v, const char **errmsg_out);
cl_int cl_advect_retrieve(cl_bundle bundle, size_t sim_size, float *h_x,
                          const char **errmsg_out);
cl_int cl_advect(cl_bundle bundle, unsigned int sim_size, float dt,
                 const bool negate_axes[2], const char **errmsg_out);
#endif
