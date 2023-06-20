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
  cl_mem d_buffers[5];
} *cl_bundle;

cl_bundle init_gpu_bundle(size_t sim_size, cl_int *status,
                          const char **errmsg_out);
cl_int free_bundle(cl_bundle bundle);
cl_int cl_dens_step_full(cl_bundle bundle, const size_t sim_size,
                         float *restrict h_x, const float *restrict h_x0,
                         const float *restrict h_u, const float *restrict h_v,
                         const float diff, const float dt,
                         const size_t iterations, const char **errmsg_out);
cl_int cl_vel_step_full(cl_bundle bundle, const size_t sim_size,
                        float *restrict h_u, float *restrict h_v,
                        const float *restrict h_u0, const float *restrict h_v0,
                        const float visc, const float dt,
                        const size_t iterations, const char **errmsg_out);
#endif
