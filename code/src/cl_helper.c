#include "cl_helper.h"

#include <CL/cl.h>
#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "cl_common.h"

static const char *OPTS_FMT = "-I%s/src -cl-std=CL2.0";
static const char *KERNEL_NAMES[] = {"jacobi", "set_bnd", "project_one",
                                     "project_two"};
#define JACOBI_IX 0
#define SET_BND_IX 1
#define PROJECT_ONE_IX 2
#define PROJECT_TWO_IX 3
#define L_SIZE 16
#define L_ACTUAL_SIZE \
  ((P_SIZE * L_SIZE + 2 * COL_BORDER) * (L_SIZE + 2 * ROW_BORDER))
#define FAIL(msg, label) \
  do {                   \
    errmsg = msg;        \
    goto label;          \
  } while (0)

// https://stackoverflow.com/questions/24326432/convenient-way-to-show-opencl-error-codes
static const char *cl_strerror(cl_int error) {
  switch (error) {
    // run-time and JIT compiler errors
    case 0:
      return "CL_SUCCESS";
    case -1:
      return "CL_DEVICE_NOT_FOUND";
    case -2:
      return "CL_DEVICE_NOT_AVAILABLE";
    case -3:
      return "CL_COMPILER_NOT_AVAILABLE";
    case -4:
      return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case -5:
      return "CL_OUT_OF_RESOURCES";
    case -6:
      return "CL_OUT_OF_HOST_MEMORY";
    case -7:
      return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case -8:
      return "CL_MEM_COPY_OVERLAP";
    case -9:
      return "CL_IMAGE_FORMAT_MISMATCH";
    case -10:
      return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case -11:
      return "CL_BUILD_PROGRAM_FAILURE";
    case -12:
      return "CL_MAP_FAILURE";
    case -13:
      return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
    case -14:
      return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
    case -15:
      return "CL_COMPILE_PROGRAM_FAILURE";
    case -16:
      return "CL_LINKER_NOT_AVAILABLE";
    case -17:
      return "CL_LINK_PROGRAM_FAILURE";
    case -18:
      return "CL_DEVICE_PARTITION_FAILED";
    case -19:
      return "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";

    // compile-time errors
    case -30:
      return "CL_INVALID_VALUE";
    case -31:
      return "CL_INVALID_DEVICE_TYPE";
    case -32:
      return "CL_INVALID_PLATFORM";
    case -33:
      return "CL_INVALID_DEVICE";
    case -34:
      return "CL_INVALID_CONTEXT";
    case -35:
      return "CL_INVALID_QUEUE_PROPERTIES";
    case -36:
      return "CL_INVALID_COMMAND_QUEUE";
    case -37:
      return "CL_INVALID_HOST_PTR";
    case -38:
      return "CL_INVALID_MEM_OBJECT";
    case -39:
      return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case -40:
      return "CL_INVALID_IMAGE_SIZE";
    case -41:
      return "CL_INVALID_SAMPLER";
    case -42:
      return "CL_INVALID_BINARY";
    case -43:
      return "CL_INVALID_BUILD_OPTIONS";
    case -44:
      return "CL_INVALID_PROGRAM";
    case -45:
      return "CL_INVALID_PROGRAM_EXECUTABLE";
    case -46:
      return "CL_INVALID_KERNEL_NAME";
    case -47:
      return "CL_INVALID_KERNEL_DEFINITION";
    case -48:
      return "CL_INVALID_KERNEL";
    case -49:
      return "CL_INVALID_ARG_INDEX";
    case -50:
      return "CL_INVALID_ARG_VALUE";
    case -51:
      return "CL_INVALID_ARG_SIZE";
    case -52:
      return "CL_INVALID_KERNEL_ARGS";
    case -53:
      return "CL_INVALID_WORK_DIMENSION";
    case -54:
      return "CL_INVALID_WORK_GROUP_SIZE";
    case -55:
      return "CL_INVALID_WORK_ITEM_SIZE";
    case -56:
      return "CL_INVALID_GLOBAL_OFFSET";
    case -57:
      return "CL_INVALID_EVENT_WAIT_LIST";
    case -58:
      return "CL_INVALID_EVENT";
    case -59:
      return "CL_INVALID_OPERATION";
    case -60:
      return "CL_INVALID_GL_OBJECT";
    case -61:
      return "CL_INVALID_BUFFER_SIZE";
    case -62:
      return "CL_INVALID_MIP_LEVEL";
    case -63:
      return "CL_INVALID_GLOBAL_WORK_SIZE";
    case -64:
      return "CL_INVALID_PROPERTY";
    case -65:
      return "CL_INVALID_IMAGE_DESCRIPTOR";
    case -66:
      return "CL_INVALID_COMPILER_OPTIONS";
    case -67:
      return "CL_INVALID_LINKER_OPTIONS";
    case -68:
      return "CL_INVALID_DEVICE_PARTITION_COUNT";

    // extension errors
    case -1000:
      return "CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR";
    case -1001:
      return "CL_PLATFORM_NOT_FOUND_KHR";
    case -1002:
      return "CL_INVALID_D3D10_DEVICE_KHR";
    case -1003:
      return "CL_INVALID_D3D10_RESOURCE_KHR";
    case -1004:
      return "CL_D3D10_RESOURCE_ALREADY_ACQUIRED_KHR";
    case -1005:
      return "CL_D3D10_RESOURCE_NOT_ACQUIRED_KHR";
    default:
      return "Unknown OpenCL error";
  }
}

static cl_int init_kernel(cl_bundle bundle, size_t kern_ix,
                          const char **errmsg_out) {
  static const char *FILENAME = "src/solver.cl";
  const char *errmsg = NULL;
  bundle->kernels[kern_ix] = NULL;
  cl_int status = CL_SUCCESS;

  FILE *fp = fopen(FILENAME, "r");
  if (!fp) FAIL(strerror(errno), fail_open);

  char buffer[32768];
  size_t len = fread(buffer, sizeof(char), 32768, fp);
  if (ferror(fp))
    FAIL((status = -1, strerror(errno)), fail_read);
  else if (!feof(fp))
    FAIL((status = -1,
          "Unable to read OpenCL proram source: buffer size exceeded.\n"),
         fail_read);
  else
    buffer[len] = '\0';

  const char *buffer_ptr = buffer;
  cl_program prog =
      clCreateProgramWithSource(bundle->ctx, 1, &buffer_ptr, NULL, &status);
  if (status != CL_SUCCESS)
    FAIL("Unable to create OpenCL program with code: %s.\n", fail_prog);

  char cwd[4096];
  if (getcwd(cwd, 4096) == NULL)
    FAIL("Unable to get current working directory.\n", fail_getcwd);
  size_t opts_len = snprintf(NULL, 0, OPTS_FMT, cwd);
  char *opts = malloc(sizeof(char) * (opts_len + 1));
  (void)snprintf(opts, opts_len + 1, OPTS_FMT, cwd);

  status = clBuildProgram(prog, 0, NULL, opts, NULL, NULL);
  if (status != CL_SUCCESS) {
    clGetProgramBuildInfo(prog, bundle->dev, CL_PROGRAM_BUILD_LOG,
                          sizeof(buffer), buffer, &len);
    fprintf(stderr, "%s\n", buffer);
    FAIL("Unable to build OpenCL program with code: %s.\n", fail_build);
  }

  bundle->kernels[kern_ix] =
      clCreateKernel(prog, KERNEL_NAMES[kern_ix], &status);
  if (status != CL_SUCCESS)
    FAIL("Failed to create kernel with code: %s.\n", fail_kernel);
fail_kernel:
fail_build:
fail_getcwd:
  // safe to do here: program will not be deleted until kernel refcount is zero
  clReleaseProgram(prog);
fail_prog:
fail_read:
  fclose(fp);
fail_open:
  *errmsg_out = errmsg;
  return status;
}

static int init_buffers(cl_bundle bundle, size_t sim_size, cl_int *status,
                        const char **errmsg_out) {
  static const size_t BUFFER_COUNT =
      sizeof(bundle->d_buffers) / sizeof(bundle->d_buffers[0]);
  const char *errmsg = NULL;
  size_t i;
  for (i = 0; i < BUFFER_COUNT; ++i) {
    bundle->d_buffers[i] =
        clCreateBuffer(bundle->ctx, CL_MEM_READ_WRITE,
                       ACTUAL_SIZE * sizeof(float), NULL, status);
    if (*status != CL_SUCCESS)
      FAIL("Failed to create OpenCL buffer with code: %s.\n", fail);
    ;
  }
  *errmsg_out = errmsg;
  return 0;
fail:
  *errmsg_out = errmsg;
  size_t j = 0;
  for (; j < i; ++j) {
    cl_int fail_status = clReleaseMemObject(bundle->d_buffers[j]);
    if (fail_status != CL_SUCCESS)
      fprintf(stderr, "An error occurred while handling a previous error: %s\n",
              cl_strerror(fail_status));
    bundle->d_buffers[j] = NULL;
  }
  for (; j < BUFFER_COUNT; ++j) bundle->d_buffers[j] = NULL;
  return -1;
}

cl_bundle init_gpu_bundle(size_t sim_size, cl_int *status,
                          const char **errmsg_out) {
  const char *errmsg = NULL;
  cl_bundle ret = malloc(sizeof(*ret));
  if (ret == NULL)
    FAIL((*status = CL_SUCCESS, "Unable to allocate memory for cl_bundle.\n"),
         fail_alloc);
  ret->d_cq = NULL;
  cl_platform_id platforms[4];
  unsigned int num_platforms;
  *status = clGetPlatformIDs(4, platforms, &num_platforms);
  if (*status != CL_SUCCESS)
    FAIL("Unable to retrieve OpenCL platform ids with code: %s.\n", fail_plat);
  else
    printf("System has %u OpenCL platforms.\n", num_platforms);

  for (size_t i = 0; i < 4 && i < num_platforms; ++i) {
    *status =
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_GPU, 1, &ret->dev, NULL);
    if (*status == CL_DEVICE_NOT_FOUND)
      continue;
    else if (*status != CL_SUCCESS)
      FAIL("Unable to retrieve OpenCL device ids with code: %s.\n", fail_dev);
    else {
      ret->plat = platforms[i];
      break;
    }
  }

  if (ret->dev == NULL) FAIL("Unable to locate GPU with code: %s.\n", fail_dev);

  ret->ctx = clCreateContext(NULL, 1, &ret->dev, NULL, NULL, status);
  if (*status != CL_SUCCESS)
    FAIL("Unable to create OpenCL context with code: %s.\n", fail_ctx);

  ret->h_cq =
      clCreateCommandQueueWithProperties(ret->ctx, ret->dev, NULL, status);
  if (*status != CL_SUCCESS)
    FAIL("Unable to create OpenCL host queue with code: %s.\n", fail_h_cq);

  const cl_queue_properties d_cq_props[] = {
      CL_QUEUE_PROPERTIES,
      CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_ON_DEVICE |
          CL_QUEUE_ON_DEVICE_DEFAULT,
      0};
  ret->d_cq = clCreateCommandQueueWithProperties(ret->ctx, ret->dev, d_cq_props,
                                                 status);
  if (*status != CL_SUCCESS)
    FAIL("Unable to create OpenCL device queue with code: %s.\n", fail_d_cq);
  *status = init_kernel(ret, JACOBI_IX, errmsg_out);
  if (*status != CL_SUCCESS) FAIL(*errmsg_out, fail_jacobi);
  *status = init_kernel(ret, SET_BND_IX, errmsg_out);
  if (*status != CL_SUCCESS) FAIL(*errmsg_out, fail_set_bnd);
  *status = init_kernel(ret, PROJECT_ONE_IX, errmsg_out);
  if (*status != CL_SUCCESS) FAIL(*errmsg_out, fail_set_bnd);
  *status = init_kernel(ret, PROJECT_TWO_IX, errmsg_out);
  if (*status != CL_SUCCESS) FAIL(*errmsg_out, fail_set_bnd);
  *status = init_buffers(ret, sim_size, status, errmsg_out);
  if (*status != CL_SUCCESS) FAIL(*errmsg_out, fail_buffers);
  *errmsg_out = errmsg;
  return ret;
  cl_int unload_status;
fail_buffers:
  clReleaseKernel(ret->kernels[SET_BND_IX]);
fail_set_bnd:
  clReleaseKernel(ret->kernels[JACOBI_IX]);
fail_jacobi:
fail_d_cq:
  clReleaseCommandQueue(ret->h_cq);
fail_h_cq:
  clReleaseContext(ret->ctx);
fail_ctx:
  unload_status = clUnloadPlatformCompiler(ret->plat);
  if (unload_status != CL_SUCCESS)
    fprintf(stderr, "Unable to unload platform compiler with code: %s.\n",
            cl_strerror(unload_status));
fail_dev:
fail_plat:
  free(ret);
fail_alloc:
  *errmsg_out = errmsg;
  return NULL;
}

cl_bundle init_cpu_bundle(cl_int *status, const char **errmsg_out) {
  const char *errmsg = NULL;
  cl_bundle ret = malloc(sizeof(*ret));
  if (ret == NULL)
    FAIL((*status = CL_SUCCESS, "Unable to allocate memory for cl_bundle.\n"),
         fail_alloc);
  ret->d_cq = NULL;
  cl_platform_id platforms[4];
  unsigned int num_platforms;
  *status = clGetPlatformIDs(4, platforms, &num_platforms);
  if (*status != CL_SUCCESS)
    FAIL("Unable to retrieve OpenCL platform ids with code: %s.\n", fail_plat);
  else
    printf("System has %u OpenCL platforms.\n", num_platforms);

  for (size_t i = 0; i < 4 && i < num_platforms; ++i) {
    *status =
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_CPU, 1, &ret->dev, NULL);
    if (*status == CL_DEVICE_NOT_FOUND)
      continue;
    else if (*status != CL_SUCCESS)
      FAIL("Unable to retrieve OpenCL device ids with code: %s.\n", fail_dev);
    else {
      ret->plat = platforms[i];
      break;
    }
  }

  if (ret->dev == NULL)
    FAIL("Unable to locate pocl with code: %s.\n", fail_dev);

  ret->ctx = clCreateContext(NULL, 1, &ret->dev, NULL, NULL, status);
  if (*status != CL_SUCCESS)
    FAIL("Unable to create OpenCL context with code: %s.\n", fail_ctx);

  ret->h_cq =
      clCreateCommandQueueWithProperties(ret->ctx, ret->dev, NULL, status);
  if (*status != CL_SUCCESS)
    FAIL("Unable to create OpenCL host queue with code: %s.\n", fail_h_cq);

  /*
  const cl_queue_properties d_cq_props[] = {CL_QUEUE_PROPERTIES,
      CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_ON_DEVICE
          | CL_QUEUE_ON_DEVICE_DEFAULT, 0};
  ret->d_cq = clCreateCommandQueueWithProperties(ret->ctx, ret->dev, d_cq_props,
  status); if (*status != CL_SUCCESS) FAIL("Unable to create OpenCL device queue
  with code: %s.\n", fail_d_cq);
  */
  *errmsg_out = errmsg;
  return ret;
  cl_int unload_status;
/*
fail_d_cq:
  clReleaseCommandQueue(ret->h_cq);
*/
fail_h_cq:
  clReleaseContext(ret->ctx);
fail_ctx:
  unload_status = clUnloadPlatformCompiler(ret->plat);
  if (unload_status != CL_SUCCESS)
    fprintf(stderr, "Unable to unload platform compiler with code: %s.\n",
            cl_strerror(unload_status));
fail_dev:
fail_plat:
  free(ret);
fail_alloc:
  *errmsg_out = errmsg;
  return NULL;
}

cl_int free_bundle(cl_bundle bundle) {
  cl_int status;
  if (!bundle) return CL_SUCCESS;
  for (size_t i = 0;
       i < sizeof(bundle->d_buffers) / sizeof(bundle->d_buffers[0]); ++i)
    if (bundle->d_buffers[i] &&
        (status = clReleaseMemObject(bundle->d_buffers[i])))
      return status;
  for (size_t i = 0; i < sizeof(bundle->kernels) / sizeof(bundle->kernels[0]);
       ++i)
    if (bundle->kernels[i] && (status = clReleaseKernel(bundle->kernels[i])))
      return status;
  if (bundle->d_cq && (status = clReleaseCommandQueue(bundle->d_cq)))
    return status;
  if (bundle->h_cq && (status = clReleaseCommandQueue(bundle->h_cq)))
    return status;
  if (bundle->ctx && (status = clReleaseContext(bundle->ctx))) return status;
  if (bundle->plat && (status = clUnloadPlatformCompiler(bundle->plat)))
    return status;
  free(bundle);
  return CL_SUCCESS;
}

static cl_int cl_setup(cl_bundle bundle, size_t sim_size,
                       const float *restrict h_x, const float *restrict h_x0,
                       const float *restrict h_x1, const char **errmsg_out) {
  cl_int status = CL_SUCCESS;
  const char *errmsg = NULL;
  if (h_x) {
    status =
        clEnqueueWriteBuffer(bundle->h_cq, bundle->d_buffers[0], CL_TRUE, 0,
                             sizeof(float) * ACTUAL_SIZE, h_x, 0, NULL, NULL);
    if (status != CL_SUCCESS)
      FAIL("Unable to write host buffer with code %s.\n", fail);
  }
  if (h_x0) {
    status =
        clEnqueueWriteBuffer(bundle->h_cq, bundle->d_buffers[1], CL_TRUE, 0,
                             sizeof(float) * ACTUAL_SIZE, h_x0, 0, NULL, NULL);
    if (status != CL_SUCCESS)
      FAIL("Unable to write host buffer with code %s.\n", fail);
  }
  if (h_x1) {
    status =
        clEnqueueWriteBuffer(bundle->h_cq, bundle->d_buffers[2], CL_TRUE, 0,
                             sizeof(float) * ACTUAL_SIZE, h_x1, 0, NULL, NULL);
    if (status != CL_SUCCESS)
      FAIL("Unable to write host buffer with code %s.\n", fail);
  }
fail:
  *errmsg_out = errmsg;
  return status;
}

static cl_int cl_retrieve(cl_bundle bundle, size_t sim_size,
                          float *restrict h_x, float *restrict h_x0,
                          float *restrict h_x1, const char **errmsg_out) {
  cl_int status = CL_SUCCESS;
  const char *errmsg = NULL;
  if (h_x) {
    status =
        clEnqueueReadBuffer(bundle->h_cq, bundle->d_buffers[0], CL_TRUE, 0,
                            sizeof(float) * ACTUAL_SIZE, h_x, 0, NULL, NULL);
    if (status != CL_SUCCESS)
      FAIL("Unable to write host buffer with code %s.\n", fail);
  }
  if (h_x0) {
    status =
        clEnqueueReadBuffer(bundle->h_cq, bundle->d_buffers[1], CL_TRUE, 0,
                            sizeof(float) * ACTUAL_SIZE, h_x0, 0, NULL, NULL);
    if (status != CL_SUCCESS)
      FAIL("Unable to write host buffer with code %s.\n", fail);
  }
  if (h_x1) {
    status =
        clEnqueueReadBuffer(bundle->h_cq, bundle->d_buffers[2], CL_TRUE, 0,
                            sizeof(float) * ACTUAL_SIZE, h_x1, 0, NULL, NULL);
    if (status != CL_SUCCESS)
      FAIL("Unable to write host buffer with code %s.\n", fail);
  }
fail:
  *errmsg_out = errmsg;
  return status;
}

/*
static float d_full_reduce(cl_command_queue h_cq, cl_kernel kern,
                           unsigned int N, buffer_pair buffers[2],
                           cl_int *status_out, const char **errmsg_out) {
  const char *errmsg;
  size_t g_size = N_TO_G_SIZE(N);
  static const size_t l_size = LOCAL_SIZE, p_size = PRIVATE_SIZE;
  size_t num_groups = final_size(N);
  fprintf(stderr,
          "GLOBAL SIZE: %zu\nLOCAL SIZE: %zu\nPRIVATE SIZE: %zu\nN: %u\n"
          "FIRST ITERATION GROUP COUNT: %zu\n"
          "FINAL GROUP COUNT AFTER FULL REDUCTION: %zu\n",
          g_size, l_size, p_size, N, g_size / l_size, num_groups);

  if (num_groups == N) {
    fprintf(stderr, "Job too small; reducing on host.\n");
    *status_out = CL_SUCCESS;
    return h_reduce(N, buffers[0].host);
  }

  *status_out = clSetKernelArg(kern, 0, sizeof(unsigned int), &N);
  if (*status_out != CL_SUCCESS)
    FAIL("Unable to set kernel argument 0 with code: %s\n", fail);
  *status_out = clSetKernelArg(kern, 1, sizeof(cl_mem), &buffers[0].dev);
  if (*status_out != CL_SUCCESS)
    FAIL("Unable to set kernel argument 1 with code: %s\n", fail);
  *status_out = clSetKernelArg(kern, 2, sizeof(cl_mem), &buffers[1].dev);
  if (*status_out != CL_SUCCESS)
    FAIL("Unable to set kernel argument 2 with code: %s\n", fail);
  *status_out = clSetKernelArg(kern, 3, sizeof(float) * LOCAL_SIZE, NULL);
  if (*status_out != CL_SUCCESS)
    FAIL("Unable to set kernel argument 3 with code: %s\n", fail);
  *status_out = clEnqueueNDRangeKernel(h_cq, kern, 1, NULL, &g_size, &l_size, 0,
                                       NULL, NULL);
  if (*status_out != CL_SUCCESS)
    FAIL("Failed to execute kernel with code: %s\n", fail);
  *status_out = clEnqueueReadBuffer(h_cq, buffers[1].dev, CL_TRUE, 0,
                                    num_groups * sizeof(float), buffers[1].host,
                                    0, NULL, NULL);
  if (*status_out != CL_SUCCESS)
    FAIL("Failed to read buffer B with code: %s\n", fail);
  return h_reduce(num_groups, buffers[1].host);
fail:
  *errmsg_out = errmsg;
  return *status_out;
}
*/

cl_int cl_solve_setup(cl_bundle bundle, size_t sim_size,
                      const float *restrict h_x, const float *restrict h_x0,
                      const char **errmsg_out) {
  return cl_setup(bundle, sim_size, h_x, h_x0, NULL, errmsg_out);
}

cl_int cl_solve_retrieve(cl_bundle bundle, size_t sim_size, float *h_x,
                         const char **errmsg_out) {
  return cl_retrieve(bundle, sim_size, h_x, NULL, NULL, errmsg_out);
}

cl_int cl_solve_step(cl_bundle bundle, unsigned int sim_size, float a, float c,
                     bool negate_axes[2], const char **errmsg_out) {
  const char *errmsg = NULL;
  cl_int status;
  float c_inv = 1.0f / c;
  static const size_t l_sizes[2] = {L_SIZE, L_SIZE};
  size_t g_sizes[2] = {sim_size / P_SIZE, sim_size};
  // size_t g_sizes[2] = {sim_size, sim_size};
  status = clSetKernelArg(bundle->kernels[JACOBI_IX], 0, sizeof(unsigned int),
                          &sim_size);
  if (status != CL_SUCCESS)
    FAIL("Unable to set kernel argument 0 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[JACOBI_IX], 1, sizeof(cl_mem),
                          &bundle->d_buffers[0]);
  if (status != CL_SUCCESS)
    FAIL("Unable to set jacobi kernel argument 1 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[JACOBI_IX], 2, sizeof(cl_mem),
                          &bundle->d_buffers[1]);
  if (status != CL_SUCCESS)
    FAIL("Unable to set jacobi kernel argument 2 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[JACOBI_IX], 3, sizeof(cl_mem),
                          &bundle->d_buffers[2]);
  if (status != CL_SUCCESS)
    FAIL("Unable to set jacobi kernel argument 3 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[JACOBI_IX], 4,
                          sizeof(float) * L_ACTUAL_SIZE, NULL);
  if (status != CL_SUCCESS)
    FAIL("Unable to set jacobi kernel argument 4 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[JACOBI_IX], 5, sizeof(float), &a);
  if (status != CL_SUCCESS)
    FAIL("Unable to set jacobi kernel argument 5 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[JACOBI_IX], 6, sizeof(float), &c_inv);
  if (status != CL_SUCCESS)
    FAIL("Unable to set jacobi kernel argument 6 with code: %s\n", fail);

  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[JACOBI_IX], 2,
                                  NULL, g_sizes, l_sizes, 0, NULL, NULL);
  if (status != CL_SUCCESS)
    FAIL("Failed to execute jacobi kernel with code: %s\n", fail);
  status = clFinish(bundle->h_cq);
  if (status != CL_SUCCESS)
    FAIL("Failed to finish jacobi kernel execution with code: %s\n", fail);

  g_sizes[0] = sim_size, g_sizes[1] = 4;
  int negate_rows = negate_axes[0], negate_cols = negate_axes[1];
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 0, sizeof(unsigned int),
                          &sim_size);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd kernel argument 0 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 1, sizeof(cl_mem),
                          &bundle->d_buffers[2]);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd kernel argument 1 with code: %s\n", fail);
  status =
      clSetKernelArg(bundle->kernels[SET_BND_IX], 2, sizeof(int), &negate_rows);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd kernel argument 2 with code: %s\n", fail);
  status =
      clSetKernelArg(bundle->kernels[SET_BND_IX], 3, sizeof(int), &negate_cols);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd kernel argument 3 with code: %s\n", fail);
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[SET_BND_IX], 2,
                                  NULL, g_sizes, NULL, 0, NULL, NULL);
  if (status != CL_SUCCESS)
    FAIL("Failed to execute set_bnd kernel with code: %s\n", fail);
  status = clFinish(bundle->h_cq);
  if (status != CL_SUCCESS)
    FAIL("Failed to finish set_bnd kernel execution with code: %s\n", fail);
  // swap d_x and d_x1; most recent value will be in d_x after swap
  cl_mem temp = bundle->d_buffers[0];
  bundle->d_buffers[0] = bundle->d_buffers[2];
  bundle->d_buffers[2] = temp;
fail:
  *errmsg_out = errmsg;
  return status;
}

cl_int cl_project_setup(cl_bundle bundle, size_t sim_size,
                        const float *restrict h_u, const float *restrict h_v,
                        const char **errmsg_out) {
  // swap device buffers 0 and 2 so that `p` is preserved. this is unnecessary
  // for the first phase of projection, but the swap is cheap and the setup for
  // both phases is otherwise identical, so we always do it
  cl_mem temp = bundle->d_buffers[2];
  bundle->d_buffers[2] = bundle->d_buffers[0];
  bundle->d_buffers[0] = temp;
  return cl_setup(bundle, sim_size, h_u, h_v, NULL, errmsg_out);
}

cl_int cl_project_retrieve(cl_bundle bundle, size_t sim_size,
                           float *restrict h_u, float *restrict h_v,
                           const char **errmsg_out) {
  return cl_retrieve(bundle, sim_size, h_u, h_v, NULL, errmsg_out);
}

cl_int cl_project_one(cl_bundle bundle, unsigned int sim_size,
                      const char **errmsg_out) {
  const char *errmsg = NULL;
  cl_int status;
  static const size_t l_sizes[2] = {L_SIZE, L_SIZE};
  size_t g_sizes[2] = {sim_size / P_SIZE, sim_size};

  status = clSetKernelArg(bundle->kernels[PROJECT_ONE_IX], 0,
                          sizeof(unsigned int), &sim_size);
  if (status != CL_SUCCESS)
    FAIL("Unable to set project_one argument 0 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[PROJECT_ONE_IX], 1, sizeof(cl_mem),
                          &bundle->d_buffers[2]);
  if (status != CL_SUCCESS)
    FAIL("Unable to set project_one argument 1 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[PROJECT_ONE_IX], 2, sizeof(cl_mem),
                          &bundle->d_buffers[0]);
  if (status != CL_SUCCESS)
    FAIL("Unable to set project_one argument 2 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[PROJECT_ONE_IX], 3, sizeof(cl_mem),
                          &bundle->d_buffers[1]);
  if (status != CL_SUCCESS)
    FAIL("Unable to set project_one argument 3 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[PROJECT_ONE_IX], 4,
                          sizeof(float) * L_ACTUAL_SIZE, NULL);
  if (status != CL_SUCCESS)
    FAIL("Unable to set project_one argument 4 with code: %s\n", fail);

  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[PROJECT_ONE_IX],
                                  2, NULL, g_sizes, l_sizes, 0, NULL, NULL);
  if (status != CL_SUCCESS)
    FAIL("Failed to execute project_one with code: %s\n", fail);
  status = clFinish(bundle->h_cq);
  if (status != CL_SUCCESS)
    FAIL("Failed to finish project_one execution with code: %s\n", fail);

  // set_bnd on div
  g_sizes[0] = sim_size, g_sizes[1] = 4;
  static const int I_FALSE = 0;
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 0, sizeof(unsigned int),
                          &sim_size);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd argument 0 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 1, sizeof(cl_mem),
                          &bundle->d_buffers[2]);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd argument 1 with code: %s\n", fail);
  status =
      clSetKernelArg(bundle->kernels[SET_BND_IX], 2, sizeof(int), &I_FALSE);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd argument 2 with code: %s\n", fail);
  status =
      clSetKernelArg(bundle->kernels[SET_BND_IX], 3, sizeof(int), &I_FALSE);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd argument 3 with code: %s\n", fail);
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[SET_BND_IX], 2,
                                  NULL, g_sizes, NULL, 0, NULL, NULL);
  if (status != CL_SUCCESS)
    FAIL("Failed to execute set_bnd kernel with code: %s\n", fail);

  // zero fill p
  static const float ZERO = 0.0f;
  status = clEnqueueFillBuffer(bundle->h_cq, bundle->d_buffers[0], &ZERO,
                               sizeof(float), 0, ACTUAL_SIZE * sizeof(float), 0,
                               NULL, NULL);
  if (status != CL_SUCCESS)
    FAIL("Failed to zero-fill device buffer 0 with code: %s\n", fail);

  // swap device buffers 1 and 2 so that div is in the x0 buffer location
  cl_mem temp = bundle->d_buffers[1];
  bundle->d_buffers[1] = bundle->d_buffers[2];
  bundle->d_buffers[2] = temp;

  // wait for set_bnd on div and zero-fill of p to finish
  status = clFinish(bundle->h_cq);
  if (status != CL_SUCCESS)
    FAIL("Failed to finish set_bnd executions with code: %s\n", fail);
fail:
  *errmsg_out = errmsg;
  return status;
}

cl_int cl_project_two(cl_bundle bundle, unsigned int sim_size,
                      const char **errmsg_out) {
  const char *errmsg = NULL;
  cl_int status;
  static const size_t l_sizes[2] = {L_SIZE, L_SIZE};
  size_t g_sizes[2] = {sim_size / P_SIZE, sim_size};

  status = clSetKernelArg(bundle->kernels[PROJECT_TWO_IX], 0,
                          sizeof(unsigned int), &sim_size);
  if (status != CL_SUCCESS)
    FAIL("Unable to set project_two argument 0 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[PROJECT_TWO_IX], 1, sizeof(cl_mem),
                          &bundle->d_buffers[0]);
  if (status != CL_SUCCESS)
    FAIL("Unable to set project_two argument 1 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[PROJECT_TWO_IX], 2, sizeof(cl_mem),
                          &bundle->d_buffers[1]);
  if (status != CL_SUCCESS)
    FAIL("Unable to set project_two argument 2 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[PROJECT_TWO_IX], 3, sizeof(cl_mem),
                          &bundle->d_buffers[2]);
  if (status != CL_SUCCESS)
    FAIL("Unable to set project_two argument 3 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[PROJECT_TWO_IX], 4,
                          sizeof(float) * L_ACTUAL_SIZE, NULL);
  if (status != CL_SUCCESS)
    FAIL("Unable to set project_two argument 4 with code: %s\n", fail);

  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[PROJECT_TWO_IX],
                                  2, NULL, g_sizes, l_sizes, 0, NULL, NULL);
  if (status != CL_SUCCESS)
    FAIL("Failed to execute project_two with code: %s\n", fail);
  status = clFinish(bundle->h_cq);
  if (status != CL_SUCCESS)
    FAIL("Failed to finish project_two execution with code: %s\n", fail);

  g_sizes[0] = sim_size, g_sizes[1] = 4;
  static const int I_TRUE = 1, I_FALSE = 0;

  // set_bnd on u
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 0, sizeof(unsigned int),
                          &sim_size);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd argument 0 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 1, sizeof(cl_mem),
                          &bundle->d_buffers[0]);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd argument 1 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 2, sizeof(int), &I_TRUE);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd argument 2 with code: %s\n", fail);
  status =
      clSetKernelArg(bundle->kernels[SET_BND_IX], 3, sizeof(int), &I_FALSE);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd argument 3 with code: %s\n", fail);
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[SET_BND_IX], 2,
                                  NULL, g_sizes, NULL, 0, NULL, NULL);
  if (status != CL_SUCCESS)
    FAIL("Failed to execute set_bnd kernel with code: %s\n", fail);

  // set_bnd on v
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 0, sizeof(unsigned int),
                          &sim_size);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd argument 0 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 1, sizeof(cl_mem),
                          &bundle->d_buffers[1]);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd argument 1 with code: %s\n", fail);
  status =
      clSetKernelArg(bundle->kernels[SET_BND_IX], 2, sizeof(int), &I_FALSE);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd argument 2 with code: %s\n", fail);
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 3, sizeof(int), &I_TRUE);
  if (status != CL_SUCCESS)
    FAIL("Failed to set set_bnd argument 3 with code: %s\n", fail);
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[SET_BND_IX], 2,
                                  NULL, g_sizes, NULL, 0, NULL, NULL);
  if (status != CL_SUCCESS)
    FAIL("Failed to execute set_bnd kernel with code: %s\n", fail);

  // wait for set_bnd on u and v to finish
  status = clFinish(bundle->h_cq);
  if (status != CL_SUCCESS)
    FAIL("Failed to finish set_bnd executions with code: %s\n", fail);
fail:
  *errmsg_out = errmsg;
  return status;
}
