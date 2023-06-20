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

#define JACOBI_IX 0
#define SET_BND_IX 1
#define PROJECT_ONE_IX 2
#define PROJECT_TWO_IX 3
#define ADVECT_IX 4
#define L_SIZE 16
#define L_ACTUAL_SIZE \
  ((P_SIZE * L_SIZE + 2 * COL_BORDER) * (L_SIZE + 2 * ROW_BORDER))
#define FAIL(msg, label) \
  do {                   \
    errmsg = msg;        \
    goto label;          \
  } while (0)
#define D_SWAP(d_x, d_x0) \
  do {                    \
    cl_mem temp = d_x0;   \
    d_x0 = d_x;           \
    d_x = temp;           \
  } while (0)
#define GEN_FAIL_READ_MSG(buffer_num) \
  "Failed to read device buffer " #buffer_num " with code: %s.\n"
#define GEN_FAIL_WRITE_MSG(buffer_num) \
  "Failed to write device buffer " #buffer_num " with code: %s.\n"
#define GEN_FAIL_ARG_MSG(kernel_name, arg_num) \
  "Failed to set " #kernel_name " argument " #arg_num " with code: %s.\n"
#define GEN_FAIL_ENQUEUE_MSG(kernel_name) \
  "Failed to enqueue " #kernel_name " with code: %s.\n"
#define GEN_FAIL_CL_CREATE_MSG(object_name) \
  "Failed to create OpenCL " #object_name " with code: %s.\n"
#define GEN_FAIL_CL_RETRIEVE_MSG(object_name) \
  "Failed to retrieve OpenCL " #object_name " with code: %s.\n"
#define GEN_FAIL_ALLOC_MSG(object_name) \
  "Failed to allocate memory for " #object_name " with code: %s.\n"

static const char *OPTS_FMT = "-I%s/src -cl-std=CL2.0";
static const char *KERNEL_NAMES[] = {"jacobi", "set_bnd", "project_one",
                                     "project_two", "advect"};
static const size_t L_SIZES[] = {L_SIZE, L_SIZE};
static const char *ERRMSG_READ_PROG =
    "Failed to read OpenCL program source: buffer size exceeded.\n";
static const char *ERRMSG_GET_CWD =
    "Failed to get current working directory.\n";
static const char *ERRMSG_BUILD_PROG =
    "Failed to build OpenCL program with code: %s.\n";
static const char *ERRMSG_HANDLE_ERR =
    "An error occurred while handling a previous error: %s\n";
static const char *ERRMSG_UNLOAD =
    "Failed to unload platform compiler with code: %s.\n";
static const char *ERRMSG_GET_GPU =
    "Failed to locate GPU device with code: %s.\n";
static const char *ERRMSG_GET_CPU =
    "Failed to locate CPU device with code: %s.\n";
static const char *ERRMSG_ALLOC_BUNDLE = GEN_FAIL_ALLOC_MSG(bundle);
static const char *ERRMSG_ALLOC_OPTS = GEN_FAIL_ALLOC_MSG(program options);
static const char *ERRMSG_CREATE_PROG = GEN_FAIL_CL_CREATE_MSG(program);
static const char *ERRMSG_CREATE_BUF = GEN_FAIL_CL_CREATE_MSG(buffer);
static const char *ERRMSG_CREATE_CTX = GEN_FAIL_CL_CREATE_MSG(context);
static const char *ERRMSG_CREATE_HQ = GEN_FAIL_CL_CREATE_MSG(host queue);
static const char *ERRMSG_CREATE_DQ = GEN_FAIL_CL_CREATE_MSG(device queue);
static const char *ERRMSG_CREATE_KERN = GEN_FAIL_CL_CREATE_MSG(kernel);
static const char *ERRMSG_GET_PLAT_IDS = GEN_FAIL_CL_RETRIEVE_MSG(platform ids);
static const char *ERRMSG_GET_DEV_IDS = GEN_FAIL_CL_RETRIEVE_MSG(device ids);

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
    FAIL((status = -1, ERRMSG_READ_PROG), fail_read);
  else
    buffer[len] = '\0';

  const char *buffer_ptr = buffer;
  cl_program prog =
      clCreateProgramWithSource(bundle->ctx, 1, &buffer_ptr, NULL, &status);
  if (status != CL_SUCCESS) FAIL(ERRMSG_CREATE_PROG, fail_prog);

  char cwd[4096];
  if (getcwd(cwd, 4096) == NULL) FAIL(ERRMSG_GET_CWD, fail_getcwd);
  size_t opts_len = snprintf(NULL, 0, OPTS_FMT, cwd);
  char *opts = malloc(sizeof(char) * (opts_len + 1));
  if (opts == NULL) FAIL(ERRMSG_ALLOC_OPTS, fail_alloc);
  (void)snprintf(opts, opts_len + 1, OPTS_FMT, cwd);

  status = clBuildProgram(prog, 0, NULL, opts, NULL, NULL);
  if (status != CL_SUCCESS) {
    clGetProgramBuildInfo(prog, bundle->dev, CL_PROGRAM_BUILD_LOG,
                          sizeof(buffer), buffer, &len);
    fprintf(stderr, "%s\n", buffer);
    FAIL(ERRMSG_BUILD_PROG, fail_build);
  }

  bundle->kernels[kern_ix] =
      clCreateKernel(prog, KERNEL_NAMES[kern_ix], &status);
  if (status != CL_SUCCESS) FAIL(ERRMSG_CREATE_KERN, fail_kernel);
fail_kernel:
fail_build:
  free(opts);
fail_alloc:
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
    if (*status != CL_SUCCESS) FAIL(ERRMSG_CREATE_BUF, fail);
  }
  *errmsg_out = errmsg;
  return 0;
fail:
  *errmsg_out = errmsg;
  size_t j = 0;
  for (; j < i; ++j) {
    cl_int fail_status = clReleaseMemObject(bundle->d_buffers[j]);
    if (fail_status != CL_SUCCESS)
      fprintf(stderr, ERRMSG_HANDLE_ERR, cl_strerror(fail_status));
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
    FAIL((*status = CL_SUCCESS, ERRMSG_ALLOC_BUNDLE), fail_alloc);
  ret->d_cq = NULL;
  cl_platform_id platforms[4];
  unsigned int num_platforms;
  *status = clGetPlatformIDs(4, platforms, &num_platforms);
  if (*status != CL_SUCCESS) FAIL(ERRMSG_GET_PLAT_IDS, fail_plat);

  for (size_t i = 0; i < 4 && i < num_platforms; ++i) {
    *status =
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_GPU, 1, &ret->dev, NULL);
    if (*status == CL_DEVICE_NOT_FOUND) {
      continue;
    } else if (*status != CL_SUCCESS) {
      FAIL(ERRMSG_GET_DEV_IDS, fail_dev);
    } else {
      ret->plat = platforms[i];
      break;
    }
  }

  if (ret->dev == NULL) FAIL(ERRMSG_GET_GPU, fail_dev);

  ret->ctx = clCreateContext(NULL, 1, &ret->dev, NULL, NULL, status);
  if (*status != CL_SUCCESS) FAIL(ERRMSG_CREATE_CTX, fail_ctx);

  ret->h_cq =
      clCreateCommandQueueWithProperties(ret->ctx, ret->dev, NULL, status);
  if (*status != CL_SUCCESS) FAIL(ERRMSG_CREATE_HQ, fail_h_cq);

  ret->d_cq = clCreateCommandQueueWithProperties(
      ret->ctx, ret->dev,
      (cl_queue_properties[]){CL_QUEUE_PROPERTIES,
                              CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE |
                                  CL_QUEUE_ON_DEVICE |
                                  CL_QUEUE_ON_DEVICE_DEFAULT,
                              0},
      status);
  if (*status != CL_SUCCESS) FAIL(ERRMSG_CREATE_DQ, fail_d_cq);

  *status = init_kernel(ret, JACOBI_IX, &errmsg);
  if (*status != CL_SUCCESS) FAIL(errmsg, fail_jacobi);
  *status = init_kernel(ret, SET_BND_IX, errmsg_out);
  if (*status != CL_SUCCESS) FAIL(errmsg, fail_set_bnd);
  *status = init_kernel(ret, PROJECT_ONE_IX, errmsg_out);
  if (*status != CL_SUCCESS) FAIL(errmsg, fail_project_one);
  *status = init_kernel(ret, PROJECT_TWO_IX, errmsg_out);
  if (*status != CL_SUCCESS) FAIL(errmsg, fail_project_two);
  *status = init_kernel(ret, ADVECT_IX, errmsg_out);
  if (*status != CL_SUCCESS) FAIL(errmsg, fail_advect);
  *status = init_buffers(ret, sim_size, status, errmsg_out);
  if (*status != CL_SUCCESS) FAIL(errmsg, fail_buffers);
  *errmsg_out = errmsg;
  return ret;
  cl_int unload_status;
fail_buffers:
  clReleaseKernel(ret->kernels[ADVECT_IX]);
fail_advect:
  clReleaseKernel(ret->kernels[PROJECT_TWO_IX]);
fail_project_two:
  clReleaseKernel(ret->kernels[PROJECT_ONE_IX]);
fail_project_one:
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
    fprintf(stderr, ERRMSG_UNLOAD, cl_strerror(unload_status));
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
    FAIL((*status = CL_SUCCESS, ERRMSG_ALLOC_BUNDLE), fail_alloc);
  ret->d_cq = NULL;
  cl_platform_id platforms[4];
  unsigned int num_platforms;
  *status = clGetPlatformIDs(4, platforms, &num_platforms);
  if (*status != CL_SUCCESS) FAIL(ERRMSG_GET_PLAT_IDS, fail_plat);

  for (size_t i = 0; i < 4 && i < num_platforms; ++i) {
    *status =
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_CPU, 1, &ret->dev, NULL);
    if (*status == CL_DEVICE_NOT_FOUND) {
      continue;
    } else if (*status != CL_SUCCESS) {
      FAIL(ERRMSG_GET_DEV_IDS, fail_dev);
    } else {
      ret->plat = platforms[i];
      break;
    }
  }

  if (ret->dev == NULL) FAIL(ERRMSG_GET_CPU, fail_dev);

  ret->ctx = clCreateContext(NULL, 1, &ret->dev, NULL, NULL, status);
  if (*status != CL_SUCCESS) FAIL(ERRMSG_CREATE_CTX, fail_ctx);

  ret->h_cq =
      clCreateCommandQueueWithProperties(ret->ctx, ret->dev, NULL, status);
  if (*status != CL_SUCCESS) FAIL(ERRMSG_CREATE_HQ, fail_h_cq);

  *errmsg_out = errmsg;
  return ret;
  cl_int unload_status;
fail_h_cq:
  clReleaseContext(ret->ctx);
fail_ctx:
  unload_status = clUnloadPlatformCompiler(ret->plat);
  if (unload_status != CL_SUCCESS)
    fprintf(stderr, "Failed to unload platform compiler with code: %s.\n",
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

static cl_int configure_jacobi(cl_bundle bundle, unsigned int sim_size,
                               const cl_mem d_x0, float a, float c,
                               const char **errmsg_out) {
  const char *errmsg = NULL;
  cl_int status = clSetKernelArg(bundle->kernels[JACOBI_IX], 0,
                                 sizeof(unsigned int), &sim_size);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(jacobi, 0), fail);
  status = clSetKernelArg(bundle->kernels[JACOBI_IX], 2, sizeof(cl_mem), &d_x0);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(jacobi, 2), fail);
  status = clSetKernelArg(bundle->kernels[JACOBI_IX], 4,
                          sizeof(float) * L_ACTUAL_SIZE, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(jacobi, 4), fail);
  status = clSetKernelArg(bundle->kernels[JACOBI_IX], 5, sizeof(float), &a);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(jacobi, 5), fail);

  const float c_inv = 1.0f / c;
  status = clSetKernelArg(bundle->kernels[JACOBI_IX], 6, sizeof(float), &c_inv);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(jacobi, 6), fail);
fail:
  *errmsg_out = errmsg;
  return status;
}

static cl_int configure_set_bnd(cl_bundle bundle, unsigned int sim_size,
                                const int bnd_opts[3],
                                const char **errmsg_out) {
  const char *errmsg = NULL;
  cl_int status = clSetKernelArg(bundle->kernels[SET_BND_IX], 0,
                                 sizeof(unsigned int), &sim_size);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(set_bnd, 0), fail);
  status =
      clSetKernelArg(bundle->kernels[SET_BND_IX], 2, sizeof(int), &bnd_opts[0]);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(set_bnd, 2), fail);
  status =
      clSetKernelArg(bundle->kernels[SET_BND_IX], 3, sizeof(int), &bnd_opts[1]);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(set_bnd, 3), fail);
  status =
      clSetKernelArg(bundle->kernels[SET_BND_IX], 4, sizeof(int), &bnd_opts[2]);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(set_bnd, 4), fail);
fail:
  *errmsg_out = errmsg;
  return status;
}

static cl_int configure_advect(cl_bundle bundle, unsigned int sim_size,
                               const cl_mem d_x, const cl_mem d_x0,
                               const cl_mem d_u, const cl_mem d_v, float dt,
                               const char **errmsg_out) {
  const char *errmsg = NULL;
  cl_int status = clSetKernelArg(bundle->kernels[ADVECT_IX], 0,
                                 sizeof(unsigned int), &sim_size);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(advect, 0), fail);
  status = clSetKernelArg(bundle->kernels[ADVECT_IX], 1, sizeof(cl_mem), &d_x);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(advect, 1), fail);
  status = clSetKernelArg(bundle->kernels[ADVECT_IX], 2, sizeof(cl_mem), &d_x0);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(advect, 2), fail);
  status = clSetKernelArg(bundle->kernels[ADVECT_IX], 3, sizeof(cl_mem), &d_u);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(advect, 3), fail);
  status = clSetKernelArg(bundle->kernels[ADVECT_IX], 4, sizeof(cl_mem), &d_v);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(advect, 4), fail);
  status = clSetKernelArg(bundle->kernels[ADVECT_IX], 5, sizeof(float), &dt);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(advect, 5), fail);
fail:
  *errmsg_out = errmsg;
  return status;
}

static cl_int enqueue_solve(cl_bundle bundle, unsigned int sim_size, cl_mem d_x,
                            cl_mem d_x1, bool set_corners, size_t iterations,
                            const char **errmsg_out) {
  const char *errmsg = NULL;
  cl_int status = CL_SUCCESS;
  const size_t jac_sizes[2] = {sim_size / P_SIZE, sim_size};
  const size_t bnd_sizes[2] = {sim_size, 4};
  for (size_t i = 0; i < iterations; ++i) {
    // set per-iteration jacobi arguments
    status =
        clSetKernelArg(bundle->kernels[JACOBI_IX], 1, sizeof(cl_mem), &d_x);
    if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(jacobi, 1), fail);
    status =
        clSetKernelArg(bundle->kernels[JACOBI_IX], 3, sizeof(cl_mem), &d_x1);
    if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(jacobi, 3), fail);

    status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[JACOBI_IX], 2,
                                    NULL, jac_sizes, L_SIZES, 0, NULL, NULL);
    if (status != CL_SUCCESS) FAIL(GEN_FAIL_ENQUEUE_MSG(jacobi), fail);

    // set per-iteration set_bnd arguments
    status =
        clSetKernelArg(bundle->kernels[SET_BND_IX], 1, sizeof(cl_mem), &d_x1);
    if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(set_bnd, 1), fail);
    // set corners on final iteration
    if (i == iterations - 1 && set_corners) {
      static const int I_TRUE = true;
      status =
          clSetKernelArg(bundle->kernels[SET_BND_IX], 4, sizeof(int), &I_TRUE);
      if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(set_bnd, 4), fail);
    }

    status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[SET_BND_IX],
                                    2, NULL, bnd_sizes, NULL, 0, NULL, NULL);
    if (status != CL_SUCCESS) FAIL(GEN_FAIL_ENQUEUE_MSG(set_bnd), fail);

    D_SWAP(d_x, d_x1);
  }
fail:
  *errmsg_out = errmsg;
  return status;
}

cl_int cl_dens_step_full(cl_bundle bundle, const size_t sim_size,
                         float *restrict h_x, const float *restrict h_x0,
                         const float *restrict h_u, const float *restrict h_v,
                         const float diff, const float dt,
                         const size_t iterations, const char **errmsg_out) {
  const char *errmsg = NULL;
  cl_int status;

  // write host arrays to device for solver
  status =
      clEnqueueWriteBuffer(bundle->h_cq, bundle->d_buffers[0], CL_FALSE, 0,
                           sizeof(float) * ACTUAL_SIZE, h_x, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_WRITE_MSG(0), fail);
  status =
      clEnqueueWriteBuffer(bundle->h_cq, bundle->d_buffers[1], CL_FALSE, 0,
                           sizeof(float) * ACTUAL_SIZE, h_x0, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_WRITE_MSG(1), fail);

  // configure and enqueue all kernels for the solver
  float a = dt * diff * sim_size * sim_size, c = 1 + 4 * a;
  status =
      configure_jacobi(bundle, sim_size, bundle->d_buffers[1], a, c, &errmsg);
  if (status != CL_SUCCESS) FAIL(errmsg, fail);

  status = configure_set_bnd(bundle, sim_size, (int[3]){false, false, false},
                             &errmsg);
  if (status != CL_SUCCESS) FAIL(errmsg, fail);
  status = enqueue_solve(bundle, sim_size, bundle->d_buffers[0],
                         bundle->d_buffers[1], true, iterations, &errmsg);
  if (status) FAIL(errmsg, fail);

  // swap d_x and d_x0
  D_SWAP(bundle->d_buffers[0], bundle->d_buffers[1]);

  // write host arrays to device for advection
  status =
      clEnqueueWriteBuffer(bundle->h_cq, bundle->d_buffers[2], CL_FALSE, 0,
                           sizeof(float) * ACTUAL_SIZE, h_u, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_WRITE_MSG(2), fail);
  status =
      clEnqueueWriteBuffer(bundle->h_cq, bundle->d_buffers[3], CL_FALSE, 0,
                           sizeof(float) * ACTUAL_SIZE, h_v, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_WRITE_MSG(3), fail);

  status = configure_advect(bundle, sim_size, bundle->d_buffers[0],
                            bundle->d_buffers[1], bundle->d_buffers[2],
                            bundle->d_buffers[3], dt, &errmsg);
  if (status != CL_SUCCESS) FAIL(errmsg, fail);
  const size_t adv_sizes[2] = {sim_size / P_SIZE, sim_size};
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[ADVECT_IX], 2,
                                  NULL, adv_sizes, L_SIZES, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ENQUEUE_MSG(advect), fail);

  // set arguments for final set_bnd enqueue. `sim_size`, `negate_rows` and
  // `negate_cols` should already have been set by the last enqueue
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 1, sizeof(cl_mem),
                          &bundle->d_buffers[0]);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(set_bnd, 1), fail);
  const int I_FALSE = false;
  status =
      clSetKernelArg(bundle->kernels[SET_BND_IX], 4, sizeof(int), &I_FALSE);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(set_bnd, 4), fail);

  const size_t bnd_sizes[2] = {sim_size, 4};
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[SET_BND_IX], 2,
                                  NULL, bnd_sizes, NULL, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ENQUEUE_MSG(set_bnd), fail);

  status = clEnqueueReadBuffer(bundle->h_cq, bundle->d_buffers[0], CL_TRUE, 0,
                               sizeof(float) * ACTUAL_SIZE, h_x, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_READ_MSG(0), fail);
fail:
  *errmsg_out = errmsg;
  return status;
}

static cl_int configure_project_one(cl_bundle bundle, unsigned int sim_size,
                                    const cl_mem d_div, const cl_mem d_u,
                                    const cl_mem d_v, const char **errmsg_out) {
  const char *errmsg = NULL;

  cl_int status = clSetKernelArg(bundle->kernels[PROJECT_ONE_IX], 0,
                                 sizeof(unsigned int), &sim_size);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(project_one, 0), fail);
  status = clSetKernelArg(bundle->kernels[PROJECT_ONE_IX], 1, sizeof(cl_mem),
                          &d_div);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(project_one, 1), fail);
  status =
      clSetKernelArg(bundle->kernels[PROJECT_ONE_IX], 2, sizeof(cl_mem), &d_u);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(project_one, 2), fail);
  status =
      clSetKernelArg(bundle->kernels[PROJECT_ONE_IX], 3, sizeof(cl_mem), &d_v);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(project_one, 3), fail);
  status = clSetKernelArg(bundle->kernels[PROJECT_ONE_IX], 4,
                          sizeof(float) * L_ACTUAL_SIZE, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(project_one, 4), fail);

fail:
  *errmsg_out = errmsg;
  return status;
}

static cl_int configure_project_two(cl_bundle bundle, unsigned int sim_size,
                                    const cl_mem d_u, const cl_mem d_v,
                                    const cl_mem d_p, const char **errmsg_out) {
  const char *errmsg = NULL;

  cl_int status = clSetKernelArg(bundle->kernels[PROJECT_TWO_IX], 0,
                                 sizeof(unsigned int), &sim_size);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(project_two, 0), fail);
  status =
      clSetKernelArg(bundle->kernels[PROJECT_TWO_IX], 1, sizeof(cl_mem), &d_u);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(project_two, 1), fail);
  status =
      clSetKernelArg(bundle->kernels[PROJECT_TWO_IX], 2, sizeof(cl_mem), &d_v);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(project_two, 2), fail);
  status =
      clSetKernelArg(bundle->kernels[PROJECT_TWO_IX], 3, sizeof(cl_mem), &d_p);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(project_two, 3), fail);
  status = clSetKernelArg(bundle->kernels[PROJECT_TWO_IX], 4,
                          sizeof(float) * L_ACTUAL_SIZE, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ARG_MSG(project_two, 4), fail);

fail:
  *errmsg_out = errmsg;
  return status;
}

static cl_int cl_project(cl_bundle bundle, unsigned int sim_size,
                         const cl_mem d_u, const cl_mem d_v, const cl_mem d_p,
                         const cl_mem d_div, const cl_mem d_scratch,
                         size_t iterations, const char **errmsg_out) {
  const char *errmsg = NULL;

  cl_int status =
      configure_project_one(bundle, sim_size, d_div, d_u, d_v, &errmsg);
  if (status) FAIL(errmsg, fail);
  const size_t proj_sizes[2] = {sim_size / P_SIZE, sim_size};
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[PROJECT_ONE_IX],
                                  2, NULL, proj_sizes, L_SIZES, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ENQUEUE_MSG(project_one), fail);

  status = configure_set_bnd(bundle, sim_size, (int[3]){false, false, false},
                             &errmsg);
  if (status) FAIL(errmsg, fail);
  status =
      clSetKernelArg(bundle->kernels[SET_BND_IX], 1, sizeof(cl_mem), &d_div);
  if (status) FAIL(GEN_FAIL_ARG_MSG(set_bnd, 1), fail);
  const size_t bnd_sizes[] = {sim_size, 4};
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[SET_BND_IX], 2,
                                  NULL, bnd_sizes, NULL, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ENQUEUE_MSG(set_bnd), fail);

  static const float ZERO = 0.0f;
  status = clEnqueueFillBuffer(bundle->h_cq, d_p, &ZERO, sizeof(float), 0,
                               ACTUAL_SIZE * sizeof(float), 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL("Failed to zero-fill p with code: %s\n", fail);

  status = configure_jacobi(bundle, sim_size, d_div, 1, 4, &errmsg);
  if (status) FAIL(errmsg, fail);
  status = enqueue_solve(bundle, sim_size, d_p, d_scratch, false, iterations,
                         &errmsg);
  if (status) FAIL(errmsg, fail);

  status = configure_project_two(bundle, sim_size, d_u, d_v, d_p, &errmsg);
  if (status) FAIL(errmsg, fail);
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[PROJECT_TWO_IX],
                                  2, NULL, proj_sizes, L_SIZES, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ENQUEUE_MSG(project_two), fail);

  status = configure_set_bnd(bundle, sim_size, (int[3]){false, true, false},
                             &errmsg);
  if (status) FAIL(errmsg, fail);
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 1, sizeof(cl_mem), &d_u);
  if (status) FAIL(GEN_FAIL_ARG_MSG(set_bnd, 1), fail);
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[SET_BND_IX], 2,
                                  NULL, bnd_sizes, NULL, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ENQUEUE_MSG(set_bnd), fail);

  status = configure_set_bnd(bundle, sim_size, (int[3]){true, false, false},
                             &errmsg);
  if (status) FAIL(errmsg, fail);
  status = clSetKernelArg(bundle->kernels[SET_BND_IX], 1, sizeof(cl_mem), &d_v);
  if (status) FAIL(GEN_FAIL_ARG_MSG(set_bnd, 1), fail);
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[SET_BND_IX], 2,
                                  NULL, bnd_sizes, NULL, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ENQUEUE_MSG(set_bnd), fail);

fail:
  *errmsg_out = errmsg;
  return status;
}

cl_int cl_vel_step_full(cl_bundle bundle, const size_t sim_size,
                        float *restrict h_u, float *restrict h_v,
                        const float *restrict h_u0, const float *restrict h_v0,
                        const float visc, const float dt,
                        const size_t iterations, const char **errmsg_out) {
  const char *errmsg = NULL;
  cl_int status;

  // write host arrays to device for diffusing u
  status =
      clEnqueueWriteBuffer(bundle->h_cq, bundle->d_buffers[0], CL_FALSE, 0,
                           sizeof(float) * ACTUAL_SIZE, h_u, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_WRITE_MSG(0), fail);
  status =
      clEnqueueWriteBuffer(bundle->h_cq, bundle->d_buffers[1], CL_FALSE, 0,
                           sizeof(float) * ACTUAL_SIZE, h_u0, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_WRITE_MSG(1), fail);

  // configure and enqueue all kernels for diffusing u
  float a = dt * visc * sim_size * sim_size, c = 1 + 4 * a;
  status =
      configure_jacobi(bundle, sim_size, bundle->d_buffers[1], a, c, &errmsg);
  if (status != CL_SUCCESS) FAIL(errmsg, fail);
  status =
      configure_set_bnd(bundle, sim_size, (int[3]){false, true, true}, &errmsg);
  if (status != CL_SUCCESS) FAIL(errmsg, fail);
  status = enqueue_solve(bundle, sim_size, bundle->d_buffers[0],
                         bundle->d_buffers[1], true, iterations, &errmsg);
  if (status) FAIL(errmsg, fail);

  // save u
  D_SWAP(bundle->d_buffers[0], bundle->d_buffers[3]);

  // write host arrays to device for diffusing v
  status =
      clEnqueueWriteBuffer(bundle->h_cq, bundle->d_buffers[0], CL_FALSE, 0,
                           sizeof(float) * ACTUAL_SIZE, h_v, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_WRITE_MSG(0), fail);
  status =
      clEnqueueWriteBuffer(bundle->h_cq, bundle->d_buffers[1], CL_FALSE, 0,
                           sizeof(float) * ACTUAL_SIZE, h_v0, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_WRITE_MSG(1), fail);

  // configure and enqueue all kernels for diffusing v
  // NOTE: this is necessary since we swapped the buffers. if we wanted to, we
  // could configure the kernel to use buffers 2, 3 and 4 instead of swapping
  // them on the host
  status =
      configure_jacobi(bundle, sim_size, bundle->d_buffers[1], a, c, &errmsg);
  if (status != CL_SUCCESS) FAIL(errmsg, fail);
  status =
      configure_set_bnd(bundle, sim_size, (int[3]){true, false, true}, &errmsg);
  if (status != CL_SUCCESS) FAIL(errmsg, fail);
  status = enqueue_solve(bundle, sim_size, bundle->d_buffers[0],
                         bundle->d_buffers[1], true, iterations, &errmsg);
  if (status) FAIL(errmsg, fail);

  // save v
  D_SWAP(bundle->d_buffers[0], bundle->d_buffers[4]);

  // at this point, the contents of the device buffers are as follows:
  // buffer 0: garbage
  // buffer 1: v0 (can be discarded)
  // buffer 2: garbage (previous iteration of v)
  // buffer 3: u
  // buffer 4: v

  status =
      cl_project(bundle, sim_size, bundle->d_buffers[3], bundle->d_buffers[4],
                 bundle->d_buffers[0], bundle->d_buffers[1],
                 bundle->d_buffers[2], iterations, &errmsg);
  if (status) FAIL(errmsg, fail);

  // at this point, the contents of the device buffers are as follows:
  // buffer 0: p
  // buffer 1: div
  // buffer 2: garbage (previous iteration of p)
  // buffer 3: u
  // buffer 4: v

  // configure and enqueue advection for p
  status = configure_advect(bundle, sim_size, bundle->d_buffers[0],
                            bundle->d_buffers[3], bundle->d_buffers[3],
                            bundle->d_buffers[4], dt, &errmsg);
  if (status) FAIL(errmsg, fail);
  const size_t adv_sizes[2] = {sim_size / P_SIZE, sim_size};
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[ADVECT_IX], 2,
                                  NULL, adv_sizes, L_SIZES, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ENQUEUE_MSG(advect), fail);

  // configure and enqueue advection for div
  status = configure_advect(bundle, sim_size, bundle->d_buffers[1],
                            bundle->d_buffers[4], bundle->d_buffers[3],
                            bundle->d_buffers[4], dt, &errmsg);
  if (status) FAIL(errmsg, fail);
  status = clEnqueueNDRangeKernel(bundle->h_cq, bundle->kernels[ADVECT_IX], 2,
                                  NULL, adv_sizes, L_SIZES, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_ENQUEUE_MSG(advect), fail);

  status =
      cl_project(bundle, sim_size, bundle->d_buffers[0], bundle->d_buffers[1],
                 bundle->d_buffers[3], bundle->d_buffers[4],
                 bundle->d_buffers[2], iterations, &errmsg);
  if (status) FAIL(errmsg, fail);

  // at this point, the contents of the device buffers are as follows:
  // buffer 0: p (modified)
  // buffer 1: div (modified)
  // buffer 2: garbage (previous iteration of p)
  // buffer 3: p (new)
  // buffer 4: div (new)

  // read contents of buffer 0 (u) and buffer 1 (v) back to the host
  status = clEnqueueReadBuffer(bundle->h_cq, bundle->d_buffers[0], CL_TRUE, 0,
                               sizeof(float) * ACTUAL_SIZE, h_u, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_READ_MSG(0), fail);
  status = clEnqueueReadBuffer(bundle->h_cq, bundle->d_buffers[1], CL_TRUE, 0,
                               sizeof(float) * ACTUAL_SIZE, h_v, 0, NULL, NULL);
  if (status != CL_SUCCESS) FAIL(GEN_FAIL_READ_MSG(1), fail);

fail:
  *errmsg_out = errmsg;
  return status;
}
