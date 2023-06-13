#include "cl_common.h"
#define TOP_BORDER 0
#define BOTTOM_BORDER 1
#define LEFT_BORDER 2
#define RIGHT_BORDER 3
#define ACTUAL_G0 \
  (get_group_id(0) * P_SIZE * get_local_size(0) + get_local_id(0))
#define G_COL_END \
  (COL_BEGIN + (get_group_id(0) + 1) * P_SIZE * get_local_size(0))
#define L_ROW_END (get_local_size(1) + ROW_BEGIN)
#define L_COL_END (get_local_size(0) * P_SIZE + COL_BEGIN)
#define L_ROW_SIZE (get_local_size(0) * P_SIZE + 2 * COL_BORDER)
#define L_IX(i, j) ((L_ROW_SIZE) * (i) + (j))

// P_SIZE defines the number of elements a single work-item is responsible for
// copying to local memory and comupting the result for. These elements will be
// located in the same row, and will be get_local_size(0) elements apart (0,
// get_local_size(0), 2 * get_local_size(0) ... (P_SIZE - 1) *
// get_local_size(0))
//
// P_SIZE and work group size have been chosen to make writing the jacobi kernel
// more convenient. currently, the kernel should be invoked as follows:
// P_SIZE = 16
// get_local_size(0) = get_local_size(1) = 16
// get_global_size(0) = sim_size / P_SIZE
// get_global_size(1) = sim_size
//
// for example, a kernel invoked with these parameters and a sim size of 256
// would have global dimensions (16 , 256), divided between 16 work groups.
//
// it is expected that all global and local arrays have 16 border columns on
// both sides, and a single border row on the top and bottom.
//
// the `x` array is copied to local memory to reduce the number of global memory
// accesses. the border elements for a region handled by a work group are also
// copied to local memory. corners are omitted, as they are not needed. the
// borders include 16 elements to the left and right of the local region, which
// should prevent work item divergence at the expense of wasting memory.
//
// the elements of `x0` and `x1` are not stored in local memory; elements in
// these arrays are only read/written once by a single work item, so there would
// be no performance benefit to doing this.
//
// NOTE: work items with consecutive ids in dimension zero are considered to be
// adjacent; therefore work item functions used to index arrays should use
// dimension zero in the second parameter and dimension one in the first
//
// NOTE: get_global_id(dim) is equivalent to get_local_size(dim)
// * get_group_id(dim) + get_local_id(dim). for our purposes, this meansss that
// get_local_size(dim) * PRIVATE_SIZE(dim) * get_group_id(dim) will give us the
// index of the first index of the work group in that dimension.
//
// NOTE: get_global_id(0) cannot be used as-is when individual work items are
// responsible for mutiple elements; the "true" id must be calculated using the
// group id, local size, and local id.

kernel void set_bnd(const unsigned int sim_size, global float *A,
                    const int negate_rows, const int negate_cols) {
  size_t offset = get_global_id(0);
  size_t border = get_global_id(1);
  switch (border) {
    case TOP_BORDER:
      A[IX(ROW_BEGIN - 1, COL_BEGIN + offset)] =
          negate_cols ? -A[IX(ROW_BEGIN, COL_BEGIN + offset)]
                      : A[IX(ROW_BEGIN, COL_BEGIN + offset)];
      break;
    case BOTTOM_BORDER:
      A[IX(ROW_END, COL_BEGIN + offset)] =
          negate_cols ? -A[IX(ROW_END - 1, COL_BEGIN + offset)]
                      : A[IX(ROW_END - 1, COL_BEGIN + offset)];
      break;
    case LEFT_BORDER:
      A[IX(ROW_BEGIN + offset, COL_BEGIN - 1)] =
          negate_rows ? -A[IX(ROW_BEGIN + offset, COL_BEGIN)]
                      : A[IX(ROW_BEGIN + offset, COL_BEGIN)];
      break;
    case RIGHT_BORDER:
      A[IX(ROW_BEGIN + offset, COL_END)] =
          negate_rows ? -A[IX(ROW_BEGIN + offset, COL_END - 1)]
                      : A[IX(ROW_BEGIN + offset, COL_END - 1)];
      break;
  }
}

// jacobi with no local memory and no assumptions about work group size
// NOTE: global work size will need to be adjusted when using this version
// of the jacobi kernel; when a work item computes multiple values it must
// know the number of work-items in the work group in order to figure out
// the stride that separates consecutive values
kernel void jacobi_simple(const unsigned int sim_size, global float *A,
                          global const float *A0, global float *A1,
                          local float *l_A, const float a, const float c_inv) {
  (void)l_A;  // unused parameter
  size_t A_begin =
      IX(get_global_id(1) + ROW_BORDER, get_global_id(0) + COL_BORDER);
  A1[A_begin] =
      (A0[A_begin] + a * (A[A_begin + 1] + A[A_begin - 1] +
                          A[A_begin + ROW_SIZE] + A[A_begin - ROW_SIZE])) *
      c_inv;
}

// jacobi with no local memory, multiple elements per work item
kernel void jacobi_multi(const unsigned int sim_size, global float *A,
                         global const float *A0, global float *A1,
                         local float *l_A, const float a, const float c_inv) {
  (void)l_A;  // unused parameter
  const size_t g_i = get_global_id(1) + ROW_BORDER;
  for (size_t g_j = ACTUAL_G0 + COL_BORDER; g_j < G_COL_END;
       g_j += get_local_size(0))
    A1[IX(g_i, g_j)] =
        (A0[IX(g_i, g_j)] + a * (A[IX(g_i + 1, g_j)] + A[IX(g_i - 1, g_j)] +
                                 A[IX(g_i, g_j + 1)] + A[IX(g_i, g_j - 1)])) *
        c_inv;
}

void local_copy(const unsigned int sim_size, global const float *A,
                local float *l_A) {
  // include left and right border elements when copying to local memory
  size_t g_i = get_global_id(1) + ROW_BORDER,
         l_i = get_local_id(1) + ROW_BORDER;
  for (size_t l_j = get_local_id(0), g_j = ACTUAL_G0; l_j < L_ROW_SIZE;
       l_j += get_local_size(0), g_j += get_local_size(0))
    l_A[L_IX(l_i, l_j)] = A[IX(g_i, g_j)];

  // TODO: replace magic number so that this function still works when private
  // and local size changes
  // (i, j) -> linear index in the range [0, 256)
  size_t col_offset = 16 * get_local_id(1) + get_local_id(0);
  // column of work item with local id (0, 0) + col_offset + COL_BEGIN
  size_t g_col =
      COL_BEGIN + get_group_id(0) * get_local_size(0) * P_SIZE + col_offset;
  size_t g_top = ROW_BEGIN + get_group_id(1) * get_local_size(1) - 1;
  size_t g_bot = ROW_BEGIN + (get_group_id(1) + 1) * get_local_size(1);

  // set top border
  l_A[L_IX(ROW_BEGIN - 1, COL_BEGIN + col_offset)] = A[IX(g_top, g_col)];
  // set bottom border
  l_A[L_IX(L_ROW_END, COL_BEGIN + col_offset)] = A[IX(g_bot, g_col)];

  work_group_barrier(CLK_LOCAL_MEM_FENCE);
}

kernel void jacobi(const unsigned int sim_size, global float *A,
                   global const float *A0, global float *A1, local float *l_A,
                   const float a, const float c_inv) {
  local_copy(sim_size, A, l_A);
  const size_t g_i = get_global_id(1) + ROW_BORDER,
               l_i = get_local_id(1) + ROW_BORDER;
  for (size_t l_j = get_local_id(0) + COL_BORDER, g_j = ACTUAL_G0 + COL_BORDER;
       l_j < L_COL_END; l_j += get_local_size(0), g_j += get_local_size(0))
    A1[IX(g_i, g_j)] =
        (A0[IX(g_i, g_j)] +
         a * (l_A[L_IX(l_i, l_j + 1)] + l_A[L_IX(l_i, l_j - 1)] +
              l_A[L_IX(l_i + 1, l_j)] + l_A[L_IX(l_i - 1, l_j)])) *
        c_inv;
}
