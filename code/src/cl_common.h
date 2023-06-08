#ifndef CL_COMMON_H
#define CL_COMMON_H
#define ROW_BORDER 1
#define ROW_BEGIN ROW_BORDER
#define ROW_END (sim_size + ROW_BORDER)
#define COL_BORDER 16
#define COL_BEGIN COL_BORDER
#define COL_END (sim_size + COL_BORDER)
#define ROW_SIZE (sim_size + 2 * COL_BORDER)
#define IX(i, j) (ROW_SIZE * (i) + (j))
#define ACTUAL_SIZE (ROW_SIZE * (sim_size + 2 * ROW_BORDER))
#define P_SIZE 16
#endif
