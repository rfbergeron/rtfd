# List of things to do
1. remove icky macros
2. write serial jacobi solver function
3. write sse3 jacobi solver function
4. write sse3 advect, project, set_bnd and add_source functions
5. maybe go back and write opencl or openmp versions of the above functions
   as well
6. consider aligning fluid matrix to 64 bytes when using simd intstructions.
   Since the length/width of the field will be a multiple of 4, it would be
   natural to divide the field into 4x4 units and work with those. If we 
   align to 64 byte boundaries, 

   Actually, no. With the way the matrix is currently stored in memory, the
   64 bytes in the cache line would all come from the same row, which would
   not be especially useful to us. It may improve the performance of
   horizontal additions done by the jacobi solver, but vertical additions
   would suffer.

   The memory storage method could be changed so that each 'row' of 16 floats
   in the matrix represent a 4x4 matrix embedded in the bigger matrix. This
   would help somewhat if we were to operate on the big matrix in 4x4 chunks,
   instead of row-by-row or column-by-column.
7. for now it might be best to require N to be a multiple of 4 so that we can
   have 64-byte alignment, and then allocate a 4-element wide border around the
   actual matrix; this way the same scheme used in the stock solver will work
   for the simd solver

   In the original solver, each row would 'actually' be N+2 cells wide, with
   the cells at indices 0 and N+1 serving as borders. So, for example, if N
   were 64, the borders would be at index 0 and 65, and the 'real' cells would
   occupy 1 through 64.

   In the case that we have 4 'border' cells on each side for alignment
   purposes, we would have a total of 64 + 8 = 72 cells, with the border cells
   occupying indices 0 through 3 and 68 through 71, with the 'real' cells
   being at indices 4 through 67.
