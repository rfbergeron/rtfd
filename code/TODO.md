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
8. the way that the vectors representing the movement of the fluid are currently
   stored prevents advection from benefitting from SIMD instructions.

   if this were to be changed, however, it would break the existing linear
   solver, and a new one would have to be written specifically for diffusing
   vectors

   the new solver would diffuse both the u and v components of the vector
   at the same time, which would actually work out pretty nicely. row-by-row
   operations would still be simple adds, while column-by-column operations
   would just require a movehl/movelh followed by an add

   the new vector scheme that would work more readily with SIMD would store
   the u and v components of vectors together, instead of in separate matrices

   the u and v components would also be compacted together, meaning that each
   16-byte block of the matrix would contain two u,v pairs

   the for loops in the advect function would store their counting variables in
   a vector to make it easier to do vector math with them

   the long one-liner arithmetic operation that occurs at the end of the for
   loops would still have to be done in serial, but

   for projection, shuffle uv[i+1,j] and uv[i,j+1] together; then shuffle
   uv[i-1,j] and uv[i,j-1] together and take the difference. another shuffle
   and addition should be enough to get the sum into the lower floating point
   so that it can be extracted and the last two operations to be performed

   the solver function could have a flag added to the parameters since quite a
   bit of the function could be shared between the density version and the
   uv vector version

   though it would be less optimal, it might be less difficult to waste space in
   the uv array and leave the upper two elements of each vector empty; this
   would have the added benefit of making the 3-dimensional version a very
   simple extension, since you'd just have to make sure that the operations
   respect the 3rd element
9. An unintuitive but easily extendable way to store the data is in the order
   {w, z, x, y}. This may actually be the way the SIMD intsructions handle it
   internally but I'm not sure. Initialization methods could still take
   arguments in the intuitive manner, but getting and setting specific values
   would require accessor methods or you'd have to remember that the stored
   order is not normal.

   When stored in this way, vector operations involving vectors of multiple
   lengths can be performed without any extra transformations since the elements
   would already be in the correct locations regardless of size.
```
function advect:
    for i from start to end:
        for j from start to end:
            counter_vec = {0,0,i,j}
            mat_vec = load(uv + IX(i,j))
            mat_vec = counter_vec - dt0 * mat_vec
            bound mat_vec to the size of the array
            ij0 = (int)mat_vec
            ij1 = ij0 + 1

function project:
    for i from start to end:
        for j from start to end:
            
            div[i,j] = 

function solve(uv version):
    for k until acceptable error:
        for j from start to end:
            for i from start to end:
                
```
