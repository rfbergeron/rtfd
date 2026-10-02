#!/usr/bin/sh
gcc -std=c17 -Wall -Wextra -Isrc -Iinclude -lGL -lglfw -lOpenCL \
    -march=x86-64-v3 -O0 -p -g -fsanitize=address,undefined \
    -DGLFW_INCLUDE_NONE -o demo src/demo.c src/common_solver.c \
    src/scalar_solver.c src/sse4_2_solver.c src/cl_solver.c \
    src/cl_helper.c src/gl.c src/renderer.c
