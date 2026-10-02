#!/usr/bin/sh
gcc -std=c17 -Wall -Wextra -Isrc -Iinclude -lGL -lglfw \
    -march=x86-64-v3 -O2 -ffast-math -DGLFW_INCLUDE_NONE -o demo \
    src/demo.c src/common_solver.c src/scalar_solver.c src/sse4_2_solver.c \
    src/renderer.c src/gl.c
strip --strip-unneeded demo
