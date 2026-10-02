#!/usr/bin/sh
clang-format -i --style=google src/demo.c src/common_solver.c \
    src/common_solver.h src/scalar_solver.c src/scalar_solver.h \
    src/sse4_2_solver.c src/sse4_2_solver.h src/renderer.c src/renderer.h \
    src/cl_solver.c src/cl_solver.h src/cl_helper.h src/cl_helper.c \
    src/solver.cl src/cl_common.h
