This directory contains supplementary material for the paper "Real-Time Fluid Dynamics for Games" by Jos Stam.

This is the content of the directories:

code/

	contains solver.c and demo.c
	solver.c is the basic fluid solver in two dimensions and code.c is a prototype that uses the
	code. It uses OpenGL and GLUT for its interface and should run on many different platforms.

Palm/

	contains fluid.prc
	Simple 2D demo that runs on Palm devices. Simply copy it over to your Palm.

PocketPC2002/

	contains four executable: WindTunnel.exe, PocketPaint.exe, PocketDensity.exe and PocketFluids.exe.
	Copy these over to your PocketPC device to run them. Note that they only run on devices with an
	ARM/Xscale processor and running PocketPC2002.

The demo included in `code/` has been updated to use the OpenGL 3.3 core profile and no longer makes use
of immediate mode. Additionally, more solvers have been added: two vector solvers targeting the SSE2 and
SSE4.2 extensions to the x86 ISA, and an OpenCL-based solver with multiple kernels variants of increasing
complexity and performance.
