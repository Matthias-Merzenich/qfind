qfind is a spaceship search program for Conway's Game of Life and other
Life-like cellular automata.  It is based on David Eppstein's gfind and
zdr's zfind search programs.

Two files are provided:

qfind.c:
  The main search program.  This program uses OpenMP to achieve rather
  basic parallelization.  It should be compiled with the OpenMP flag
  for your compiler.  For GCC use

  gcc -O3 -fopenmp qfind.c -o qfind


get-rows.py:
  This is a Golly Python script to help with extending partial results.
  Usage instructions are provided in the source code.  This script
  requires the Life simulation application Golly.


--
Matthias Merzenich