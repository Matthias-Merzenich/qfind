qfind v1.2.1

qfind is a spaceship search program for Conway's Game of Life and other
Life-like cellular automata.  It is based on David Eppstein's gfind and
zdr's zfind search programs.

Three scripts are provided:

qfind.cpp:
  The main search program.  This program uses OpenMP to achieve rather
  basic parallelization.  It should be compiled with the OpenMP flag
  for your compiler.  I use

  g++ qfind.cpp -O3 -fopenmp -march=native -o qfind


qfind-s.cpp
  A simplified version of qfind with fewer features that is slightly faster.
  You must change the period, offset, and width of the desired search in the
  code before compiling.  This version only allows gcd(period,offset) = 1.
  Compilation is otherwise the same.


get-rows.py:
  This is a Golly Python script to help with extending partial results.
  Usage instructions are provided in the source code.  This script
  requires the Life application Golly.

-----------------------------------------------------------------------------
Version History:
   0.1, 19 June 2017
      Initial release
   0.2, July 2017
      Add mimimum deepening increment parameter
      Add ability to extend partial results
      Make parallel loop scheduling dynamic
   1.0, 3 January 2020
      Add support for non-totalistic rules
      Add lookahead caching
      Add memory limit parameter
      Make table generation dynamic
      Reduce memory usage
      Allow searches at widths greater than 10
   1.0.1, 7 January 2020
      Clean up code
      Make memlimit and cachemem into proper parameters
   1.1, 7 January 2020
      Add option to disable output during deepening step
   1.2, 24 January 2021
      Make options POSIX-like
      Add option to specify dump file prefix
      Add some simple checks for parameter validity
   1.2.1, 27 January 2021
      Put functions common to both qfind and qfind-s into a header file

-----------------------------------------------------------------------------
Matthias Merzenich