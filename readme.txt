qfind v2.1

qfind is a spaceship search program for Conway's Game of Life and other
Life-like cellular automata written by Matthias Merzenich.  It is based on
David Eppstein's gfind and zdr's zfind search programs.  Additional code and
suggestions were provided by Paul Tooke, Tomas Rokicki, Aidan F. Pierce, and
Adam P. Goucher.

Three scripts are included:

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


get-rows.lua:
   This is a Golly Lua script to help with extending partial results.  Usage
   instructions are provided in the source code.  This script requires the
   Life application Golly.

------------------------------------------------------------------------------
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
      Add option to specify dump file prefix (replaces -d)
      Add some simple checks for parameter validity
   1.2.1, 27 January 2021
      Put functions common to both qfind and qfind-s into a header file
   1.3, 29 January 2021
      Add optional longest partial output at end of search
      Make loading from saved state behave like other options
      Add dump-and-exit option
      Add option to choose depth level of first deepening step
      Replace initial rows parameter with dedicated global variable
      Fix bug with free() of non-malloced pointer in success()
   1.4, 31 January 2021
      Add option to set maximum number of ships in output
      Add option to split the search state (replaces dump-and-exit)
      Change longest partial suppression option from -v to -a
      Fix bug when changing queue size after loading state
      Fix bug due to non-initialization of variables
      Fix bug due to missing return type of loadState()
   2.0, 2 March 2021
      Add ability to save and restore depth-first extensions
      Remove naive search order option
   2.1, 3 May 2021
      Automate whether to enable lookahead caching based on speed
      Change -a and -z options to toggles
      Replace Python script with Lua script
      Fix bug causing redundant output

------------------------------------------------------------------------------
Matthias Merzenich