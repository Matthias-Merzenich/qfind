/* qfind-simple
** A Simplified version of qfind that is more limited but slightly faster.
** You must set the desired period, offset, and width before compiling.
*/

/*  Lookahead caching seems to speed up searches for spaceships with    */
/*  speeds exceeding c/5 and slow down searches for spaceships with     */
/*  speeds at or below c/5.  If appropriate for your input parameters,  */
/*  uncomment the following line to disable lookahead caching.          */

//#define NOCACHE

/*  Change the following three values before compiling.  */
/*  You must have gcd(PERIOD,OFFSET) = 1.                */
#define PERIOD 5
#define OFFSET 2
#define WIDTH 8

#define QSIMPLE

#include "common.hpp"

#if WIDTH < 1 || PERIOD < 1 || OFFSET < 1
   #error "Invalid value for PERIOD, OFFSET, or WIDTH."
#endif

#if PERIOD > MAXPERIOD
   #error "maximum allowed PERIOD exceeded."
#endif

#if OFFSET > PERIOD && PERIOD > 0
   #error "OFFSET cannot exceed PERIOD."
#endif

#if OFFSET == PERIOD && PERIOD > 0
   #error "Photons are not supported."
#endif

#if 2 * OFFSET > PERIOD && PERIOD > 0
   #warning "Searches for speeds exceeding c/2 may not work correctly."
#endif

#ifdef NOCACHE
   #if 5 * OFFSET > PERIOD && PERIOD > 0
      #warning "Searches for speeds exceeding c/5 may be slower without caching. It is recommended that you recompile with NOCACHE undefined."
   #endif
#else
   #if 5 * OFFSET <= PERIOD && OFFSET > 0
      #warning "Searches for speeds at or below c/5 may be slower with caching. It is recommended that you recompile with NOCACHE defined."
   #endif
#endif

int lookAhead(row *pRows, int a){
/* indices: first digit represents vertical offset,      */
/*          second digit represents generational offset  */
   int ri11, ri12, ri13, ri22, ri23;
   uint16_t *riStart11, *riStart12, *riStart13, *riStart22, *riStart23;
   int numRows11, numRows12, numRows13, numRows22, numRows23;
   int row11, row12, row13, row22, row23;
   int k;
   
   getoffsetcount(pRows[a - PERIOD - OFFSET],
                  pRows[a - OFFSET],
                  pRows[a], riStart11, numRows11);
   if (!numRows11)
      return 0;
   
   getoffsetcount(pRows[a - PERIOD - 2*OFFSET],
                  pRows[a - 2*OFFSET],
                  pRows[a - OFFSET], riStart12, numRows12);
   
   #if 3*OFFSET >= PERIOD
      riStart13 = pRows + (a + PERIOD - 3*OFFSET);
      numRows13 = 1;
#ifndef NOCACHE
      k = getkey(riStart11, riStart12, (uint16_t*)(gcount + riStart13[0]),
         (pRows[a-2*OFFSET] << width) + pRows[a-3*OFFSET]);
#endif
   #else
      getoffsetcount(pRows[a - PERIOD - 3*OFFSET],
                     pRows[a - 3*OFFSET],
                     pRows[a - 2*OFFSET], riStart13, numRows13);
#ifndef NOCACHE
      k = getkey(riStart11, riStart12, riStart13,
         (pRows[a-2*OFFSET] << width) + pRows[a-3*OFFSET]);
#endif
   #endif
#ifndef NOCACHE
   if (k < 0)
      return k+2;
#endif
   
   for(ri11 = 0; ri11 < numRows11; ++ri11){
      row11 = riStart11[ri11];
      for(ri12 = 0; ri12 < numRows12; ++ri12){
         row12 = riStart12[ri12];
         getoffsetcount(pRows[a - 2*OFFSET],
                        row12, row11, riStart22, numRows22);
         if(!numRows22) continue;
         
         for(ri13 = 0; ri13 < numRows13; ++ri13){
            row13 = riStart13[ri13];
            getoffsetcount(pRows[a - 3*OFFSET],
                           row13, row12, riStart23, numRows23);
            if(!numRows23) continue;
            
            for(ri23 = 0; ri23 < numRows23; ++ri23){
               row23 = riStart23[ri23];
               uint16_t *p = getoffset(row13, row23);
               for(ri22 = 0; ri22 < numRows22; ++ri22){
                  row22 = riStart22[ri22];
                  if (p[row22+1]!=p[row22]) {
#ifndef NOCACHE
                     setkey(k, 1);
#endif
                     return 1;
                  }
               }
            }
         }
      }
   }
#ifndef NOCACHE
   setkey(k, 0);
#endif
   return 0;
}

void process(node theNode)
{
   long long int i;
   int firstRow = 0;
   int numRows;
   uint32_t newRowSet;
   node x = theNode;
   int pPhase = peekPhase(x);
   row *riStart;
   row pRows[3*MAXPERIOD];
   int currRow = 2*period + pPhase + 1;
   for(i = currRow - 1; i >= 0; --i){
      pRows[i] = ROW(x);
      x = PARENT(x);
   }

   ++pPhase;
   if(pPhase == period) pPhase = 0;
   
   getoffsetcount(pRows[currRow - 2 * PERIOD],
                     pRows[currRow - PERIOD],
                     pRows[currRow - PERIOD + OFFSET],
                     riStart, numRows);

   if(theNode == 0){
      firstRow = 1;
   }
   
   for(i = firstRow; i < numRows; ++i){
      pRows[currRow] = riStart[i];
      if (!isVisited(theNode, pRows[currRow]) && lookAhead(pRows, currRow)){
         enqueue(theNode, pRows[currRow]);
         if(currentDepth() > longest){
            if(params[P_LONGEST]) bufferPattern(qTail-1, NULL, 0, 0, 0);
            longest = currentDepth();
         }
         if (terminal(qTail-1)) success(qTail-1, NULL, 0, 0);
         setVisited(qTail - 1);
      }
   }
}

int depthFirst(node theNode, long howDeep, uint16_t **pInd, int *pRemain, row *pRows){
   node x = theNode;
   uint32_t startRow = 2*PERIOD + 1;
   uint32_t currRow = startRow;
   
   int i;
   for(i = currRow - 1; i >= 0; --i){
      pRows[i] = ROW(x);
      x = PARENT(x);
   }
   
   getoffsetcount(pRows[currRow - 2 * PERIOD],
                  pRows[currRow - PERIOD],
                  pRows[currRow - PERIOD + OFFSET],
                  pInd[currRow], pRemain[currRow]);
   pInd[currRow] += pRemain[currRow];
   
   
   
   for(;;){

      /* back up if there are no rows left to check at this depth */
      if(!pRemain[currRow]){
         --currRow;
         if(currRow < startRow)
            return 0;
         
         continue;
      }
      pRows[currRow] = *(pInd[currRow] - pRemain[currRow]);
      --pRemain[currRow];
      if(!lookAhead(pRows, currRow)) continue;

      ++currRow;
      /* Check if we reached the desired depth. If so,  
         check if the result is a complete spaceship */
      if(currRow > startRow + howDeep){
         if(params[P_PRINTDEEP] == 0) return 1;
         for(i = 1; i <= PERIOD; ++i){
            if(pRows[currRow - i]) return 1;
         }
         currRow -= PERIOD;
         for(i = 1; i<= PERIOD; ++i){
            if(causesBirth[pRows[currRow - i]]) return 1;
         }
         /* If we got here, then we found a spaceship! */
         #pragma omp critical(printWhileDeepening)
         {
            success(theNode, pRows, startRow - 1, currRow + PERIOD - 1);
         }
         return 1;
      }
      
      getoffsetcount(pRows[currRow - 2 * PERIOD],
                     pRows[currRow - PERIOD],
                     pRows[currRow - PERIOD + OFFSET],
                     pInd[currRow], pRemain[currRow]);
      pInd[currRow] += pRemain[currRow];
   }
}

int main(int argc, char *argv[]){
   printf("%s\n",BANNER);
   printf("Input:");
   for (int i=0; i<argc; i++)
      printf(" %s", argv[i]);
   printf("\n\n");
   
   parseRule(rule, nttable); /* pick up default rule */
   
   params[P_WIDTH] = WIDTH;
   params[P_PERIOD] = PERIOD;
   params[P_OFFSET] = OFFSET;
   params[P_SYMMETRY] = 0;
   params[P_REORDER] = 1;
   params[P_CHECKPOINT] = 0;
   params[P_BASEBITS] = 4;
   params[P_QBITS] = QBITS;
   params[P_HASHBITS] = HASHBITS;
   params[P_NUMTHREADS] = 1;
   params[P_MINDEEP] = 0;
   params[P_CACHEMEM] = 32;
   params[P_MEMLIMIT] = -1;
   params[P_PRINTDEEP] = 1;
   params[P_LONGEST] = 1;
   
   parseOptions(argc, argv);
   
   searchSetup();
   
   printf("Starting search\n");
   fflush(stdout);
   
   breadthFirst();
   
   finalReport();
   
   return 0;
}