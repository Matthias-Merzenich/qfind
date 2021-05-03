/* qfind-simple
** A Simplified version of qfind that is more limited but slightly faster.
** You must set the desired period, offset, and width before compiling.
*/

/*  Lookahead caching seems to speed up searches for spaceships with    */
/*  speeds exceeding c/5 and slow down searches for spaceships with     */
/*  speeds at or below c/5.  Lookahead caching will automatically be    */
/*  enabled or disabled depending on your choice of PERIOD and OFFSET.  */
/*  To override the automatic choice, uncomment one of the following    */
/*  two lines.                                                          */

//#define NOCACHE
//#define FORCECACHE

/*  Change the following three values before compiling.  */
/*  You must have gcd(PERIOD,OFFSET) = 1.                */
#define PERIOD 5
#define OFFSET 2
#define WIDTH 8

#define QSIMPLE

#if WIDTH < 1 || PERIOD < 1 || OFFSET < 1
   #error "Invalid value for PERIOD, OFFSET, or WIDTH."
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

#if defined (NOCACHE)
   #if 5 * OFFSET > PERIOD && PERIOD > 0
      #warning "Searches for speeds exceeding c/5 may be slower without caching. It is recommended that you recompile with NOCACHE undefined."
   #endif
#elif defined (FORCECACHE)
   #if 5 * OFFSET <= PERIOD && OFFSET > 0
      #warning "Searches for speeds at or below c/5 may be slower with caching. It is recommended that you recompile with NOCACHE defined."
   #endif
#endif

#ifndef FORCECACHE
   #if 5 * OFFSET <= PERIOD
      #define NOCACHE
   #endif
#endif

#include "common.hpp"

#if PERIOD > MAXPERIOD
   #error "maximum allowed PERIOD exceeded."
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
   int matchFlag = 1;
   node x = theNode;
   row *riStart;
   row pRows[2*MAXPERIOD + 2];
   int currRow = 2*period + 1;
   for(i = currRow - 1; i >= 0; --i){
      pRows[i] = ROW(x);
      x = PARENT(x);
   }
   
   getoffsetcount(pRows[currRow - 2 * PERIOD],
                     pRows[currRow - PERIOD],
                     pRows[currRow - PERIOD + OFFSET],
                     riStart, numRows);
   
   /* we just ran dequeue() so we need to look at the previous head location */
   uint32_t deepIndex = deepRowIndices[oldDeepQHead];
   
   if(theNode == 0){
      firstRow = 1;
   }
   else if(deepIndex > 1){  /* This means we found an extension previously */
      
      /* Sanity check: do the extension rows match the node rows? */
      node y = theNode;
      for(i = 0; i < 2*period; ++i){
         uint16_t startRow = deepRows[deepIndex][1] + 1;
         if(deepRows[deepIndex][startRow - i] != ROW(y)){
            fprintf(stderr, "Warning: non-matching rows detected at node %u in process()\n",theNode);
            matchFlag = 0;
            free(deepRows[deepIndex]);
            deepRows[deepIndex] = 0;
            break;
         }
         y = PARENT(y);
      }
      
      if(matchFlag){
         uint16_t deepStart = deepRows[deepIndex][1] + 2;
         ++deepRows[deepIndex][1];
         
         while(riStart[firstRow] != deepRows[deepIndex][deepStart]) ++firstRow;
         if (!isVisited(theNode, riStart[firstRow])){
            enqueue(theNode, riStart[firstRow]);
            deepRowIndices[deepQTail - 1] = deepIndex;
            if(currentDepth() > longest){
               if(params[P_LONGEST]) bufferPattern(qTail-1, NULL, 0, 0, 0);
               longest = currentDepth();
            }
            if (terminal(qTail-1) && !terminal(PARENT(qTail-1))) success(qTail-1, NULL, 0, 0);
            setVisited(qTail - 1);
            if(deepRows[deepIndex][1] > deepRows[deepIndex][0]){
               deepRowIndices[deepQTail - 1] = 0;
            }
         }
         
         /* eliminate extension if it gets too short */
         if(deepRows[deepIndex][1] > deepRows[deepIndex][0]){
            free(deepRows[deepIndex]);
            deepRows[deepIndex] = 0;
         }
         ++firstRow;
      }
   }
   
   /* This is no longer in the queue, so we can clear it */
   deepRowIndices[oldDeepQHead] = 0;
   
   for(i = firstRow; i < numRows; ++i){
      pRows[currRow] = riStart[i];
      if (!isVisited(theNode, pRows[currRow]) && lookAhead(pRows, currRow)){
         enqueue(theNode, pRows[currRow]);
         if(currentDepth() > longest){
            if(params[P_LONGEST]) bufferPattern(qTail-1, NULL, 0, 0, 0);
            longest = currentDepth();
         }
         if (terminal(qTail-1) && !terminal(PARENT(qTail-1))) success(qTail-1, NULL, 0, 0);
         setVisited(qTail - 1);
      }
   }
}

int reloadDepthFirst(uint16_t startRow, uint16_t howDeep, row *pRows, uint16_t **pIndGen, int *pRemainGen, row *pRowsGen){
   uint16_t currRow = startRow;
   
   /* Return value if length of extension is greater than deepening amount */
   if(pRows[0] >= howDeep + pRows[1]) return 1;
   
   memcpy(pRowsGen + startRow, pRows + 2 + pRows[1], (pRows[0] - pRows[1] + 1) * (long long)sizeof(*pRowsGen));
   
   for(currRow = startRow; currRow <= startRow + (pRows[0] - pRows[1]); ++currRow){
      getoffsetcount(pRowsGen[currRow - 2 * PERIOD],
                  pRowsGen[currRow - PERIOD],
                  pRowsGen[currRow - PERIOD + OFFSET],
                  pIndGen[currRow], pRemainGen[currRow]);
      
      pIndGen[currRow] += pRemainGen[currRow];
      
      while(*(pIndGen[currRow] - pRemainGen[currRow]) != pRowsGen[currRow]){
         --pRemainGen[currRow];
      }
      --pRemainGen[currRow];
   }
   return 0;
}

int depthFirst(node theNode, uint16_t howDeep, uint16_t **pInd, int *pRemain, row *pRows){
   node x = theNode;
   uint32_t startRow = 2*PERIOD + 1;
   uint32_t currRow = startRow;
   uint32_t theDeepIndex = deepRowIndices[deepQHead + theNode - qHead];
   int matchFlag = 1;
   
   int i;
   for(i = currRow - 1; i >= 0; --i){
      pRows[i] = ROW(x);
      x = PARENT(x);
   }
   
   /* Reload state if we have a previous extension */
   if(theDeepIndex > 1){
      if(reloadDepthFirst( startRow,
                           howDeep,
                           deepRows[theDeepIndex],
                           pInd,
                           pRemain,
                           pRows ))
      {
         return 1;   /* return if howDeep is less than the length of the previous extension */
      }
      
      /* Sanity check: do the extension rows match the node rows? */
      node y = theNode;
      for(i = 0; i < 2*PERIOD; ++i){
         if(deepRows[theDeepIndex][deepRows[theDeepIndex][1] + 1 - i] != ROW(y)){
            fprintf(stderr, "Warning: non-matching rows detected at node %u in depthFirst()\n",theNode);
            matchFlag = 0;
            break;
         }
         y = PARENT(y);
      }
      
      if(matchFlag){
         currRow = startRow   + deepRows[theDeepIndex][0]
                              - deepRows[theDeepIndex][1]
                              + 1;
         
         #pragma omp critical(findDeepIndex)
         {
            free(deepRows[theDeepIndex]);
            deepRows[theDeepIndex] = 0;
         }
      }
   }
   
   deepRowIndices[deepQHead + theNode - qHead] = 0;
   
   getoffsetcount(pRows[currRow - 2 * PERIOD],
                  pRows[currRow - PERIOD],
                  pRows[currRow - PERIOD + OFFSET],
                  pInd[currRow], pRemain[currRow]);
   pInd[currRow] += pRemain[currRow];
   
   for(;;){
      /* back up if there are no rows left to check at this depth */
      if(!pRemain[currRow]){
         --currRow;
         if(currRow < startRow) return 0;
         
         continue;
      }
      pRows[currRow] = *(pInd[currRow] - pRemain[currRow]);
      --pRemain[currRow];
      if(!lookAhead(pRows, currRow)) continue;

      ++currRow;
      /* Check if we reached the desired depth. If so,  
         check if the result is a complete spaceship */
      if(currRow > startRow + howDeep){
         /* Flag that an extension was found. This value will be changed by saveDepthFirst() */
         deepRowIndices[deepQHead + theNode - qHead] = 1;
         
         /* Save the extension if it is long enough */
         if(howDeep >= params[P_MINEXTENSION]){
            saveDepthFirst(theNode, startRow, howDeep, pRows);
         }
         
         /* Check if the extension represents a spaceship */
         if(params[P_PRINTDEEP] == 0) return 1;
         for(i = 1; i <= PERIOD; ++i){
            if(pRows[currRow - i]) return 1;
         }
         currRow -= PERIOD;
         for(i = 1; i <= PERIOD; ++i){
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
   
   setDefaultParams();
   
   parseOptions(argc, argv);
   
   searchSetup();
   
   printf("Starting search\n");
   fflush(stdout);
   
   breadthFirst();
   
   finalReport();
   
   return 0;
}