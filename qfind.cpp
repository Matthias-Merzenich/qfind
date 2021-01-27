/* qfind
** A spaceship search program by Matthias Merzenich.
** Based on code by David Eppstein, "zdr", Paul Tooke, and Tomas Rokicki.
** Thanks also to Aidan F. Pierce and Adam P. Goucher for code and suggestions.
**
** This is an attempt at combining the functionality of gfind and zfind.
*/

/* Lookahead caching seems to speed up searches for spaceships with  */
/* speeds exceeding c/5 and slow down searches for spaceships with   */
/* speeds at or below c/5.  It is recommended that you compile two   */
/* versions of this program, one unchanged and one with following    */
/* line uncommented to disable lookahead caching.  You can then use  */
/* the version of the program most appropriate for your input.       */

//#define NOCACHE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include <time.h>
#include "common.hpp"

int fwdOff[MAXPERIOD], backOff[MAXPERIOD], doubleOff[MAXPERIOD], tripleOff[MAXPERIOD];

void makePhases(){
   int i;
   for (i = 0; i < period; i++) backOff[i] = -1;
   i = 0;
   for (;;) {
      int j = offset;
      while (backOff[(i+j)%period] >= 0 && j < period) j++;
      if (j == period) {
         backOff[i] = period-i;
         break;
      }
      backOff[i] = j;
      i = (i+j)%period;
   }
   for (i = 0; i < period; i++)
      fwdOff[(i+backOff[i])%period] = backOff[i];
   for (i = 0; i < period; i++) {
      int j = (i - fwdOff[i]);
      if (j < 0) j += period;
      doubleOff[i] = fwdOff[i] + fwdOff[j];
   }
   for (i = 0; i <  period; i++){
      int j = (i - fwdOff[i]);
      if (j < 0) j += period;
      tripleOff[i] = fwdOff[i] + doubleOff[j];
   }
}

int lookAhead(row *pRows, int a, int pPhase){
// indices: first number represents vertical offset,
//          second number represents generational offset
   int ri11, ri12, ri13, ri22, ri23;
   uint16_t *riStart11, *riStart12, *riStart13, *riStart22, *riStart23;
   int numRows11, numRows12, numRows13, numRows22, numRows23;
   int row11, row12, row13, row22, row23;
   int k;
   
   getoffsetcount(pRows[a - params[P_PERIOD] - fwdOff[pPhase]],
                  pRows[a - fwdOff[pPhase]],
                  pRows[a], riStart11, numRows11);
   if (!numRows11)
      return 0;
   
   getoffsetcount(pRows[a - params[P_PERIOD] - doubleOff[pPhase]],
                  pRows[a - doubleOff[pPhase]],
                  pRows[a - fwdOff[pPhase]], riStart12, numRows12);
   
   if(tripleOff[pPhase] >= params[P_PERIOD]){
      riStart13 = pRows + (a + params[P_PERIOD] - tripleOff[pPhase]);
      numRows13 = 1;
#ifndef NOCACHE
      k = getkey(riStart11, riStart12, (uint16_t*)(gcount + riStart13[0]),
         (pRows[a-doubleOff[pPhase]] << width) + pRows[a-tripleOff[pPhase]]);
#endif
   }
   else{
      getoffsetcount(pRows[a - params[P_PERIOD] - tripleOff[pPhase]],
                     pRows[a - tripleOff[pPhase]],
                     pRows[a - doubleOff[pPhase]], riStart13, numRows13);
#ifndef NOCACHE
      k = getkey(riStart11, riStart12, riStart13,
         (pRows[a-doubleOff[pPhase]] << width) + pRows[a-tripleOff[pPhase]]);
#endif
   }
#ifndef NOCACHE
   if (k < 0)
      return k+2;
#endif
   
   for(ri11 = 0; ri11 < numRows11; ++ri11){
      row11 = riStart11[ri11];
      for(ri12 = 0; ri12 < numRows12; ++ri12){
         row12 = riStart12[ri12];
         getoffsetcount(pRows[a - doubleOff[pPhase]],
                        row12, row11, riStart22, numRows22);
         if(!numRows22) continue;
         
         for(ri13 = 0; ri13 < numRows13; ++ri13){
            row13 = riStart13[ri13];
            getoffsetcount(pRows[a - tripleOff[pPhase]],
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
   
   getoffsetcount(pRows[currRow - 2 * period],
                     pRows[currRow - period],
                     pRows[currRow - period + backOff[pPhase]],
                     riStart, numRows) ;

   if(theNode == 0){
      firstRow = 1;
   }
   
   for(i = firstRow; i < numRows; ++i){
      pRows[currRow] = riStart[i];
      if (!isVisited(theNode, pRows[currRow]) && lookAhead(pRows, currRow, pPhase)){
         enqueue(theNode, pRows[currRow]);
         if (terminal(qTail-1)) success(qTail-1, NULL, 0, 0);
         setVisited(qTail - 1);
      }
   }
}

int depthFirst(node theNode, long howDeep, uint16_t **pInd, int *pRemain, row *pRows){
   int pPhase;
   pPhase = peekPhase(theNode);

   node x = theNode;
   uint32_t startRow = 2*period + pPhase + 1;
   uint32_t currRow = startRow;
   
   int i;
   for(i = currRow - 1; i >= 0; --i){
      pRows[i] = ROW(x);
      x = PARENT(x);
   }

   ++pPhase;
   if(pPhase == period) pPhase = 0;
   
   getoffsetcount(pRows[currRow - 2 * period],
                  pRows[currRow - period],
                  pRows[currRow - period + backOff[pPhase]],
                  pInd[currRow], pRemain[currRow]) ;
   pInd[currRow] += pRemain[currRow];
   
   
   
   for(;;){

      /* back up if there are no rows left to check at this depth */
      if(!pRemain[currRow]){
         --currRow;
         if(pPhase == 0) pPhase = period;
         --pPhase;
         if(currRow < startRow)
            return 0;
         
         continue;
      }
      pRows[currRow] = *(pInd[currRow] - pRemain[currRow]);
      --pRemain[currRow];
      if(!lookAhead(pRows, currRow, pPhase)) continue;

      ++currRow;
      ++pPhase;
      if(pPhase == period) pPhase = 0;
      /* Check if we reached the desired depth. If so,  
         check if the result is a complete spaceship */
      if(currRow > startRow + howDeep){
         if(params[P_PRINTDEEP] == 0) return 1;
         for(i = 1; i <= period; ++i){
            if(pRows[currRow - i]) return 1;
         }
         currRow -= period;
         for(i = 1; i<= period; ++i){
            if(causesBirth[pRows[currRow - i]]) return 1;
         }
         /* If we got here then we found a spaceship! */
         #pragma omp critical(printWhileDeepening)
         {
            success(theNode, pRows, startRow - 1, currRow + period - 1);
         }
         
         return 1;
         
      }
      
      getoffsetcount(pRows[currRow - 2 * period],
                     pRows[currRow - period],
                     pRows[currRow - period + backOff[pPhase]],
                     pInd[currRow], pRemain[currRow]) ;
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
   
   params[P_WIDTH] = 0;
   params[P_PERIOD] = 0;
   params[P_OFFSET] = 0;
   params[P_SYMMETRY] = 0;
   params[P_REORDER] = 1;
   params[P_CHECKPOINT] = 0;
   params[P_BASEBITS] = 4;
   params[P_QBITS] = QBITS;
   params[P_HASHBITS] = HASHBITS;
   params[P_NUMTHREADS] = 1;
   params[P_INITROWS] = 0;
   params[P_MINDEEP] = 0;
   params[P_CACHEMEM] = 32;
   params[P_MEMLIMIT] = -1;
   params[P_PRINTDEEP] = 1;
   
   parseOptions(argc, argv);
   
   searchSetup();
   
   printf("Starting search\n");
   fflush(stdout);
   
   breadthFirst();
   
   printf("Search complete.\n");
   
   return 0;
}