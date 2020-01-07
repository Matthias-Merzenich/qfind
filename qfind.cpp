/* qfind v1.0.1
** A spaceship search program by Matthias Merzenich.
** Based on code by David Eppstein, "zdr", Paul Tooke, Tomas Rokicki,.
** Thanks to Aidan F. Pierce, and Adam P. Goucher for code and suggestions.
**
** This is an attempt at combining the functionality of gfind and zfind.
**
** Version History:
** 0.1, 19 June 2017
**    Initial release
** 0.2, July 2017
**    Add mimimum deepening increment parameter
**    Add ability to extend partial results
**    Make parallel loop scheduiing dynamic
** 1.0, 3 January 2020
**    Add support for non-totalistic rules
**    Add lookahead caching
**    Add memory limit parameter
**    Make table generation dynamic
**    Reduce memory usage
**    Allow searches of width greater than 10 
** 1.0.1, 7 January 2020
**    Clean up code
**    Make memlimit and cachemem into proper parameters
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include <time.h>
#include "tab.cpp"

//#define NOCACHE

#define BANNER "qfind v1.0.1 by Matthias Merzenich, 7 January 2020"
#define FILEVERSION ((unsigned long) 2020010701)  /* yyyymmddnn */

#define CHUNK_SIZE 64

#define QBITS 23

#define HASHBITS 21
#define DEFAULT_DEPTHLIMIT (qBits-3)

#define P_WIDTH 0
#define P_PERIOD 1
#define P_OFFSET 2
#define P_SYMMETRY 3
#define P_REORDER 4
#define P_CHECKPOINT 5
#define P_BASEBITS 6
#define P_QBITS 7
#define P_HASHBITS 8 
#define P_DEPTHLIMIT 9 
#define P_NUMTHREADS 10
#define P_INIT_ROWS 11
#define P_MINDEEP 12
#define P_MEMLIMIT 13
#define P_CACHEMEM 14

#define NUM_PARAMS 15

#define SYM_ASYM 1
#define SYM_ODD 2
#define SYM_EVEN 3
#define SYM_GUTTER 4

const char *rule = "B3/S23";
char loadRule[256]; /* used for loading rule from file */

int params[NUM_PARAMS];
int width;
int deepeningAmount;
int lastdeep;
int nRowsInState;
int phase;

int period;
#define MAXPERIOD 30

int fwdOff[MAXPERIOD], backOff[MAXPERIOD], doubleOff[MAXPERIOD], tripleOff[MAXPERIOD];

int offset;

int aborting;
int nFound;

enum Mode {
   asymmetric,          /* basic orthogonal or diagonal pattern */
   odd, even,           /* orthogonal with bilateral symmetry */
   gutter,              /* orthogonal bilateral symmetry with empty column in middle */
} mode;

/* the big data structures */
#define qBits params[P_QBITS]
#define QSIZE (1<<qBits)

#define hashBits params[P_HASHBITS]
#define HASHSIZE (1<<hashBits)
#define HASHMASK (HASHSIZE - 1)

typedef uint32_t node;
typedef uint16_t row;

row * rows;
node * base;
node * hash;

int nttable[512] ;
uint16_t **gInd3 ;
uint32_t *gcount ;
uint16_t *gRows;
long long memusage = 0;
long long memlimit = 0;

#ifndef NOCACHE
long long cachesize ;
struct cacheentry {
   uint16_t *p1, *p2, *p3 ;
   int abn, r ;
} *totalCache ;

struct cacheentry **cache;
#endif

/*
** Representation of vertices.
**
** Each vertex is represented by an entry in the rows[] array.
** That entry consists of the bits of the last row of the vertex's pattern
** concatenated with a number, the "offset" of the vertex.
** The parent of vertex i is formed by adding the offset of i
** with the number in base[i/BASEFACTOR].
**
** If node i has offset -1, the entry is considered to be empty and is skipped.
** This is to deal with the case that base[i/BASEFACTOR] is already set to
** something too far away to deal with via the offset.
**
** qIsEmpty() checks the size of the queue
** enqueue(n,r) adds a new node with last row r and parent n
** dequeue() returns the index of the next unprocessed node
** pop() undoes the previous enqueue operation
** resetQ() clears the queue
**
** ROW(b) returns the last row of b
** PARENT(b) returns the index of the parent of b
*/

#define MAXWIDTH (12)

#define ROWBITS ((1<<width)-1)
#define BASEBITS (params[P_BASEBITS])
#define BASEFACTOR (1<<BASEBITS)
#define MAXOFFSET ((((row) -1) >> width) - 1)

#define ROW(i) (rows[i] & ROWBITS)
#define OFFSET(i) (rows[i] >> width)
#define EMPTY(i) (rows[i] == (row)-1)
#define MAKEEMPTY(i) rows[i] = (row)-1
#define PARENT(i) (base[(i)>>BASEBITS]+OFFSET(i))
#define FIRSTBASE(i) (((i) & ((1<<BASEBITS) - 1)) == 0)

#define MINDEEP  ((params[P_MINDEEP]>0) ? params[P_MINDEEP] : period)

void printRow(row theRow){
   int i;
   for(i = width - 1; i >= 0; i--) printf("%c",(theRow & 1 << i ? 'o' : '.'));
   printf("\n");
}

/* =========================================== */
/*  Lookup Tables to determine successor rows  */
/* =========================================== */

unsigned char *causesBirth;

char nttable2[512] ;

void error(const char *s) {
   fprintf(stderr, "%s\n", s) ;
   exit(10) ;
}

int slowEvolveBit(int row1, int row2, int row3, int bshift){
   return nttable[(((row2>>bshift) & 2)<<7) | (((row1>>bshift) & 2)<<6)
                | (((row1>>bshift) & 4)<<4) | (((row2>>bshift) & 4)<<3)
                | (((row3>>bshift) & 7)<<2) | (((row2>>bshift) & 1)<<1)
                |  ((row1>>bshift) & 1)<<0];
}

void fasterTable() {
   int p = 0 ;
   for (int row1=0; row1<8; row1++)
      for (int row2=0; row2<8; row2++)
         for (int row3=0; row3<8; row3++)
            nttable2[p++] = slowEvolveBit(row1, row2, row3, 0) ;
}

int evolveBit(int row1, int row2, int row3, int bshift) {
   return nttable2[
      (((row1 << 6) >> bshift) & 0700) +
      (((row2 << 3) >> bshift) &  070) +
      (( row3       >> bshift) &   07)] ;
}

int evolveBit(int row1, int row2, int row3) {
   return nttable2[
      ((row1 << 6) & 0700) +
      ((row2 << 3) &  070) +
      ( row3       &   07)] ;
}

int evolveRow(int row1, int row2, int row3){
   int row4;
   int row1_s,row2_s,row3_s;
   int j,s = 0;
   if(params[P_SYMMETRY] == SYM_ODD) s = 1;
   if(evolveBit(row1, row2, row3, width - 1)) return -1;
   if(params[P_SYMMETRY] == SYM_ASYM && evolveBit(row1 << 2, row2 << 2, row3 << 2)) return -1;
   if(params[P_SYMMETRY] == SYM_ODD || params[P_SYMMETRY] == SYM_EVEN){
      row1_s = (row1 << 1) + ((row1 >> s) & 1);
      row2_s = (row2 << 1) + ((row2 >> s) & 1);
      row3_s = (row3 << 1) + ((row3 >> s) & 1);
   }
   else{
      row1_s = (row1 << 1);
      row2_s = (row2 << 1);
      row3_s = (row3 << 1);
   }
   row4 = evolveBit(row1_s, row2_s, row3_s);
   for(j = 1; j < width; j++)row4 += evolveBit(row1, row2, row3, j - 1) << j;
   return row4;
}

int evolveRowHigh(int row1, int row2, int row3, int bits){
   int row4=0;
   int row1_s,row2_s,row3_s;
   int j ;
   if(evolveBit(row1, row2, row3, width - 1)) return -1;
   row1_s = (row1 << 1);
   row2_s = (row2 << 1);
   row3_s = (row3 << 1);
   for(j = width-bits; j < width; j++)row4 += evolveBit(row1, row2, row3, j - 1) << j;
   return row4;
}

int evolveRowLow(int row1, int row2, int row3, int bits){
   int row4;
   int row1_s,row2_s,row3_s;
   int j,s = 0;
   if(params[P_SYMMETRY] == SYM_ODD) s = 1;
   if(params[P_SYMMETRY] == SYM_ASYM && evolveBit(row1 << 2, row2 << 2, row3 << 2)) return -1;
   if(params[P_SYMMETRY] == SYM_ODD || params[P_SYMMETRY] == SYM_EVEN){
      row1_s = (row1 << 1) + ((row1 >> s) & 1);
      row2_s = (row2 << 1) + ((row2 >> s) & 1);
      row3_s = (row3 << 1) + ((row3 >> s) & 1);
   }
   else{
      row1_s = (row1 << 1);
      row2_s = (row2 << 1);
      row3_s = (row3 << 1);
   }
   row4 = evolveBit(row1_s, row2_s, row3_s);
   for(j = 1; j < bits; j++)row4 += evolveBit(row1, row2, row3, j - 1) << j;
   return row4;
}

void sortRows(uint16_t *theRow, uint32_t totalRows) {
   uint32_t i;
   int64_t j;
   uint16_t t;
   for(i = 1; i < totalRows; ++i){
      t = theRow[i];
      j = i - 1;
      while(j >= 0 && gcount[theRow[j]] < gcount[t]){
         theRow[j+1] = theRow[j];
         --j;
      }
      theRow[j+1] = t;
   }
}
uint16_t *makeRow(int row1, int row2) ;

uint16_t *getoffset(int row12) {
   uint16_t *r = gInd3[row12] ;
   if (r == 0)
      r = makeRow(row12 >> width, row12 & ((1 << width) - 1)) ;
   return r ;
}
uint16_t *getoffset(int row1, int row2) {
   return getoffset((row1 << width) + row2) ;
}

void getoffsetcount(int row1, int row2, int row3, uint16_t* &p, int &n) {
   uint16_t *theRow = getoffset(row1, row2) ;
   p = theRow + theRow[row3] ;
   n = theRow[row3+1] - theRow[row3] ;
}
int getcount(int row1, int row2, int row3) {
   uint16_t *theRow = getoffset(row1, row2) ;
   return theRow[row3+1] - theRow[row3] ;
}
int *gWorkConcat ;      /* gWorkConcat to be parceled out between threads */
int *rowHash ;
uint16_t *valorder ;
void genStatCounts() ;
void makeTables() {
   causesBirth = (unsigned char*)malloc((long long)sizeof(*causesBirth)<<width);
   gInd3 = (uint16_t **)calloc(sizeof(*gInd3),(1LL<<(width*2))) ;
   rowHash = (int *)calloc(sizeof(int),(2LL<<(width*2))) ;
   for (int i=0; i<1<<(2*width); i++)
      gInd3[i] = 0 ;
   for (int i=0; i<2<<(2*width); i++)
      rowHash[i] = -1 ;
   gcount = (uint32_t *)calloc(sizeof(*gcount), (1LL << width));
   memusage += (sizeof(*gInd3)+2*sizeof(int)) << (width*2) ;
   uint32_t i;
   for(i = 0; i < 1 << width; ++i) causesBirth[i] = (evolveRow(i,0,0) ? 1 : 0);
   for(i = 0; i < 1 << width; ++i) gcount[i] = 0 ;
   gWorkConcat = (int *)calloc(sizeof(int), (3LL*params[P_NUMTHREADS])<<width);
   if (params[P_REORDER] == 1)
      genStatCounts() ;
   if (params[P_REORDER] == 2)
      for (int i=1; i<1<<width; i++)
         gcount[i] = 1 + gcount[i & (i - 1)] ;
   gcount[0] = 0xffffffff;  /* Maximum value so empty row is chosen first */
   valorder = (uint16_t *)calloc(sizeof(uint16_t), 1LL << width) ;
   for (int i=0; i<1<<width; i++)
      valorder[i] = (1<<width)-1-i ;
   if (params[P_REORDER] != 0)
      sortRows(valorder, 1<<width) ;
   for (int row2=0; row2<1<<width; row2++)
      makeRow(0, row2) ;
}
uint16_t *bbuf ;
int bbuf_left = 0 ;
/* reduce fragmentation by allocating chunks larger than needed and */
/* parceling out the small pieces.                                  */
uint16_t *bmalloc(int siz) {
   if (siz > bbuf_left) {
      bbuf_left = 1 << (2 * width) ;
      memusage += 2*bbuf_left ;
      if (params[P_MEMLIMIT] >= 0 && memusage > memlimit) {
         printf("Aborting due to excessive memory usage\n") ;
         exit(0) ;
      }
      bbuf = (uint16_t *)calloc(sizeof(uint16_t), bbuf_left) ;
   }
   uint16_t *r = bbuf ;
   bbuf += siz ;
   bbuf_left -= siz ;
   return r ;
}
void unbmalloc(int siz) {
   bbuf -= siz ;
   bbuf_left += siz ;
}
unsigned int hashRow(uint16_t *theRow, int siz) {
   unsigned int h = 0 ;
   for (int i=0; i<siz; i++)
      h = h * 3 + theRow[i] ;
   return h ;
}

uint16_t *makeRow(int row1, int row2) {
   int good = 0 ;
   /* Set up gWork for this particular thread */
   int *gWork = gWorkConcat + ((3LL * omp_get_thread_num()) << width);
   int *gWork2 = gWork + (1 << width) ;
   int *gWork3 = gWork2 + (1 << width) ;
   if (width < 4) {
      for (int row3=0; row3<1<<width; row3++)
         gWork3[row3] = evolveRow(row1, row2, row3) ;
   } else {
      int lowbitcount = (width >> 1) + 1 ;
      int hibitcount = ((width + 1) >> 1) + 1 ;
      int hishift = lowbitcount - 2 ;
      int lowcount = 1 << lowbitcount ;
      for (int row3=0; row3<1<<lowbitcount; row3++)
         gWork2[row3] = evolveRowLow(row1, row2, row3, lowbitcount-1) ;
      for (int row3=0; row3<1<<width; row3 += 1<<hishift)
         gWork2[lowcount+(row3>>hishift)] =
                        evolveRowHigh(row1, row2, row3, hibitcount-1) ;
      for (int row3=0; row3<1<<width; row3++)
         gWork3[row3] = gWork2[row3 & ((1<<lowbitcount) - 1)] |
                        gWork2[lowcount+(row3 >> hishift)] ;
   }
   for (int row3i = 0; row3i < 1<<width; row3i++) {
      int row3 = valorder[row3i] ;
      int row4 = gWork3[row3] ;
      if (row4 < 0)
         continue ;
      gWork2[good] = row3 ;
      gWork[good++] = row4 ;
   }
   
   /* bmalloc, unbmalloc, and all operations that read or write to row, */
   /* rowHash, and gInd3 must be included in a critical region.         */
   uint16_t *theRow;
  #pragma omp critical(updateTable)
  {
      theRow = bmalloc((1+(1<<width)+good)) ;
      for (int row3=0; row3 < 1<<width; row3++)
         theRow[row3] = 0 ;
      theRow[0] = 1 + (1 << width) ;
      for (int row3=0; row3 < good; row3++)
         theRow[gWork[row3]]++ ;
      theRow[1<<width] = 0 ;
      for (int row3=0; row3 < (1<<width); row3++)
         theRow[row3+1] += theRow[row3] ;
      for (int row3=good-1; row3>=0; row3--) {
         int row4 = gWork[row3] ;
         theRow[--theRow[row4]] = gWork2[row3] ;
      }
      unsigned int h = hashRow(theRow, 1+(1<<width)+good) ;
      h &= (2 << (2 * width)) - 1 ;
      while (1) {
         if (rowHash[h] == -1) {
            rowHash[h] = (row1 << width) + row2 ;
            break ;
         }
         /* Maybe two different row12s result in the exact same rows for the */
         /* lookup table. This prevents two different threads from trying to */
         /* build the same part of the lookup table.                         */
         if (memcmp(theRow, gInd3[rowHash[h]], 2*(1+(1<<width)+good)) == 0) {
            theRow = gInd3[rowHash[h]] ;
            unbmalloc(1+(1<<width)+good) ;
            break ;
         }
         h = (h + 1) & ((2 << (2 * width)) - 1) ;
      }
      
      gInd3[(row1<<width)+row2] = theRow ;
  }
   
/*
 *   For debugging:
 *
   printf("R") ;
   for (int i=0; i<1+(1<<width)+good; i++)
      printf(" %d", theRow[i]) ;
   printf("\n") ;
   fflush(stdout) ;
 */
 
   return theRow ;
}

/*
**   We calculate the stats using a 2 * 64 << width array.  We use a
**   leading 1 to separate them.  Index 1 aaa bb cc dd represents
**   the count for a result of aaa when the last two bits of row1, row2,
**   and row3 were bb, cc, and dd, respectively.  We have to manage
**   the edge conditions appropriately.
*/
void genStatCounts() {
   int *cnt = (int*)calloc((128 * sizeof(int)), 1LL << width) ;
   for (int i=0; i<128<<width; i++)
      cnt[i] = 0 ;
   int s = 0 ;
   if (params[P_SYMMETRY] == SYM_ODD)
      s = 2 ;
   else if (params[P_SYMMETRY] == SYM_EVEN)
      s = 1 ;
   else
      s = width + 2 ;
   /* left side: never permit generation left of row4 */
   for (int row1=0; row1<2; row1++)
      for (int row2=0; row2<2; row2++)
         for (int row3=0; row3<2; row3++)
            if (evolveBit(row1, row2, row3) == 0)
               cnt[(1<<6) + (row1 << 4) + (row2 << 2) + row3]++ ;
   for (int nb=0; nb<width; nb++) {
      for (int row1=0; row1<8; row1++)
         for (int row2=0; row2<8; row2++)
            for (int row3=0; row3<8; row3++) {
               if (nb == width-1)
                  if ((((row1 >> s) ^ row1) & 1) ||
                      (((row2 >> s) ^ row2) & 1) ||
                      (((row3 >> s) ^ row3) & 1))
                     continue ;
               int row4b = evolveBit(row1, row2, row3) ;
               for (int row4=0; row4<1<<nb; row4++)
                  cnt[(((((1<<nb) + row4) << 1) + row4b) << 6) +
                    ((row1 & 3) << 4) + ((row2 & 3) << 2) + (row3 & 3)] +=
                     cnt[(((1<<nb) + row4) << 6) +
                       ((row1 >> 1) << 4) + ((row2 >> 1) << 2) + (row3 >> 1)] ;
            }
   }
   /* right side; check left, and accumulate into gcount */
   for (int row1=0; row1<4; row1++)
      for (int row2=0; row2<4; row2++)
         for (int row3=0; row3<4; row3++)
            if (params[P_SYMMETRY] != SYM_ASYM ||
                evolveBit(row1<<1, row2<<1, row3<<1) == 0)
               for (int row4=0; row4<1<<width; row4++)
                  gcount[row4] +=
                     cnt[(((1<<width) + row4) << 6) +
                       (row1 << 4) + (row2 << 2) + row3] ;
   free(cnt) ;
}

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

/* ====================================================== */
/*  Hash table for detecting equivalent partial patterns  */
/* ====================================================== */

void resetHash() { if (hash != 0) memset(hash,0,4*HASHSIZE); }

int hashPhase = 0;

static inline long hashFunction(node b, row r) {
   long h = r;
   int i;
   for (i = 0; i < nRowsInState; i++) {
      h = (h * 269) + ROW(b);
      b = PARENT(b);
   }
   h += (h>>16)*269;
   h += (h>>8)*269;
   return h & HASHMASK;
}

/* test if q+r is same as p */
static inline int same(node p, node q, row r) {
   int i;
   for (i = 0; i < nRowsInState; i++) {
      if (p >= QSIZE || q >= QSIZE || EMPTY(p) || EMPTY(q)) return 0;   /* sanity check */
      if (ROW(p) != r) return 0;
      p = PARENT(p);
      r = ROW(q);
      q = PARENT(q);
   }
   return 1;
}

/* test if we've seen this child before */
static inline int isVisited(node b, row r) {
   if (same(0,b,r)) return 1;
   if (hash != 0) {
      int hashVal = hashFunction(b,r);
      node hashNode = hash[hashVal];
      if (hashNode == 0) return 0;
      if (same(hashNode,b,r)){
         return 1;
      }
   }
   return 0;
}

/* set node (NOT child) to be visited */
static inline void setVisited(node b) {
   if (hash != 0) hash[ hashFunction(PARENT(b),ROW(b)) ] = b;
}



/* ============================================== */
/*  Output patterns found by successful searches  */
/* ============================================== */

#define MAXRLELINEWIDTH 63
int RLEcount = 0;
int RLElineWidth = 0;
char RLEchar;

void sendRLE(char c) {
  if (RLEcount > 0 && c != RLEchar) {
      if (RLElineWidth++ >= MAXRLELINEWIDTH) {
         if (RLEchar != '\n') putchar('\n');
         RLElineWidth = 0;
      }
    if (RLEcount == 1) putchar(RLEchar);
    else {
       printf("%d%c", RLEcount, RLEchar);
       RLElineWidth ++;
       if (RLEcount > 9) RLElineWidth++;
    }
    RLEcount = 0;
    if (RLEchar == '\n') RLElineWidth = 0;
  }
  if (c != '\0') {
    RLEcount++;
    RLEchar = c;
  } else RLElineWidth = 0;
}

int outputParity = 0;

void putRow(unsigned long rr, unsigned long r, int shift) {
   while (r | rr) {
      if (shift == 0)
         sendRLE(r & 1 ? 'o' : 'b');
      else shift--;
      r >>= 1;
      if (rr & 1) r |= (1<<31);
      rr >>= 1;
   }
   sendRLE('$');
}

int modeWidth() {
   switch(mode) {
      case asymmetric:
         return width;
      case odd:
         return 2*width-1;
      case even:
         return 2*width;
      case gutter:
         return 2*width+1;
   }
   return 0;
}

/* Avoid Intel shift bug */
static inline unsigned long
safeShift(unsigned long r, int i)
{
    unsigned long rr = r;
    while (i>16) { rr >>= 16; i-=16;}
    return (rr>>i);
}


void success(node b, row *pRows, int nodeRow, uint32_t lastRow){
   node c;
   int nrows = 0;
   int skewAmount = 0;
   int swidth;
   int p, i, j, margin;
   unsigned long *srows, *ssrows, *drows, *ddrows;
   static unsigned long *oldsrows = 0, *oldssrows = 0;
   static unsigned long *olddrows = 0, *oldddrows = 0;
   static int oldnrows = 0;

   uint32_t currRow = lastRow;
   int nDeepRows = 0;
   int nodeDiff;

   /* check if output disabled while deepening */
   /*if (perdidor) {
      perdidor = 2;
      return;
   }*/

   if(pRows != NULL){
      while(pRows[currRow] == 0){
         if(currRow == 0){
            printf("Success called on search root!\n");
            aborting = 1;
            return;
         }
         currRow--;
      }
      nDeepRows = (currRow / period) - 1;
      nodeDiff = nodeRow - period - (currRow%period);
      nodeRow -= nodeDiff;

      for(j = 0; j < nodeDiff; j++){
         b = PARENT(b);
      }
      currRow = currRow - period + 1;
      nrows = nDeepRows;
      
   }
   else{
      /* shift until we find nonzero row.
         then shift PERIOD-1 more times to get leading edge of glider */
      while (ROW(b) == 0) {
         b = PARENT(b);
         if (b == 0) {
            printf("Success called on search root!\n");
            aborting = 1;
            return;
         }
      }
   }
   if(nrows < 0) nrows = 0;
   
   for (p = 0; p < period-1; p++) b = PARENT(b);
   if (b == 0) {
      printf("Success called on search root!\n");
      aborting = 1;
      return;
   }
   
   /* count rows */
   c = b;
   while (c != 0) {
      for (p = 0; p < period; p++)
         c = PARENT(c);
      nrows++;
   }
   
   /* build data structure of rows so we can reduce width etc */
   srows = (unsigned long*)malloc((nrows+MAXWIDTH+1) * sizeof(unsigned long));
   ssrows = (unsigned long*)malloc((nrows+MAXWIDTH+1) * sizeof(unsigned long));
   drows = (unsigned long*)srows; ddrows = (unsigned long*)ssrows; /* save orig ptr for free() */
   for (i = 0; i <= nrows+MAXWIDTH; i++) srows[i]=ssrows[i]=0;
   for (i = nrows - 1; i >= 0; i--) {
      row r;
      if(nDeepRows > 0){
         r = pRows[currRow];
         currRow -= period;
         nDeepRows--;
      }
      else{
         r = ROW(b);
         for (p = 0; p < period; p++) {
            b = PARENT(b);
         }
      }
      //row rx;
      switch(mode) {
         case asymmetric:
            srows[i] = r;
            break;

         case odd:
            srows[i] = r << (MAXWIDTH - 1);
            ssrows[i] = r >> (32 - (MAXWIDTH - 1));
            for (j = 1; j < MAXWIDTH; j++)
               if (r & (1<<j))
                  srows[i] |= 1 << (MAXWIDTH - 1 - j);
            break;

         case even:
            srows[i] = r << MAXWIDTH;
            ssrows[i] = r >> (32 - MAXWIDTH);
            for (j = 0; j < MAXWIDTH; j++)
               if (r & (1<<j))
                  srows[i] |= 1 << (MAXWIDTH - 1 - j);
            break;

         case gutter:
            srows[i] = r << (MAXWIDTH + 1);
            ssrows[i] = r >> (32 - (MAXWIDTH + 1));
            for (j = 0; j < MAXWIDTH; j++)
               if (r & (1<<j))
                  srows[i+skewAmount] |= 1 << (MAXWIDTH - 1 - j);
            break;

         default:
            printf("Unexpected mode in success!\n");
            aborting = 1;
            return;
      }
   }
   
   /* normalize nrows to only include blank rows */
   nrows += MAXWIDTH;
   while (srows[nrows-1] == 0 && ssrows[nrows-1] == 0 && nrows>0) nrows--;
   while (srows[0] == 0 && ssrows[0] == 0 && nrows>0) {
      srows++;
      ssrows++;
      nrows--;
   }
   
   /* sanity check: are all rows empty? */
   int allEmpty = 1;
   for(i = 0; i < nrows; i++){
      if(srows[i]){
         allEmpty = 0;
         break;
      }
   }
   if(allEmpty) return;
   
   /* make at least one row have nonzero first bit */
   i = 0;
   while ((srows[i] & 1) == 0) {
      for (i = 0; (srows[i] & 1) == 0 && i < nrows; i++) { }
      if (i == nrows) {
         for (i = 0; i < nrows; i++) {
            srows[i] >>= 1;
            if (ssrows[i] & 1) srows[i] |= (1<<31);
            ssrows[i] >>= 1;
         }
         i = 0;
      }
   }
   
   swidth = 0;
   for (i = 0; i < nrows; i++)
      while (safeShift(ssrows[i],swidth))
         swidth++;
   if (swidth) swidth += 32;
   for (i = 0; i < nrows; i++)
      while (safeShift(srows[i],swidth))
         swidth++;
   
      
   /* compute margin on other side of width */
   margin = 0;

   /* make sure we didn't just output the exact same pattern (happens a lot for puffer) */
   if (nrows == oldnrows) {
      int different = 0;
      for (i = 0; i < nrows && !different; i++)
         different = (srows[i] != oldsrows[i] || ssrows[i] != oldssrows[i]);
      if (!different) {
         free(drows);
         free(ddrows);
         return;
      }
   }
   if (olddrows != 0) free(olddrows);
   if (oldddrows != 0) free(oldddrows);
   oldsrows = srows;
   oldssrows = ssrows;
   olddrows = drows;
   oldddrows = ddrows;
   oldnrows = nrows;

   /* output it all */
   printf("\nx = %d, y = %d, rule = %s", swidth - margin, nrows, rule);
   putchar('\n');

   while (nrows-- > 0) {
      if (margin > nrows) putRow(ssrows[nrows], srows[nrows], margin - nrows);
      else putRow(ssrows[nrows], srows[nrows], 0);
   }
   RLEchar = '!';
   sendRLE('\0');
   printf("\n\n");
   fflush(stdout);
   //if (++nFound >= findLimit) aborting = 1;
}


/* Is this a node at which we can stop? */
int terminal(node n){
   int p;

   for (p = 0; p < period; p++) {   /* last row in each phase must be zero */
      if (ROW(n) != 0) return 0;
      n = PARENT(n);
   }

   for (p = 0; p < period; p++) {
      if(causesBirth[ROW(n)]) return 0;
      n = PARENT(n);
   }
   return 1;
}




/* ================================================ */
/*  Queue of partial patterns still to be examined  */
/* ================================================ */

/* SU patch */
node qHead,qTail;

/* PATCH queue dimensions required during save/restore */
node qStart; /* index of first node in queue */
node qEnd;   /* index of first unused node after end of queue */

/* Maintain phase of queue nodes. After dequeue(), the global variable phase
   gives the phase of the dequeued item.  If the queue is compacted, this information
   needs to be reinitialized by a call to rephase(), after which phase will not be
   valid until the next call to dequeue().  Variable nextRephase points to the next
   node for which dequeue will need to increment the phase. Phase is not maintained
   when treating queue as a stack (using pop()) -- caller must do it in that case.
   It's ok to change phase since we maintain a separate copy in queuePhase. */

int queuePhase = 0;
node nextRephase = 0;
void rephase() {
   node x, y;
   while (qHead < qTail && EMPTY(qHead)) qHead++;   /* skip past empty queue cells */
   x = qHead;   /* find next item in queue */
   queuePhase = period - 1;
   while (x != 0) {
      x = PARENT(x);
      queuePhase++;
   }
   queuePhase %= period;

   /* now walk forward through queue finding breakpoints between each generation
      invariants: y is always the first in its generation */
   x = 0; y = 0;
   while (y <= qHead) {
      ++x;
      if (x >= qTail || (!EMPTY(x) && PARENT(x) >= y)) y = x;
   }
   nextRephase = y;
}

/* phase of an item on the queue */
int peekPhase(node i) {
   return (i < nextRephase? queuePhase : (queuePhase+1)%period);
}

/* Test queue status */
static inline int qIsEmpty() {
   while (qHead < qTail && EMPTY(qHead)) qHead++;
   return (qTail == qHead);
}

void qFull() {
    if (aborting != 2) {
      printf("Exceeded %d node limit, search aborted\n", QSIZE);
      fflush(stdout);
      aborting = 2;
   }
}

static inline void enqueue(node b, row r) {
   node i = qTail++;
   if (i >= QSIZE) qFull();
   else if (FIRSTBASE(i)) {
      base[i>>BASEBITS] = b;
      rows[i] = r;
   } else {
      long o = b - base[i>>BASEBITS];
      if (o < 0 || o >(long) MAXOFFSET) {   /* offset out of range */
         while (!FIRSTBASE(i)) {
            rows[i] = -1;
            i = qTail++;
            if (i >= QSIZE) qFull();
         }
         base[i>>BASEBITS] = b;
         rows[i] = r;
      } else rows[i] = (o << width) + r;
   }
}

static inline node dequeue() {
   while (qHead < qTail && EMPTY(qHead)) qHead++;
   if (qHead >= nextRephase) {
      queuePhase = (queuePhase+1)%period;
      nextRephase = qTail;
   }
   phase = queuePhase;
   return qHead++;
}

static inline void pop() {
   qTail--;
   while (qTail > qHead && EMPTY(qTail-1)) qTail--;
}

void resetQ() { qHead = qTail = 0; }

static inline int qTop() { return qTail - 1; }


/* ================================= */
/* PATCH08 - dump state              */
/* ================================= */

int dumpNum = 1;
char dumpFile[12];
#define DUMPROOT "dump"
int dumpFlag = 0; /* Dump status flags, possible values follow */
#define DUMPPENDING (1)
#define DUMPFAILURE (2)
#define DUMPSUCCESS (3)

int dumpandexit = 0;

FILE * openDumpFile()
{
    FILE * fp;

    while (dumpNum < 10000)
    {
        sprintf(dumpFile,"%s%04d",DUMPROOT,dumpNum++);
        if ((fp=fopen(dumpFile,"r")))
            fclose(fp);
        else
            return fopen(dumpFile,"w");
    }
    return (FILE *) 0;
}

void dumpState()
{
    FILE * fp;
    int i;
    dumpFlag = DUMPFAILURE;
    if (!(fp = openDumpFile())) return;
    fprintf(fp,"%lu\n",FILEVERSION);
    fprintf(fp,"%s\n",rule);
    for (i = 0; i < NUM_PARAMS; i++)
        fprintf(fp,"%d\n",params[i]);
    fprintf(fp,"%d\n",width);
    fprintf(fp,"%d\n",period);
    fprintf(fp,"%d\n",offset);
    fprintf(fp,"%d\n",mode);
    fprintf(fp,"%d\n",lastdeep);

    fprintf(fp,"%u\n",qHead-qStart);
    fprintf(fp,"%u\n",qEnd-qStart);
    for (i = qStart; i < qEnd; i++)
        fprintf(fp,"%u\n",rows[i]);
    fclose(fp);
    dumpFlag = DUMPSUCCESS;
}





/* ================================= */
/*  Compaction of nearly full queue  */
/* ================================= */

void putnum(long n) {
   char suffix;
   if (n >= 1000000) {
      n /= 100000;
      suffix = 'M';
   } else if (n >= 1000) {
      n /= 100;
      suffix = 'k';
   } else {
      printf("%ld", n);
      return;
   }

   if (n >= 100) printf("%ld", n/10);
   else printf("%ld.%ld", n/10, n%10);
   putchar(suffix);
}

long currentDepth() {
   long i;
   node x;
   x = qTail - 1;
   i = 1;
   while (x != 0) {
      x = PARENT(x);
      i++;
   }
   return i;
}



/*
** doCompact() now has two parts.  The first part compresses
** the queue.  The second part consists of the last loop which
** converts parent bits to back parent pointers. The search
** state may be saved in between.  The queue dimensions, which
** were previously saved in local variables are saved in globals.
*/

void doCompactPart1()
{
   node x,y;
   qEnd = qTail;
   
   /* make a pass backwards from the end finding unused nodes at or before qHead */
   x = qTail - 1;
   y = qHead - 1;
   while (y > 0) {
      /* invariants: everything after y is still active.
                     everything after x points to something after y.
                     x is nonempty and points to y or something before y.
                     so, if x doesnt point to y, y must be unused and can be removed. */
      if (!EMPTY(y)) {
         if (y > PARENT(x)) rows[y] = -1;
         else while (EMPTY(x) || PARENT(x) == y) x--;
      }
      y--;
   }
   
   /* make a pass forwards converting parent pointers to offset from prev parent ptr.
      note that after unused nodes are eliminated, all these offsets are zero or one. */
   y = 0;
   for (x = 0; x < qTail; x++) if (!EMPTY(x)) {
      if (PARENT(x) == y) rows[x] = ROW(x);
      else {
         y = PARENT(x);
         rows[x] = (1<<width) + ROW(x);
      }
   }
   
   /*
      Make a pass backwards compacting gaps.
   
      For most times we run this, it could be combined with the next phase, but
      every once in a while the process of repacking the remaining items causes them
      to use *MORE* space than they did before they were repacked (because of the need
      to leave empty space when OFFSET gets too big) and without this phase the repacked
      stuff overlaps the not-yet-repacked stuff causing major badness.
      
      For this phase, y points to the current item to be repacked, and x points
      to the next free place to pack an item.
    */
   x = y = qTail-1;
   for (;;) {
      if (qHead == y) qHead = x;
      if (!EMPTY(y)) {
         rows[x] = rows[y];
         x--;
      }
      if (y-- == 0) break;   /* circumlocution for while (y >= 0) because x is unsigned */
   }
    qStart = ++x;    /* mark start of queue */
}

void doCompactPart2()
{
    node x,y;

   /*
      Make a pass forwards converting parent bits back to parent pointers.
      
      For this phase, x points to the current item to be repacked, and y points
      to the parent of the previously repacked item.
      After the previous pass, x is initialized to first nonempty item,
      and all items after x are nonempty. 
    */
   qTail = 0; y = 0;
   resetHash();
   for (x = qStart; x < qEnd; x++) {
      if (OFFSET(x)) {   /* skip forward to next parent */
         y++;
         while (EMPTY(y)) y++;
      }
      enqueue(y,ROW(x));
      if (aborting) return;
      if (qHead == x) qHead = qTail - 1;
      setVisited(qTail - 1);
   }
   rephase();
}

void doCompact()
{
   /* make sure we still have something left in the queue */
   if (qIsEmpty()) {
      qTail = qHead = 0;   /* nothing left, make an extremely compact queue */
      return;
   }
    /* First loop of part 1 requires qTail-1 to be non-empty. Make it so */
    while(EMPTY(qTail-1))
        qTail--;

    doCompactPart1();
    if (dumpFlag == DUMPPENDING) dumpState();
    doCompactPart2();
}



/* =================================== */
/* PATCH08 - preview partial results   */
/* =================================== */

static void preview(int allPhases)
{
    node i,j,k;
    int ph;

    for (i = qHead; (i<qTail) && EMPTY(i); i++);
    for (j = qTail-1; (j>i) && EMPTY(j); j--);
    if (j<i) return;
    
    while (j>=i && !aborting)
    {
        if (!EMPTY(j))
        {
            success(j, NULL, 0, 0);
            //success(j);
            if (allPhases == 0)
            {
                k=j;
                for (ph = 1; ph < period; ph++)
                {
                    k=PARENT(k);
                    success(k, NULL, 0, 0);
                    //success(k);
                }
            }
        }
        j--;
    }
}


/* ================================= */
/* PATCH08 - resume from saved state */
/* ================================= */

char * loadFile;

void loadFail()
{
    printf("Load from file %s failed\n",loadFile);
    exit(1);
}

signed int loadInt(FILE *fp)
{
    signed int v;
    if (fscanf(fp,"%d\n",&v) != 1) loadFail();
    return v;
}

unsigned int loadUInt(FILE *fp)
{
    unsigned int v;
    if (fscanf(fp,"%u\n",&v) != 1) loadFail();
    return v;
}

void loadState(char * cmd, char * file)
{
   FILE * fp;
   int i;

   loadFile = file;
   fp = fopen(loadFile, "r");
   if (!fp) loadFail();
   if (loadUInt(fp) != FILEVERSION)
   {
       printf("Incompatible file version\n");
       exit(1);
   }
   
   if (fscanf(fp,"%255s\n",loadRule) != 1) loadFail();
   rule = loadRule;
   
   if (parseRule(rule, nttable) != 0) {
      fprintf(stderr, "Failed to parse rule %s\n", rule) ;
      exit(10) ;
   }
   
   /* Load parameters and set stuff that can be derived from them */
   for (i = 0; i < NUM_PARAMS; i++)
       params[i] = loadInt(fp);

   /* Load / initialise globals */
   width          = loadInt(fp);
   period         = loadInt(fp);
   offset         = loadInt(fp);
   mode           = (Mode)(loadInt(fp));
   lastdeep       = loadInt(fp);
   
   deepeningAmount = period; /* Currently redundant, since it's recalculated */
   //perdidor        = 0;
   aborting        = 0;
   nRowsInState = period+period;   /* how many rows needed to compute successor graph? */

    /* Allocate space for the data structures */
   base = (node*)malloc((QSIZE>>BASEBITS)*sizeof(node));
   rows = (row*)malloc(QSIZE*sizeof(row));
   if (base == 0 || rows == 0) {
      printf("Unable to allocate BFS queue!\n");
      exit(0);
   }
   
   if (hashBits == 0) hash = 0;
   else {
      hash = (node*)malloc(HASHSIZE*sizeof(node));
      if (hash == 0) printf("Unable to allocate hash table, duplicate elimination disabled\n");
   }
   
    /* Load up BFS queue and complete compaction */
    qHead  = loadUInt(fp);
    qEnd   = loadUInt(fp);
    qStart = QSIZE - qEnd;
    qEnd   = QSIZE;
    qHead += qStart;
    if (qStart > QSIZE || qStart < QSIZE/16)
    {
      printf("BFS queue is too small for saved state\n");
      exit(0);
   }
   for (i = qStart; i < qEnd; i++)
        rows[i] = (row) loadUInt(fp);
    fclose(fp);
/*
   printf("qHead:  %d qStart: %d qEnd: %d\n",qHead,qStart,qEnd);
   printf("rows[0]: %d\n",rows[qStart]);
   printf("rows[1]: %d\n",rows[qStart+1]);
   printf("rows[2]: %d\n",rows[qStart+2]);
   fflush(stdout);
   exit(0);
*/
    doCompactPart2();


    /* Let the user know that we got this far */
   printf("State successfully loaded from file %s\n",loadFile);
   
   if(!strcmp(cmd,"p") || !strcmp(cmd,"P")){
      preview(1);
      exit(0);
   }
   
   fflush(stdout);
   
   omp_set_num_threads(params[P_NUMTHREADS]);
}

#ifndef NOCACHE
int getkey(uint16_t *p1, uint16_t *p2, uint16_t *p3, int abn) {
   unsigned long long h = (unsigned long long)p1 +
      17 * (unsigned long long)p2 + 257 * (unsigned long long)p3 +
      513 * abn ;
   h = h + (h >> 15) ;
   h &= (cachesize-1) ;
   struct cacheentry &ce = cache[omp_get_thread_num()][h] ;
   if (ce.p1 == p1 && ce.p2 == p2 && ce.p3 == p3 && ce.abn == abn)
      return -2 + ce.r ;
   ce.p1 = p1 ;
   ce.p2 = p2 ;
   ce.p3 = p3 ;
   ce.abn = abn ;
   return h ;
}

void setkey(int h, int v) {
   cache[omp_get_thread_num()][h].r = v ;
}
#endif

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


static void deepen(){
   node i;
   //node j;

//   if (findLimit > 1) perdidor = 1;   /* disable success if want mult pattern output */

   /* compute amount to deepen, apply reduction if too deep */
#ifdef PATCH07
   timeStamp();
#endif
   printf("Queue full");
   i = currentDepth();
   if (i >= lastdeep) deepeningAmount = MINDEEP;
   else deepeningAmount = lastdeep + MINDEEP - i;   /* go at least MINDEEP deeper */

   lastdeep = i + deepeningAmount;

   /* start report of what's happening */
   printf(", depth %ld, deepening %d, ", (long int) i, deepeningAmount);
   putnum(qTail - qHead);
   printf("/");
   putnum(qTail);
   fflush(stdout);


   /* go through queue, deepening each one */
   
   #pragma omp parallel
   {
   node j;
   uint16_t **pInd;
   int *pRemain;
   row *pRows;
   
   pInd = (uint16_t**)malloc((long long)sizeof(*pInd) * (deepeningAmount + 4 * params[P_PERIOD]));
   pRemain = (int*)malloc((long long)sizeof(*pRemain) * (deepeningAmount + 4 * params[P_PERIOD]));
   pRows = (row*)malloc((long long)sizeof(*pRows) * (deepeningAmount + 4 * params[P_PERIOD]));
   
   #pragma omp for schedule(dynamic, CHUNK_SIZE)
   for (j = qHead; j < qTail; j++) {
      if (!EMPTY(j) && !depthFirst(j, deepeningAmount, pInd, pRemain, pRows))
         MAKEEMPTY(j);
   }
   free(pInd);
   free(pRemain);
   free(pRows);
   }
   
   if (deepeningAmount > period) deepeningAmount--; /* allow for gradual depth reduction */
   
   /* before reporting new queue size, shrink tree back down */
   printf(" -> ");
   fflush(stdout);
   //perdidor = 0;
   
   /* signal time for dump */
   if (params[P_CHECKPOINT]) dumpFlag = DUMPPENDING;\
   
   doCompact();
   
   /* now finish report */
   putnum(qTail - qHead);
   printf("/");
   putnum(qTail);
   printf("\n");
   
   /* Report successful/unsuccessful dump */
   if (dumpFlag == DUMPSUCCESS)
   {
#ifdef PATCH07
      timeStamp();
#endif
       printf("State dumped to %s\n",dumpFile);
       /*analyse();
       if (chainWidth)
           printf("[%d/%d]\n",chainDepth,chainWidth+1);
       else
           printf("[%d/-]\n",chainDepth);*/
   }
   else if (dumpFlag == DUMPFAILURE)
   {
#ifdef PATCH07
      timeStamp();
#endif
      printf("State dump unsuccessful\n");
   }
   
   fflush(stdout);
}


static void breadthFirst()
{
   while (!aborting && !qIsEmpty()) {
      if (qTail - qHead >= (1<<params[P_DEPTHLIMIT]) || qTail >= QSIZE - QSIZE/16 ||
          qTail >= QSIZE - (deepeningAmount << 2)) deepen();
      else process(dequeue());
   }
}


int gcd(int a, int b) {
   if (a > b) return gcd(b,a);
   else if (a == 0) return b;
   else return gcd(b-a,a);
}

void echoParams(){
   printf("Rule: %s\n",rule);
   printf("Period: %d\n",params[P_PERIOD]);
   printf("Offset: %d\n",params[P_OFFSET]);
   printf("Width:  %d\n",params[P_WIDTH]);
   if(params[P_SYMMETRY] == SYM_ASYM) printf("Symmetry: asymmetric\n");
   else if(params[P_SYMMETRY] == SYM_ODD) printf("Symmetry: odd\n");
   else if(params[P_SYMMETRY] == SYM_EVEN) printf("Symmetry: even\n");
   else if(params[P_SYMMETRY] == SYM_GUTTER) printf("Symmetry: gutter\n");
   if(params[P_CHECKPOINT]) printf("Dump state after queue compaction\n");
   if(!params[P_REORDER]) printf("Use naive search order.\n");
   printf("Queue size: 2^%d\n",params[P_QBITS]);
   printf("Hash table size: 2^%d\n",params[P_HASHBITS]);
   printf("Minimum deepening increment: %d\n",MINDEEP);
#ifndef NOCACHE
   printf("Cache memory per thread: %d megabytes\n", params[P_CACHEMEM]);
#endif
   if(params[P_MEMLIMIT] >= 0) printf("Memory limit: %d megabytes\n",params[P_MEMLIMIT]);
   printf("Number of threads: %d\n",params[P_NUMTHREADS]);
}


void usage(){
   printf("%s\n",BANNER);
   printf("\n");
   printf("Usage: \"qfind options\"\n");
   printf("  e.g. \"qfind B3/S23 p3 k1 w6 v\" searches Life (rule B3/S23) for\n");
   printf("  c/3 orthogonal spaceships with even bilateral symmetry and a\n");
   printf("  search width of 6 (full width 12).\n");
   printf("\n");
   printf("Available options:\n");
   printf("  bNN/sNN searches for spaceships in the specified rule (default: b3/s23)\n");
   printf("          Non-totalistic rules can be entered using Hensel notation.\n");
   printf("\n");
   printf("  pNN  searches for spaceships with period NN\n");
   printf("  kNN  searches for spaceships that travel NN cells every period\n");
   printf("  wNN  searches for spaceships with search width NN\n");
   printf("       (full width depends on symmetry type)\n");
   printf("\n");
   printf("  a    searches for asymmetric spaceships\n");
   printf("  u    searches for odd bilaterally symmetric spaceships\n");
   printf("  v    searches for even bilaterally symmetric spaceships\n");
   printf("  g    searches for symmetric spaceships with gutters (empty center column)\n");
   printf("\n");
   printf("  tNN  runs search using NN threads during deepening step (default: 1)\n");
   printf("  hNN  sets the hash table size to 2^NN (default: %d)\n",HASHBITS);
   printf("       Use h0 to disable duplicate elimination.\n");
   printf("  qNN  sets the BFS queue size to 2^NN (default: %d)\n",QBITS);
   printf("  iNN  groups 2^NN queue entries to an index node (default: 4)\n");
#ifndef NOCACHE
   printf("  cNN  allocates NN megabytes per thread for lookahead cache (default: 32)\n");
#endif
   printf("  rNN  limits memory usage to NN megabytes (default: no limit)\n");
   printf("\n");
   printf("  mNN  sets minimum deepening increment to NN (default: period)\n");
   printf("\n");
   printf("  d    dumps the search state after each queue compaction\n");
   //printf("  j    dumps the state at start of search\n");
   printf("\n");
   printf("  o    uses naive search order (not recommended)\n");
   printf("\n");
   printf("  e FF uses rows in the file FF as the initial rows for the search\n");
   printf("       (use the companion Golly python script to easily generate the\n");
   printf("       initial row file)\n");
   printf("\n");
   printf("\"qfind command file\" reloads the state from the specified file\n");
   printf("and performs the command. Available commands: \n");
   printf("  s    resumes search from the loaded state\n");
   printf("  p    previews partial results\n");
}

void loadInitRows(char * file){
   FILE * fp;
   int i,j;
   char rowStr[MAXWIDTH];
   row theRow = 0;
   
   loadFile = file;
   fp = fopen(loadFile, "r");
   if (!fp) loadFail();
   
   printf("Starting search from rows in %s:\n",loadFile);
   
   for(i = 0; i < 2 * period; i++){
      fscanf(fp,"%s",rowStr);
      for(j = 0; j < width; j++){
         theRow |= ((rowStr[width - j - 1] == '.') ? 0:1) << j;
      }
      printRow(theRow);
      enqueue(dequeue(),theRow);
      theRow = 0;
   }
   fclose(fp);
}

int main(int argc, char *argv[]){
   printf("%s\n",BANNER);
   printf("-") ;
   for (int i=0; i<argc; i++)
      printf(" %s", argv[i]) ;
   printf("\n");
   
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
   params[P_INIT_ROWS] = 0;
   params[P_MINDEEP] = 0;
   params[P_CACHEMEM] = 32;
   params[P_MEMLIMIT] = -1;
   
   int loadDumpFlag = 0;
   const char *err ;

   //int dumpandexit = 0;
   int s;
   if(argc == 2 && !strcmp(argv[1],"c")){
      usage();
      return 0;
   }
   parseRule(rule, nttable) ; /* pick up default rule */
   if(argc == 3 && (!strcmp(argv[1],"s") || !strcmp(argv[1],"S") || !strcmp(argv[1],"p") || !strcmp(argv[1],"P"))) loadDumpFlag = 1;
   else{
      for(s = 1; s < argc; s++){    /* read input parameters */
         switch(argv[s][0]){
            case 'b': case 'B':     /* read rule */
               rule = argv[s] ;
               err = parseRule(argv[s], nttable) ;
               if (err != 0) {
                  fprintf(stderr, "Failed to parse rule %s\n", argv[s]) ;
                  exit(10) ;
               }
            break;
            case 'w': case 'W': sscanf(&argv[s][1], "%d", &params[P_WIDTH]); break;
            case 'p': case 'P': sscanf(&argv[s][1], "%d", &params[P_PERIOD]); break;
            case 'k': case 'K': sscanf(&argv[s][1], "%d", &params[P_OFFSET]); break;
            case 'u': case 'U': params[P_SYMMETRY] = SYM_ODD; mode = odd; break;
            case 'v': case 'V': params[P_SYMMETRY] = SYM_EVEN; mode = even; break;
            case 'a': case 'A': params[P_SYMMETRY] = SYM_ASYM; mode = asymmetric; break;
            case 'g': case 'G': params[P_SYMMETRY] = SYM_GUTTER; mode = gutter; break;
            case 'd': case 'D': params[P_CHECKPOINT] = 1; break;
            //case 'j': case 'J': dumpandexit = 1; break;
            case 'e': case 'E': params[P_INIT_ROWS] = ++s; break;
            case 'm': case 'M': sscanf(&argv[s][1], "%d", &params[P_MINDEEP]); break;
            case 't': case 'T': sscanf(&argv[s][1], "%d", &params[P_NUMTHREADS]); break;
            case 'o': case 'O': params[P_REORDER] = 0; break;
            case 'q': case 'Q': sscanf(&argv[s][1], "%d", &params[P_QBITS]); break;
            case 'h': case 'H': sscanf(&argv[s][1], "%d", &params[P_HASHBITS]); break;
            case 'i': case 'I': sscanf(&argv[s][1], "%d", &params[P_BASEBITS]); break;
            case 'c': case 'C': sscanf(&argv[s][1], "%d", &params[P_CACHEMEM]); break;
            case 'r': case 'R': sscanf(&argv[s][1], "%d", &params[P_MEMLIMIT]); break;
            default:
               printf("Unrecognized option %s\n", argv[s]) ;
               exit(10);
         }
      }
   }
   
   if(loadDumpFlag) loadState(argv[1],argv[2]);
   else{
      width = params[P_WIDTH];
      period = params[P_PERIOD];
      offset = params[P_OFFSET];
      deepeningAmount = period;
      lastdeep = 0;
      hashPhase = (gcd(period,offset)>1);
      
      nRowsInState = period+period;
      
      params[P_DEPTHLIMIT] = DEFAULT_DEPTHLIMIT;
      
      base = (node*)malloc((QSIZE>>BASEBITS)*sizeof(node));
      rows = (row*)malloc(QSIZE*sizeof(row));
      if (base == 0 || rows == 0) {
         printf("Unable to allocate BFS queue!\n");
         exit(0);
      }
   
      if (hashBits == 0) hash = 0;
      else {
         hash = (node*)malloc(HASHSIZE*sizeof(node));
         if (hash == 0) printf("Unable to allocate hash table, duplicate elimination disabled\n");
      }
      
      omp_set_num_threads(params[P_NUMTHREADS]);
      
      enqueue(0,0);
      
      if(params[P_INIT_ROWS]) loadInitRows(argv[params[P_INIT_ROWS]]);
   }
   
   memlimit = ((long long)params[P_MEMLIMIT]) << 20;
   
   /* Allocate lookahead cache */
#ifndef NOCACHE
   cachesize = 32768 ;
   while (cachesize * sizeof(cacheentry) < 550000 * params[P_CACHEMEM])
      cachesize <<= 1 ;
   memusage += sizeof(cacheentry) * (cachesize + 5) * params[P_NUMTHREADS];
   if(params[P_MEMLIMIT] >= 0 && memusage > memlimit){
      printf("Not enough memory to allocate lookahead cache\n");
      exit(0);
   }
   totalCache = (struct cacheentry *)calloc(sizeof(cacheentry),
         (cachesize + 5) * params[P_NUMTHREADS]) ;
   cache = (struct cacheentry **)calloc(sizeof(**cache), params[P_NUMTHREADS]);
   
   for(int i = 0; i < params[P_NUMTHREADS]; i++)
      cache[i] = totalCache + (cachesize + 5) * i;
#endif
   
   echoParams();
   
   makePhases();
   fasterTable();
   makeTables();
   
   rephase();

   printf("Starting search\n");
   fflush(stdout);
   
   breadthFirst();
   
   printf("Search complete.\n");

   return 0;
}