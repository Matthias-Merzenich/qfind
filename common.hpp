/* This header file contains functions common to both qfind and qfind-s.
** Some functions contained in this file work slightly differently depending
** on whether qfind or or qfind-s is being compiled.  Such differences are
** determined by the presence of the macro QSIMPLE defined in qfind-s.cpp.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include <time.h>

#ifdef QSIMPLE
   #define WHICHPROGRAM qfind-simple
#else
   #define WHICHPROGRAM qfind
#endif

#define STR(x) #x
#define XSTR(x) STR(x)

#define BANNER XSTR(WHICHPROGRAM)" v1.4 by Matthias Merzenich, 31 January 2021"

#define FILEVERSION ((unsigned long) 2021013101)  /* yyyymmddnn */

#define MAXPERIOD 30
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
#define P_MINDEEP 11
#define P_MEMLIMIT 12
#define P_CACHEMEM 13
#define P_PRINTDEEP 14
#define P_LONGEST 15
#define P_LASTDEEP 16
#define P_NUMSHIPS 17

#define NUM_PARAMS 18

#define SYM_ASYM 1
#define SYM_ODD 2
#define SYM_EVEN 3
#define SYM_GUTTER 4

const char *rule = "B3/S23";
char loadRule[256]; /* used for loading rule from file */

char *initRows;

int params[NUM_PARAMS];
int width;
int deepeningAmount;
int nRowsInState;
int phase;

int period;

int offset;

int aborting;
int numFound = 0;       /* number of spaceships found so far */
int longest = 0;        /* length of current longest partial result */

enum Mode {
   asymmetric,          /* basic orthogonal pattern */
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

#define MAXWIDTH (14)

#define ROWBITS ((1<<width)-1)
#define BASEBITS (params[P_BASEBITS])
#define BASEFACTOR (1<<BASEBITS)
#define MAXOFFSET ((((row) -1) >> width) - 1)

#define ROW(i) (rows[i] & ROWBITS)
#define ROFFSET(i) (rows[i] >> width)
#define EMPTY(i) (rows[i] == (row)-1)
#define MAKEEMPTY(i) rows[i] = (row)-1
#define PARENT(i) (base[(i)>>BASEBITS]+ROFFSET(i))
#define FIRSTBASE(i) (((i) & ((1<<BASEBITS) - 1)) == 0)

#define MINDEEP ((params[P_MINDEEP]>0) ? params[P_MINDEEP] : period)

int gcd(int a, int b) {
   if (a > b) return gcd(b,a);
   else if (a == 0) return b;
   else return gcd(b-a,a);
}

/* =========================================== */
/*  Lookup tables to determine successor rows  */
/* =========================================== */

const char *rulekeys[] = {
   "0", "1c", "1e", "2a", "1c", "2c", "2a", "3i",
   "1e", "2k", "2e", "3j", "2a", "3n", "3a", "4a",
   "1c", "2n", "2k", "3q", "2c", "3c", "3n", "4n",
   "2a", "3q", "3j", "4w", "3i", "4n", "4a", "5a",
   "1e", "2k", "2i", "3r", "2k", "3y", "3r", "4t",
   "2e", "3k", "3e", "4j", "3j", "4k", "4r", "5n",
   "2a", "3q", "3r", "4z", "3n", "4y", "4i", "5r",
   "3a", "4q", "4r", "5q", "4a", "5j", "5i", "6a",
   "1c", "2c", "2k", "3n", "2n", "3c", "3q", "4n",
   "2k", "3y", "3k", "4k", "3q", "4y", "4q", "5j",
   "2c", "3c", "3y", "4y", "3c", "4c", "4y", "5e",
   "3n", "4y", "4k", "5k", "4n", "5e", "5j", "6e",
   "2a", "3n", "3r", "4i", "3q", "4y", "4z", "5r",
   "3j", "4k", "4j", "5y", "4w", "5k", "5q", "6k",
   "3i", "4n", "4t", "5r", "4n", "5e", "5r", "6i",
   "4a", "5j", "5n", "6k", "5a", "6e", "6a", "7e",
   "1e", "2a", "2e", "3a", "2k", "3n", "3j", "4a",
   "2i", "3r", "3e", "4r", "3r", "4i", "4r", "5i",
   "2k", "3q", "3k", "4q", "3y", "4y", "4k", "5j",
   "3r", "4z", "4j", "5q", "4t", "5r", "5n", "6a",
   "2e", "3j", "3e", "4r", "3k", "4k", "4j", "5n",
   "3e", "4j", "4e", "5c", "4j", "5y", "5c", "6c",
   "3j", "4w", "4j", "5q", "4k", "5k", "5y", "6k",
   "4r", "5q", "5c", "6n", "5n", "6k", "6c", "7c",
   "2a", "3i", "3j", "4a", "3q", "4n", "4w", "5a",
   "3r", "4t", "4j", "5n", "4z", "5r", "5q", "6a",
   "3n", "4n", "4k", "5j", "4y", "5e", "5k", "6e",
   "4i", "5r", "5y", "6k", "5r", "6i", "6k", "7e",
   "3a", "4a", "4r", "5i", "4q", "5j", "5q", "6a",
   "4r", "5n", "5c", "6c", "5q", "6k", "6n", "7c",
   "4a", "5a", "5n", "6a", "5j", "6e", "6k", "7e",
   "5i", "6a", "6c", "7c", "6a", "7e", "7c", "8" } ;
/*
 *   Parses the rule.  If there's an error, return a string describing the
 *   error.  Fills in the 512-element array pointed to by tab.
 */
const char *parseRule(const char *rule, int *tab) {
   const char *p = rule ;
   for (int i=0; i<512; i++)
      tab[i] = 0 ;
   for (int bs=0; bs<512; bs += 256) {
      if (bs == 0) {
         if (*p != 'B' && *p != 'b')
            return "Expected B at start of rule" ;
      } else {
         if (*p != 'S' && *p != 's')
            return "Expected S after slash" ;
      }
      p++ ;
      while (*p != '/' && *p != 0) {
         if (!('0' <= *p || *p <= '9'))
            return "Missing number in rule" ;
         char dig = *p++ ;
         int neg = 0 ;
         if (*p == '/' || *p == 0 || *p == '-' || ('0' <= *p && *p <= '8'))
            for (int i=0; i<256; i++)
               if (rulekeys[i][0] == dig)
                  tab[bs+i] = 1 ;
         for (; *p != '/' && *p != 0 && !('0' <= *p && *p <= '8'); p++) {
            if (*p == '-') {
               if (neg)
                  return "Can't have multiple negation signs" ;
               neg = 1 ;
            } else if ('a' <= *p && *p <= 'z') {
               int used = 0 ;
               for (int i=0; i<256; i++)
                  if (rulekeys[i][0] == dig && rulekeys[i][1] == *p) {
                     tab[bs+i] = 1-neg ;
                     used++ ;
                  }
               if (!used)
                  return "Unexpected character in rule" ;
            } else
               return "Unexpected character in rule" ;
         }
      }
      if (bs == 0) {
         if (*p++ != '/')
            return "Missing expected slash between b and s" ;
      } else {
         if (*p++ != 0)
            return "Extra unparsed junk at end of rule string" ;
      }
   }
   return 0 ;
}

#ifndef QSIMPLE
void makePhases();
#endif

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
char * patternBuf;

void bufRLE(char c) {
  if (RLEcount > 0 && c != RLEchar) {
      if (RLElineWidth++ >= MAXRLELINEWIDTH) {
         if (RLEchar != '\n') sprintf(patternBuf+strlen(patternBuf),"\n");
         RLElineWidth = 0;
      }
    if (RLEcount == 1) strncat(patternBuf,&RLEchar,1);
    else {
       sprintf(patternBuf+strlen(patternBuf),"%d%c", RLEcount, RLEchar);
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

void bufRow(unsigned long rr, unsigned long r, int shift) {
   while (r | rr) {
      if (shift == 0)
         bufRLE(r & 1 ? 'o' : 'b');
      else shift--;
      r >>= 1;
      if (rr & 1) r |= (1<<31);
      rr >>= 1;
   }
   bufRLE('$');
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

int sxsAllocRows =  0;
unsigned long * sxsAllocData;
unsigned long * sxsAllocData2;
int oldnrows = 0;
unsigned long * oldsrows;
unsigned long * oldssrows;

/* Buffers RLE into patternBuf; returns 1 if pattern successfully buffered */ 
int bufferPattern(node b, row *pRows, int nodeRow, uint32_t lastRow, int printExpected){
   node c;
   int nrows = 0;
   int skewAmount = 0;
   int swidth;
   int sxsNeeded;
   int p, i, j, margin;
   unsigned long *srows, *ssrows;

   uint32_t currRow = lastRow;
   int nDeepRows = 0;
   int nodeDiff;

   if(pRows != NULL){
      while(pRows[currRow] == 0){
         if(currRow == 0){
            if(!printExpected) return 0;
            printf("Success called on search root!\n");
            aborting = 1;
            return 0;
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
            if(!printExpected) return 0;
            printf("Success called on search root!\n");
            aborting = 1;
            return 0;
         }
      }
   }
   if(nrows < 0) nrows = 0;
   
   for (p = 0; p < period-1; p++) b = PARENT(b);
   if (b == 0) {
      if(!printExpected) return 0;
      printf("Success called on search root!\n");
      aborting = 1;
      return 0;
   }
   
   /* count rows */
   c = b;
   while (c != 0) {
      for (p = 0; p < period; p++)
         c = PARENT(c);
      nrows++;
   }
   
   /* build data structure of rows so we can reduce width etc */
   sxsNeeded = nrows+MAXWIDTH+1;
   if (!sxsAllocRows)
   {
      sxsAllocRows = sxsNeeded;
      sxsAllocData = (unsigned long*)malloc(sxsAllocRows * sizeof(unsigned long));
      sxsAllocData2 = (unsigned long*)malloc(sxsAllocRows * sizeof(unsigned long));
      oldsrows = (unsigned long*)calloc(sxsAllocRows, sizeof(unsigned long));
      oldssrows = (unsigned long*)calloc(sxsAllocRows, sizeof(unsigned long));
      patternBuf = (char*)malloc(((2 * MAXWIDTH + 4) * sxsAllocRows + 300) * sizeof(char));
   }
   else if (sxsAllocRows < sxsNeeded)
   {
      sxsAllocRows = sxsNeeded;
      sxsAllocData = (unsigned long*)realloc(sxsAllocData, sxsAllocRows * sizeof(unsigned long));
      sxsAllocData2 = (unsigned long*)realloc(sxsAllocData2, sxsAllocRows * sizeof(unsigned long));
   }
   srows  = sxsAllocData;
   ssrows = sxsAllocData2;
   
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
            return 0;
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
   if(allEmpty) return 0;
   
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
   if(printExpected){
      if (nrows == oldnrows) {
         int different = 0;
         for (i = 0; i < nrows && !different; i++)
            different = (srows[i] != oldsrows[i] || ssrows[i] != oldssrows[i]);
         if (!different) return 0;
      }
      
      /* replace previous saved rows with new rows */
      oldnrows = nrows;
      oldsrows = (unsigned long*)realloc(oldsrows, sxsAllocRows * sizeof(unsigned long));
      oldssrows = (unsigned long*)realloc(oldssrows, sxsAllocRows * sizeof(unsigned long));
      memcpy(oldsrows, sxsAllocData, sxsAllocRows * sizeof(unsigned long));
      memcpy(oldssrows, sxsAllocData2, sxsAllocRows * sizeof(unsigned long));
   }
   
   /* Buffer output */
   patternBuf = (char*)realloc(patternBuf, ((2 * MAXWIDTH + 4) * sxsAllocRows + 300) * sizeof(char));
   
   sprintf(patternBuf,"x = %d, y = %d, rule = %s\n", swidth - margin, nrows, rule);
      
   while (nrows-- > 0) {
      if (margin > nrows) bufRow(ssrows[nrows], srows[nrows], margin - nrows);
      else bufRow(ssrows[nrows], srows[nrows], 0);
   }
   RLEchar = '!';
   bufRLE('\0');
   sprintf(patternBuf+strlen(patternBuf),"\n");
   
   if(printExpected){
      numFound++;
      if(params[P_NUMSHIPS] > 0){
         if(--params[P_NUMSHIPS] == 0) aborting = 3;  /* use 3 to flag that we reached ship limit */
      }
   }
   
   return 1;
}

void success(node b, row *pRows, int nodeRow, uint32_t lastRow){
   if(bufferPattern(b, pRows, nodeRow, lastRow, 1))
      printf("\n%s\n",patternBuf);
   fflush(stdout);
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

/* Maintain phase of queue nodes.  After dequeue(), the global variable phase
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

/* =================== */
/*  Dump search state  */
/* =================== */

int dumpNum = 1;
char dumpFile[256];
const char *dumpRoot = "dump";
char loadDumpRoot[251];    /* used for loading dump root from file */
int dumpFlag = 0;    /* Dump status flags, possible values follow */
#define DUMPPENDING (1)
#define DUMPFAILURE (2)
#define DUMPSUCCESS (3)

FILE * openDumpFile()
{
   FILE * fp;

   while (dumpNum < 100000)
   {
      sprintf(dumpFile,"%s%05d",dumpRoot,dumpNum++);
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
   fprintf(fp,"%s\n",dumpRoot);
   for (i = 0; i < NUM_PARAMS; i++)
      fprintf(fp,"%d\n",params[i]);
   fprintf(fp,"%d\n",width);
   fprintf(fp,"%d\n",period);
   fprintf(fp,"%d\n",offset);
   fprintf(fp,"%d\n",mode);

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
** converts parent bits to back parent pointers.  The search
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
                     so, if x doesn't point to y, y must be unused and can be removed. */
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
      if (y-- == 0) break;    /* circumlocution for while (y >= 0) because x is unsigned */
   }
   qStart = ++x;     /* mark start of queue */
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
      if (ROFFSET(x)) {   /* skip forward to next parent */
         y++;
         while (EMPTY(y)) y++;
      }
      enqueue(y,ROW(x));
      //if (aborting) return;    /* why is this here? value of aborting is not changed by enqueue(). */
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
   /* First loop of part 1 requires qTail-1 to be non-empty.  Make it so */
   while(EMPTY(qTail-1))
      qTail--;
   
   doCompactPart1();
   if (dumpFlag == DUMPPENDING) dumpState();
   doCompactPart2();
}

/* ================= */
/*  Lookahead cache  */
/* ================= */

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

/* ========================== */
/*  Primary search functions  */
/* ========================== */

void process(node theNode);
int depthFirst(node theNode, long howDeep, uint16_t **pInd, int *pRemain, row *pRows);

static void deepen(){
   node i;

   /* compute amount to deepen, apply reduction if too deep */
#ifdef PATCH07
   timeStamp();
#endif
   printf("Queue full");
   i = currentDepth();
   if (i >= params[P_LASTDEEP]) deepeningAmount = MINDEEP;
   else deepeningAmount = params[P_LASTDEEP] + MINDEEP - i;   /* go at least MINDEEP deeper */

   params[P_LASTDEEP] = i + deepeningAmount;

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

/* ========================== */
/*  Print usage instructions  */
/* ========================== */

/* Note: currently reserving -v for potentially editing an array of extra variables */
void usage(){
#ifndef QSIMPLE
   printf("Usage: \"qfind options\"\n");
   printf("  e.g. \"qfind -r B3/S23 -p 3 -y 1 -w 6 -s even\" searches Life (rule B3/S23)\n");
   printf("  for c/3 orthogonal spaceships with even bilateral symmetry and a\n");
   printf("  logical width of 6 (full width 12).\n");
#else
   printf("Three required parameters, the period, offset, and width, must be\n");
   printf("set within the code before it is compiled. You have compiled with\n");
   printf("\n");
   printf("Period: %d\n",PERIOD);
   printf("Offset: %d\n",OFFSET);
   printf("Width:  %d\n",WIDTH);
#endif
   printf("\n");
   printf("Available options:\n");
   printf("  -r bNN/sNN  searches for spaceships in the specified rule (default: b3/s23)\n");
   printf("              Non-totalistic rules can be entered using Hensel notation.\n");
   printf("\n");
#ifndef QSIMPLE
   printf("  -p NN  searches for spaceships with period NN\n");
   printf("  -y NN  searches for spaceships that travel NN cells every period\n");
   printf("  -w NN  searches for spaceships with logical width NN\n");
   printf("         (full width depends on symmetry type)\n");
#endif
   printf("  -s FF  searches for spaceships with symmetry type FF\n");
   printf("         Valid symmetry types are asymmetric, odd, even, and gutter.\n");
   printf("\n");
   printf("  -t NN  runs search using NN threads during deepening step (default: 1)\n");
   printf("  -i NN  sets minimum deepening increment to NN (default: period)\n");
   printf("  -n NN  deepens to total depth at least NN during first deepening step\n");
   printf("         (total depth includes depth of BFS queue)\n");
   printf("  -f NN  stops search if NN ships are found (default: no limit)\n");
   printf("  -q NN  sets the BFS queue size to 2^NN (default: %d)\n",QBITS);
   printf("  -h NN  sets the hash table size to 2^NN (default: %d)\n",HASHBITS);
   printf("         Use -h 0 to disable duplicate elimination.\n");

   printf("  -b NN  groups 2^NN queue entries to an index node (default: 4)\n");
   printf("  -m NN  limits memory usage to NN megabytes (default: no limit)\n");
#ifndef NOCACHE
   printf("  -c NN  allocates NN megabytes per thread for lookahead cache (default: 32)\n");
#endif
   printf("  -z     disables output during deepening step\n");
   printf("         (useful for searches that find many spaceships)\n");
   printf("  -a     Suppresses output of longest partial result at end of search\n");
   printf("\n");
   printf("  -e FF  uses rows in the file FF as the initial rows for the search\n");
   printf("         (use the companion Golly python script to easily generate the\n");
   printf("         initial row file)\n");
   printf("  -d FF  dumps the search state after each queue compaction using\n");
   printf("         file name prefix FF\n");
   printf("  -l FF  loads the search state from the file FF\n");
   printf("  -j NN  splits the search state into at most NN files\n");
   printf("         (uses the file name prefix defined by -d option\n");
   printf("  -u     previews partial results from the loaded state\n");
   printf("\n");
   printf("  -o     uses naive search order (not recommended)\n");
   printf("\n");
   printf("  --help  prints usage instructions and exits\n");
}

/* ======================== */
/*  Echo loaded parameters  */
/* ======================== */

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
   if(!params[P_REORDER]) printf("Use naive search order\n");
   printf("Queue size: 2^%d\n",params[P_QBITS]);
   printf("Hash table size: 2^%d\n",params[P_HASHBITS]);
   printf("Minimum deepening increment: %d\n",MINDEEP);
   if(params[P_PRINTDEEP] == 0)printf("Output disabled while deepening\n");
#ifndef NOCACHE
   printf("Cache memory per thread: %d megabytes\n", params[P_CACHEMEM]);
#endif
   if(params[P_MEMLIMIT] >= 0) printf("Memory limit: %d megabytes\n",params[P_MEMLIMIT]);
   printf("Number of threads: %d\n",params[P_NUMTHREADS]);
   if(params[P_LONGEST] == 0)printf("Printing of longest partial result disabled\n");
}

/* ========================= */
/*  Preview partial results  */
/* ========================= */

static void preview(int allPhases) {
   node i,j,k;
   int ph;

   for (i = qHead; (i<qTail) && EMPTY(i); i++);
   for (j = qTail-1; (j>i) && EMPTY(j); j--);
   if (j<i) return;
   
   while (j>=i && !aborting) {
      if (!EMPTY(j)) {
         success(j, NULL, 0, 0);
         if (allPhases == 0) {
            k=j;
            for (ph = 1; ph < period; ph++) {
               k=PARENT(k);
               success(k, NULL, 0, 0);
            }
         }
      }
      j--;
   }
}

/* ============================================================ */
/*  Check parameters for validity and exit if there are errors  */
/* ============================================================ */

int splitNum = 0;
int loadDumpFlag = 0;
int previewFlag = 0;
int initRowsFlag = 0;
int newLastDeep = 0;

void checkParams(){
   int exitFlag = 0;
   const char *ruleError;
   
   /* Errors */
   ruleError = parseRule(rule, nttable);
   if (ruleError != 0){
      fprintf(stderr, "Error: failed to parse rule %s\n", rule);
      fprintf(stderr, "       %s\n", ruleError);
      exitFlag = 1;
   }
#ifdef QSIMPLE
   if(gcd(PERIOD,OFFSET) > 1){
      fprintf(stderr, "Error: qfind-s does not support gcd(PERIOD,OFFSET) > 1. Use qfind instead.\n");
      exitFlag = 1;
   }
#else
   if(params[P_WIDTH] < 1 || params[P_PERIOD] < 1 || params[P_OFFSET] < 1){
      fprintf(stderr, "Error: period (-p), translation (-y), and width (-w) must be positive integers.\n");
      exitFlag = 1;
   }
   if(params[P_PERIOD] > MAXPERIOD){
      fprintf(stderr, "Error: maximum allowed period (%d) exceeded.\n", MAXPERIOD);
      exitFlag = 1;
   }
   if(params[P_OFFSET] > params[P_PERIOD] && params[P_PERIOD] > 0){
      fprintf(stderr, "Error: translation (-y) cannot exceed period (-p).\n");
      exitFlag = 1;
   }
   if(params[P_OFFSET] == params[P_PERIOD] && params[P_PERIOD] > 0){
      fprintf(stderr, "Error: photons are not supported.\n");
      exitFlag = 1;
   }
#endif
   if(params[P_SYMMETRY] == 0){
      fprintf(stderr, "Error: you must specify a symmetry type (-s).\n");
      exitFlag = 1;
   }
   if(previewFlag && !loadDumpFlag){
      fprintf(stderr, "Error: the search state must be loaded from a file to preview partial results.\n");
      exitFlag = 1;
   }
   if(initRowsFlag && loadDumpFlag){
      fprintf(stderr, "Error: Initial rows file cannot be used when the search state is loaded from a\n       saved state.\n");
      exitFlag = 1;
   }
   
   /* Warnings */
   if(2 * params[P_OFFSET] > params[P_PERIOD] && params[P_PERIOD] > 0){
      fprintf(stderr, "Warning: searches for speeds exceeding c/2 may not work correctly.\n");
   }
#ifdef NOCACHE
   if(5 * params[P_OFFSET] > params[P_PERIOD] && params[P_PERIOD] > 0){
      fprintf(stderr, "Warning: Searches for speeds exceeding c/5 may be slower without caching.\n         It is recommended that you recompile with NOCACHE undefined.\n");
   }
#else
   if(5 * params[P_OFFSET] <= params[P_PERIOD] && params[P_OFFSET] > 0){
      fprintf(stderr, "Warning: Searches for speeds at or below c/5 may be slower with caching.\n         It is recommended that you recompile with NOCACHE defined.\n");
   }
#endif
   
   /* exit if there are errors */
   if(exitFlag){
      fprintf(stderr, "\nUse --help for a list of available options.\n");
      exit(1);
   }
   fprintf(stderr, "\n");
}

/* ============================ */
/*  Load saved state from file  */
/* ============================ */

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

void loadParams() {
   FILE * fp;
   int i;
   
   /* reset flag to prevent modification of params[P_LASTDEEP] at start of search */
   newLastDeep = 0;
   
   fp = fopen(loadFile, "r");
   if (!fp) loadFail();
   if (loadUInt(fp) != FILEVERSION)
   {
      printf("Incompatible file version\n");
      exit(1);
   }
   
   /* Load rule */
   if (fscanf(fp,"%255s\n",loadRule) != 1) loadFail();
   rule = loadRule;
   
   /* Load dump root */
   if (fscanf(fp,"%250s\n",loadDumpRoot) != 1) loadFail();
   dumpRoot = loadDumpRoot;
   
   /* Load parameters */
   for (i = 0; i < NUM_PARAMS; i++)
      params[i] = loadInt(fp);
}

loadState(){
   FILE * fp;
   int i;
   
   fp = fopen(loadFile, "r");
   if (!fp) loadFail();
   
   loadUInt(fp);                                   /* skip file version */
   fscanf(fp, "%*[^\n]\n");                        /* skip rule */
   fscanf(fp, "%*[^\n]\n");                        /* skip dump root */
   for (i = 0; i < NUM_PARAMS; i++) loadInt(fp);   /* skip parameters */
   
   /* Load / initialise globals */
   width          = loadInt(fp);
   period         = loadInt(fp);
   offset         = loadInt(fp);
   mode           = (Mode)(loadInt(fp));
   
   deepeningAmount = period; /* Currently redundant, since it's recalculated */
   aborting        = 0;
   nRowsInState    = period+period;   /* how many rows needed to compute successor graph? */
   
   params[P_DEPTHLIMIT] = DEFAULT_DEPTHLIMIT;

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
   if (qStart > QSIZE || qStart < QSIZE/16) {
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
   
   /* Let the user know that we got this far (suppress if splitting) */
   if(!splitNum) printf("State successfully loaded from file %s\n",loadFile);
   
   fflush(stdout);
}

/* ================================================= */
/*  Load initial rows for extending partial results  */
/* ================================================= */

void printRow(row theRow){
   int i;
   for(i = width - 1; i >= 0; i--) printf("%c",(theRow & 1 << i ? 'o' : '.'));
   printf("\n");
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

/* ================================= */
/*  Parse options and set up search  */
/* ================================= */

void setDefaultParams(){
#ifdef QSIMPLE
   params[P_WIDTH] = WIDTH;
   params[P_PERIOD] = PERIOD;
   params[P_OFFSET] = OFFSET;
#else
   params[P_WIDTH] = 0;
   params[P_PERIOD] = 0;
   params[P_OFFSET] = 0;
#endif
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
   params[P_LASTDEEP] = 0;
   params[P_NUMSHIPS] = 0;
}

/* Note: currently reserving -v for potentially editing an array of extra variables */
void parseOptions(int argc, char *argv[]){
   while(--argc > 0){               /* read input parameters */
      if ((*++argv)[0] == '-'){
         switch ((*argv)[1]){
            case 'r': case 'R':
               --argc;
               rule = *++argv;
               break;
#ifndef QSIMPLE
            case 'p': case 'P':
               --argc;
               sscanf(*++argv, "%d", &params[P_PERIOD]);
               break;
            case 'y': case 'Y':
               --argc;
               sscanf(*++argv, "%d", &params[P_OFFSET]);
              break;
            case 'w': case 'W':
               --argc;
               sscanf(*++argv, "%d", &params[P_WIDTH]);
               break;
#endif
            case 's': case 'S':
               --argc;
               switch((*++argv)[0]) {
                  case 'a': case 'A':
                     params[P_SYMMETRY] = SYM_ASYM; mode = asymmetric; break;
                  case 'o': case 'O':
                     params[P_SYMMETRY] = SYM_ODD; mode = odd; break;
                  case 'e': case 'E':
                     params[P_SYMMETRY] = SYM_EVEN; mode = even; break;
                  case 'g': case 'G':
                     params[P_SYMMETRY] = SYM_GUTTER; mode = gutter; break;
                  default:
                     fprintf(stderr, "Error: unrecognized symmetry type %s\n", *argv);
                     fprintf(stderr, "\nUse --help for a list of available options.\n");
                     exit(1);
                     break;
               }
               break;
            case 'm': case 'M':
               --argc;
               sscanf(*++argv, "%d", &params[P_MEMLIMIT]);
               break;
            case 'n': case 'N':
               --argc;
               sscanf(*++argv, "%d", &params[P_LASTDEEP]);
               newLastDeep = 1;
               break;
            case 'c': case 'C':
               --argc;
               sscanf(*++argv, "%d", &params[P_CACHEMEM]);
               break;
            case 'i': case 'I':
               --argc;
               sscanf(*++argv, "%d", &params[P_MINDEEP]);
               break;
            case 'q': case 'Q':
               --argc;
               sscanf(*++argv, "%d", &params[P_QBITS]);
               break;
            case 'h': case 'H':
               --argc;
               sscanf(*++argv, "%d", &params[P_HASHBITS]);
               break;
            case 'b': case 'B':
               --argc;
               sscanf(*++argv, "%d", &params[P_BASEBITS]);
               break;
            case 't': case 'T':
               --argc;
               sscanf(*++argv, "%d", &params[P_NUMTHREADS]);
               break;
            case 'f': case 'F':
               --argc;
               sscanf(*++argv, "%d", &params[P_NUMSHIPS]);
               break;
            case 'z': case 'Z':
               params[P_PRINTDEEP] = 0;
               break;
            case 'a': case 'A':
               params[P_LONGEST] = 0;
               break;
            case 'o': case 'O':
               params[P_REORDER] = 0;
               break;
            case 'u': case 'U':
               previewFlag = 1;
               break;
            case 'd': case 'D':
               --argc;
               dumpRoot = *++argv;
               params[P_CHECKPOINT] = 1;
               break;
            case 'j': case 'J':
               --argc;
               sscanf(*++argv, "%d", &splitNum);
               if(splitNum < 0) splitNum = 0;
               break;
            case 'e': case 'E':
               --argc;
               initRows = *++argv;
               initRowsFlag = 1;
               break;
            case 'l': case 'L':
               --argc;
               loadFile = *++argv;
               loadDumpFlag = 1;
               loadParams();
               break;
            case '-':
               if(!strcmp(*argv,"--help") || !strcmp(*argv,"--Help")){
                  usage();
                  exit(0);
               }
               else{
                  fprintf(stderr, "Error: unrecognized option %s\n", *argv);
                  fprintf(stderr, "\nUse --help for a list of available options.\n");
                  exit(1);
               }
               break;
           default:
              fprintf(stderr, "Error: unrecognized option %s\n", *argv);
              fprintf(stderr, "\nUse --help for a list of available options.\n");
              exit(1);
              break;
         }
      }
   }
}

void searchSetup(){
   checkParams();  /* Exit if parameters are invalid */
   
   if(loadDumpFlag) loadState();
   else {
      width = params[P_WIDTH];
      period = params[P_PERIOD];
      offset = params[P_OFFSET];
      deepeningAmount = period;
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
      
      resetQ();
      resetHash();
      
      enqueue(0,0);
      
      if(initRowsFlag) loadInitRows(initRows);
   }
   
   if(previewFlag){
      preview(1);
      exit(0);
   }
   
   /* correction of params[P_LASTDEEP] after modification */
   if(newLastDeep){
      params[P_LASTDEEP] -= MINDEEP;
      if(params[P_LASTDEEP] < 0) params[P_LASTDEEP] = 0;
   }
   
   /* split queue across multiple files */
   if(splitNum > 0){
      node x;
      int firstDumpNum = 0;
      int totalNodes = 0;
      
      echoParams();
      printf("\n");
      
      if(!loadDumpFlag || qHead == 0 || splitNum == 1){
         dumpFlag = DUMPPENDING;
         if(qHead == 0){      /* can't use doCompact() here, because it tries to access rows[-1] */
            qStart = qHead;
            qEnd = qTail;
            dumpState();
         }
         else doCompact();
         if(dumpFlag == DUMPSUCCESS){
            printf("State dumped to %s\n",dumpFile);
            exit(0);
         }
         else{
            fprintf(stderr, "Error: dump failed.\n");
            exit(1);
         }
      }
      
      if(splitNum >= 100000){
         fprintf(stderr, "Warning: queue cannot be split into more than 99999 files.\n");
         splitNum = 99999;
      }
      
      /* count nodes in queue */
      for(x = qHead; x < qTail; x++){
         if(!EMPTY(x)) totalNodes++;
      }
      
      /* nodes per file is rounded up */
      int nodesPerFile = (totalNodes - 1) / splitNum + 1;
      
      printf("Splitting search state with %d queue nodes per file\n",nodesPerFile);
      
      /* save qHead and qTail, as creating the pieces will change their values */
      node fixedQHead = qHead;
      node fixedQTail = qTail;
      
      /* delete the queue; we will reload it as needed */
      free(base);
      free(rows);
      free(hash);
      
      node currNode = fixedQHead;
      
      while(currNode < fixedQTail){
         
         /* load the queue */
         loadState();
         
         /* empty everything before currNode */
         for(x = fixedQHead; x < currNode; x++) MAKEEMPTY(x);
         
         /* skip the specified number of nonempty nodes */
         int i = 0;
         while(i < nodesPerFile && x < fixedQTail){
            if(!EMPTY(x))
               i++;
            x++;
         }
         
         /* update currNode */
         currNode = x;
         
         /* empty everything after specified nonempty nodes */
         while(x < fixedQTail){
            MAKEEMPTY(x);
            x++;
         }
         
         /* save the piece */
         dumpFlag = DUMPPENDING;
         doCompact();
         
         if(!firstDumpNum) firstDumpNum = dumpNum - 1;
         
         if (dumpFlag != DUMPSUCCESS){
            printf("Failed to save %s\n",dumpFile);
            exit(1);
         }
         
         if(dumpNum >= 100000){
            fprintf(stderr, "Error: dump file number limit reached.\n");
            fprintf(stderr, "       Try splitting the queue in a new directory.\n");
            exit(1);
         }
         
         /* prevent memory leak */
         free(base);
         free(rows);
         free(hash);
      }
      
      printf("Saved pieces in files %s%05d to %s\n",dumpRoot,firstDumpNum,dumpFile);
      exit(0);
   }
   
   omp_set_num_threads(params[P_NUMTHREADS]);
   
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
   
#ifndef QSIMPLE
   makePhases();
#endif
   fasterTable();
   makeTables();
   
   rephase();
}

void finalReport(){
   printf("Search complete.\n\n");
   
   printf("%d spaceship%s found.\n",numFound,(numFound == 1) ? "" : "s");
   printf("Maximum depth reached: %d\n",longest);
   if(params[P_LONGEST] && aborting != 3){ /* aborting == 3 means we reached ship limit */
      if(patternBuf) printf("Longest partial result:\n\n%s",patternBuf);
      else printf("No partial results found.\n");
   }
}