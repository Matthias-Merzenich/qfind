/* This header file contains functions common to both qfind and qfind-s.
** Some functions contained in this file work slightly differently depending
** on whether qfind or qfind-s is being compiled.  Such differences are
** determined by the presence of the macro QSIMPLE defined in qfind-s.c.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <stdatomic.h>

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
   #define omp_set_num_threads(x)
#endif

#ifdef QSIMPLE
   #define WHICHPROGRAM qfind-simple
#else
   #define WHICHPROGRAM qfind
#endif

#define STR(x) #x
#define XSTR(x) STR(x)
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define BANNER XSTR(WHICHPROGRAM)" v2.4b by Matthias Merzenich, 3 September 2025"

#define FILEVERSION ((unsigned long) 2025090301)  /* yyyymmddnn */

#define MAXPERIOD 30
#define MAXDUMPROOT 50     /* maximum allowed length of dump root */
#define DUMPLIMIT 100000   /* maximum allowed number of sequential dumps */
#define CHUNK_SIZE 1
#define QBITS 20
#define HASHBITS 20
#define DEFAULT_DEPTHLIMIT (qBits-3)
#define DEFAULT_CACHEMEM 32

#define P_WIDTH 0
#define P_PERIOD 1
#define P_OFFSET 2
#define P_SYMMETRY 3
#define P_REORDER 4           /* currently cannot be set at runtime */
#define P_DUMPMODE 5
#define P_BASEBITS 6
#define P_QBITS 7
#define P_HASHBITS 8
#define P_DEPTHLIMIT 9        /* currently cannot be set at runtime */
#define P_NUMTHREADS 10
#define P_MINDEEP 11
#define P_MEMLIMIT 12
#define P_CACHEMEM 13
#define P_PRINTDEEP 14
#define P_LONGEST 15
#define P_FIRSTDEEP 16
#define P_NUMSHIPS 17
#define P_MINEXTENSION 18
#define P_FULLPERIOD 19
#define P_BOUNDARYSYM 20
#define P_DUMPINTERVAL 21
#define P_EVERYDEPTH 22
#define P_EARLYEXIT 23

#define NUM_PARAMS 24U

#define SYM_UNDEF 0
#define SYM_ASYM 1
#define SYM_ODD 2
#define SYM_EVEN 3
#define SYM_GUTTER 4

const char *rule = "B3/S23";     /* Default rule set to B3/S23 (Life) */
char loadRule[151];              /* Used for loading rule from file */
char baseRule[151];              /* Used in case of forbidden conditions */

char *initRows;

int params[NUM_PARAMS];
int width;
int nRowsInState;                /* Could be replaced with 2*period.  Should be removed. */
int phase;

int period;
int offset;
int lastDeep = 0;
int numFound = 0;       /* number of spaceships found so far */
int longest = 0;        /* length of current longest partial result */

int aborting = 0;       /* A flag to indicate that we are stopping the search.
                        ** valid aborting values:
                        **   0: not aborting
                        **   1: fatal error
                        **   2: queue size limit reached
                        **   3: desired number of ships found
                        */

int gutterSkew = 0;     /* number of cells to skew halves in gutter symmetric search */

/* the big data structures */
#define qBits params[P_QBITS]
#define QSIZE (1LU<<qBits)

#define hashBits params[P_HASHBITS]
#define HASHSIZE (1LU<<hashBits)
#define HASHMASK (HASHSIZE - 1)

typedef uint32_t node;
typedef uint16_t row;

row * rows;
node * base;
node * hash;

int nttable[512];
uint16_t **gInd3;
uint32_t *gcount;
uint16_t *gRows;
long long memusage = 0;
long long memlimit = 0;

row **deepRows = 0;
uint32_t *deepRowIndices;
uint32_t deepQHead, deepQTail, oldDeepQHead;

#ifndef NOCACHE
long long cachesize;
typedef struct {
   uint16_t *p1, *p2, *p3;
   int abn, r;
} cacheentry;

cacheentry *totalCache;

cacheentry **cache;
#endif

/* Representation of vertices.
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

#define MAXWIDTH (14)   /* hard limit as long as rows are of type uint16_t */

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

#define MINDEEP ((params[P_MINDEEP]>0) ? params[P_MINDEEP] : 3)

void optError(const char *errorMsg, const char *opt);
void printError(const char *errorMsg);

int gcd(int a, int b) {
   if (a > b) return gcd(b,a);
   else if (a == 0) return b;
   else return gcd(b-a,a);
}

/* =============================== */
/*  Display current date and time  */
/* =============================== */

static char timeStr[20] = "00/00/00 00:00:00 ";

static void timeStamp() {
   time_t t;
   
   time(&t);
   strftime(timeStr, 20, "%d/%m/%y %H:%M:%S ", localtime(&t));
   printf("%s", timeStr);
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
   "5i", "6a", "6c", "7c", "6a", "7e", "7c", "8" };

/*   Parses the rule.  If there's an error, return a string describing the error.
**   Fills in the 512-element array pointed to by tab.
*/
const char *parseRule(const char *rule, int *tab) {
   const char *p = rule;
   int tempTab[256];    /* needed to keep track of values when setting forbidden conditions */
   for (int i=0; i<512; i++)
      tab[i] = 0;
   for (int bs=0; bs<512; bs += 256) {
      if (bs == 0) {
         if (*p != 'B' && *p != 'b')
            return "Expected B at start of rule";
      } else {
         if (*p != 'S' && *p != 's')
            return "Expected S after slash";
      }
      p++;
      int allowed = 1;
      while (*p != '/' && *p != '\0') {
         if (*p == '~'){
            p++;
            if (allowed == -1 || *p == '~'){
               if (bs)
                  return "Can't have multiple tildes in survival conditions";
               else
                  return "Can't have multiple tildes in birth conditions";
            }
            if (*p == '/' || *p == '\0')
               continue;
            allowed = -1;  /* table entry will be -1 if condition is forbidden */
         }
         if (!('0' <= *p && *p <= '9'))
            return "Missing number in rule";
         if (*p == '9')
            return "Unexpected character in rule";
         char dig = *p++;
         int neg = 0;
         if (*p == '/' || *p == '\0' || (*p == '-' && allowed == 1) || *p == '~' || ('0' <= *p && *p <= '8'))
            for (int i=0; i<256; i++)
               if (rulekeys[i][0] == dig)
                  tab[bs+i] = 1*allowed;
         int forbiddenCount = 0;
         if (*p == '-'){
            neg = 1;
            for (int i=0; i<256; i++)
               tempTab[i] = 0;
            p++;
         }
         for (; *p != '/' && *p != '\0' && *p != '~' && !('0' <= *p && *p <= '8'); p++) {
            if (*p == '-')
               return "Improperly placed negation sign";
            if ('a' <= *p && *p <= 'z') {
               int used = 0;
               for (int i=0; i<256; i++){
                  if (rulekeys[i][0] == dig){
                     if (rulekeys[i][1] == *p){
                        if (allowed == 1)
                           tab[bs+i] = (1-neg);
                        else if (!neg)
                           tab[bs+i] = -1;
                        used++;
                     }
                     else if (neg && allowed == -1)
                        tempTab[i]++;
                  }
               }
               if (neg && allowed == -1)
                  forbiddenCount++;
               if (!used)
                  return "Unexpected character in rule";
            }
            else
               return "Unexpected character in rule";
         }
         if (neg && allowed == -1)
            for (int i=0; i<256; i++)
               if (tempTab[i] == forbiddenCount)
                  tab[bs+i] = -1;
      }
      if (bs == 0) {
         if (*p++ != '/')
            return "Missing expected slash between B and S";
      } else {
         if (*p++ != 0)
            return "Extra unparsed junk at end of rule string";
      }
   }
   return 0;
}

#ifndef QSIMPLE
void makePhases(void);
#endif

char nttable2[512];

int slowEvolveBit(int row1, int row2, int row3, int bshift) {
   return nttable[(((row2>>bshift) & 2)<<7) | (((row1>>bshift) & 2)<<6)
                | (((row1>>bshift) & 4)<<4) | (((row2>>bshift) & 4)<<3)
                | (((row3>>bshift) & 7)<<2) | (((row2>>bshift) & 1)<<1)
                |  ((row1>>bshift) & 1)<<0];
}

void fasterTable() {
   int p = 0;
   for (int row1=0; row1<8; row1++)
      for (int row2=0; row2<8; row2++)
         for (int row3=0; row3<8; row3++)
            nttable2[p++] = slowEvolveBit(row1, row2, row3, 0);
}

int evolveBitShift(int row1, int row2, int row3, int bshift) {
   return nttable2[
      (((row1 << 6) >> bshift) & 0700) +
      (((row2 << 3) >> bshift) &  070) +
      (( row3       >> bshift) &   07)];
}

int evolveBit(int row1, int row2, int row3) {
   return nttable2[
      ((row1 << 6) & 0700) +
      ((row2 << 3) &  070) +
      ( row3       &   07)];
}

int evolveRow(int row1, int row2, int row3) {
   int row4;
   int row1_s,row2_s,row3_s;
   int j,s = 0;
   int t = 0;
   int theBit;
   if (params[P_BOUNDARYSYM] == SYM_GUTTER && !gutterSkew){ /* Need to check the gutter for forbidden births */
      theBit = (row1 >> (width-1)) + ((row2 >> (width-1)) << 1) + ((row3 >> (width-1)) << 2);
      if (evolveBit(theBit, 0, theBit))
         return -1;
   }
   if (params[P_SYMMETRY] == SYM_GUTTER && !gutterSkew){    /* Need to check the gutter for forbidden births */
      theBit = (row1 & 1) + ((row2 & 1) << 1) + ((row3 & 1) << 2);
      if (evolveBit(theBit, 0, theBit))
         return -1;
   }
   if (params[P_SYMMETRY] == SYM_ODD) s = 1;
   if (params[P_BOUNDARYSYM] == SYM_UNDEF && evolveBitShift(row1, row2, row3, width - 1)) return -1;
   if (params[P_BOUNDARYSYM] == SYM_ODD) t = 1;
   if (params[P_SYMMETRY] == SYM_ASYM && evolveBit(row1 << 2, row2 << 2, row3 << 2)) return -1;
   if (params[P_SYMMETRY] == SYM_ODD || params[P_SYMMETRY] == SYM_EVEN){
      row1_s = (row1 << 1) + ((row1 >> s) & 1);
      row2_s = (row2 << 1) + ((row2 >> s) & 1);
      row3_s = (row3 << 1) + ((row3 >> s) & 1);
   }
   else {
      row1_s = (row1 << 1);
      row2_s = (row2 << 1);
      row3_s = (row3 << 1);
   }
   if (params[P_BOUNDARYSYM] == SYM_ODD || params[P_BOUNDARYSYM] == SYM_EVEN){
      row1 += ((row1 >> (width-1-t)) & 1) << (width);
      row2 += ((row2 >> (width-1-t)) & 1) << (width);
      row3 += ((row3 >> (width-1-t)) & 1) << (width);
   }
   row4 = evolveBit(row1_s, row2_s, row3_s);
   if (row4 == -1) return -1;
   for (j = 1; j < width; j++){
      theBit = evolveBitShift(row1, row2, row3, j - 1);
      if (theBit == -1) return -1;
      row4 += theBit << j;
   }
   return row4;
}

int evolveRowHigh(int row1, int row2, int row3, int bits) {
   int row4=0;
   int j,t = 0;
   int theBit;
   if (params[P_BOUNDARYSYM] == SYM_GUTTER && !gutterSkew){ /* Need to check the gutter for forbidden births */
      theBit = (row1 >> (width-1)) + ((row2 >> (width-1)) << 1) + ((row3 >> (width-1)) << 2);
      if (evolveBit(theBit, 0, theBit))
         return -1;
   }
   if (params[P_BOUNDARYSYM] == SYM_UNDEF && evolveBitShift(row1, row2, row3, width - 1)) return -1;
   if (params[P_BOUNDARYSYM] == SYM_ODD) t = 1;
   if (params[P_BOUNDARYSYM] == SYM_ODD || params[P_BOUNDARYSYM] == SYM_EVEN){
      row1 += ((row1 >> (width-1-t)) & 1) << (width);
      row2 += ((row2 >> (width-1-t)) & 1) << (width);
      row3 += ((row3 >> (width-1-t)) & 1) << (width);
   }
   for (j = width-bits; j < width; j++){
      theBit = evolveBitShift(row1, row2, row3, j - 1);
      if (theBit == -1)
         return -1;
      row4 += theBit << j;
   }
   return row4;
}

int evolveRowLow(int row1, int row2, int row3, int bits) {
   int row4;
   int row1_s,row2_s,row3_s;
   int j,s = 0;
   int theBit;
   if (params[P_SYMMETRY] == SYM_GUTTER && !gutterSkew){ /* Need to check the gutter for forbidden births */
      theBit = (row1 & 1) + ((row2 & 1) << 1) + ((row3 & 1) << 2);
      if (evolveBit(theBit, 0, theBit))
         return -1;
   }
   if (params[P_SYMMETRY] == SYM_ODD) s = 1;
   if (params[P_SYMMETRY] == SYM_ASYM && evolveBit(row1 << 2, row2 << 2, row3 << 2)) return -1;
   if (params[P_SYMMETRY] == SYM_ODD || params[P_SYMMETRY] == SYM_EVEN){
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
   if (row4 == -1) return -1;
   for (j = 1; j < bits; j++){
      theBit = evolveBitShift(row1, row2, row3, j - 1);
      if (theBit == -1) return -1;
      row4 += theBit << j;
   }
   return row4;
}

/* sortRows() is used to sort the global array valorder. 
** This array determines the order in which new rows are 
** added in the search.  The idea is that certain search
** orders give slightly quicker solutions during depthFirst()
** and LookAhead() than the naive order.
*/
void sortRows(uint16_t *theRow, uint32_t totalRows) {
   uint32_t i;
   int64_t j;
   uint16_t t;
   for (i = 1; i < totalRows; ++i){
      t = theRow[i];
      j = i - 1;
      while (j >= 0 && gcount[theRow[j]] < gcount[t]){
         theRow[j+1] = theRow[j];
         --j;
      }
      theRow[j+1] = t;
   }
}

uint16_t *makeRow(int row1, int row2);

/* getoffset() returns a pointer to a lookup table where further information */
/* is found.  It is used in getoffsetcount() and in lookAhead().             */
uint16_t *getoffset(int row12) {
   uint16_t *r = gInd3[row12];
   if (r == 0)
      r = makeRow(row12 >> width, row12 & ((1 << width) - 1));
   return r;
}
uint16_t *getoffset2(int row1, int row2) {
   return getoffset((row1 << width) + row2);
}

/* Given rows row1, row2, and row3, getoffsetcount() gives the     *\
** location (p) in the lookup table containing rows XXXX such that **
**                                                                 **
**    row1   evolves into                                          **
**    row2 ----------------> row3                                  **
**    XXXX                                                         **
**                                                                 **
\* as well as the number (n) of such rows.                         */
void getoffsetcount(int row1, int row2, int row3, uint16_t** p, int *n) {
   uint16_t *theRow = getoffset2(row1, row2);
   *p = theRow + theRow[row3];
   *n = theRow[row3+1] - theRow[row3];
}

/* Like getoffsetcount(), but only gives the number of rows.  Currently unused. */
int getcount(int row1, int row2, int row3) {
   uint16_t *theRow = getoffset2(row1, row2);
   return theRow[row3+1] - theRow[row3];
}

unsigned char *causesBirth;
row *flip;
int *gWorkConcat;       /* gWorkConcat to be parceled out between threads */
int *rowHash;
uint16_t *valorder;
void genStatCounts();

void makeFlip() {
	row theRow;
	int i;
	for (theRow = 0; theRow < (1<<width); theRow++) {
		row flippedRow = 0;
		for (i = 0; i < width; i++)
			if (theRow & (1<<i))
				flippedRow |= 1 << (width - i - 1);
		flip[theRow] = flippedRow;
	}
}

void makeTables() {
   flip = (row*)malloc(sizeof(*flip)<<width);
   makeFlip();
   causesBirth = (unsigned char*)malloc(sizeof(*causesBirth)<<width);
   gInd3 = (uint16_t **)calloc(sizeof(*gInd3),(1LL<<(width*2)));
   rowHash = (int *)calloc(sizeof(int),(2LL<<(width*2)));
   for (int i=0; i<1<<(2*width); i++)
      gInd3[i] = 0;
   for (int i=0; i<2<<(2*width); i++)
      rowHash[i] = -1;
   gcount = (uint32_t *)calloc(sizeof(*gcount), (1LL << width));
   memusage += (sizeof(*gInd3)+2*sizeof(int)) << (width*2);
   
   uint32_t i;
   for (i = 0; i < 1LLU << width; ++i) causesBirth[i] = (evolveRow(i,0,0) ? 1 : 0);
   for (i = 0; i < 1LLU << width; ++i) gcount[i] = 0;
   gWorkConcat = (int *)calloc(sizeof(int), (3LL*params[P_NUMTHREADS])<<width);
   if (params[P_REORDER] == 1)
      genStatCounts();
   if (params[P_REORDER] == 2)   /* this option currently cannot be set at runtime */
      for (int i=1; i<1<<width; i++)
         gcount[i] = 1 + gcount[i & (i - 1)];
   gcount[0] = 0xffffffff;    /* Maximum value so empty row is chosen first */
   valorder = (uint16_t *)calloc(sizeof(uint16_t), 1LL << width);
   for (int i=0; i<1<<width; i++)
      valorder[i] = (1<<width)-1-i;
   if (params[P_REORDER] != 0)
      sortRows(valorder, 1<<width);
   for (int row2=0; row2<1<<width; row2++)
      makeRow(0, row2);
}

uint16_t *bbuf;
long long int bbuf_left = 0;

/* reduce fragmentation by allocating chunks larger */
/* than needed and parceling out the small pieces.  */
uint16_t *bmalloc(int siz) {
   if (siz + (1<<width) > bbuf_left) {
      bbuf_left = (1LL << (2 * width)) + (1LL<<width);
      memusage += 2*bbuf_left;
      if (params[P_MEMLIMIT] >= 0 && memusage > memlimit) {
         printf("Aborting due to excessive memory usage\n");
         exit(1);
      }
      bbuf = (uint16_t *)calloc(sizeof(uint16_t), bbuf_left);
   }
   uint16_t *r = bbuf;
   bbuf += siz;
   bbuf_left -= siz;
   return r;
}

void unbmalloc(int siz) {
   bbuf -= siz;
   bbuf_left += siz;
}

/* hashRow() is used to identify if we've already */
/* built an identical part of the lookup table.   */
unsigned int hashRow(uint16_t *theRow, int siz) {
   unsigned int h = 0;
   for (int i=0; i<siz; i++)
      h = h * 3 + theRow[i];
   return h;
}

uint16_t *makeRow(int row1, int row2) {
   int good = 0;
   /* Set up gWork for this particular thread */
   int *gWork = gWorkConcat + ((3LL * omp_get_thread_num()) << width);
   int *gWork2 = gWork + (1 << width);
   int *gWork3 = gWork2 + (1 << width);
   
   /* For each row3, find all row4 such that   *\
   **                                          **
   **      row1   evolves into                 **
   **      row2 ----------------> row4         **
   **      row3                                **
   **                                          **
   ** list of row4s is stored in gwork and     **
   \* corresponding row3s are stored in gwork2 */
   if (width < 4) {
      for (int row3=0; row3<1<<width; row3++)
         gWork3[row3] = evolveRow(row1, row2, row3);
   } else {
      int lowbitcount = (width >> 1) + 1;
      int hibitcount = ((width + 1) >> 1) + 1;
      int hishift = lowbitcount - 2;
      int lowcount = 1 << lowbitcount;
      for (int row3=0; row3<1<<lowbitcount; row3++)
         gWork2[row3] = evolveRowLow(row1, row2, row3, lowbitcount-1);
      for (int row3=0; row3<1<<width; row3 += 1<<hishift)
         gWork2[lowcount+(row3>>hishift)] =
                        evolveRowHigh(row1, row2, row3, hibitcount-1);
      for (int row3=0; row3<1<<width; row3++)
         gWork3[row3] = gWork2[row3 & ((1<<lowbitcount) - 1)] |
                        gWork2[lowcount+(row3 >> hishift)];
   }
   for (int row3i = 0; row3i < 1<<width; row3i++) {
      int row3 = valorder[row3i];
      int row4 = gWork3[row3];
      if (row4 < 0)
         continue;
      gWork2[good] = row3;
      gWork[good++] = row4;
   }
   
   /* bmalloc, unbmalloc, and all operations that read from or write to */
   /* theRow, rowHash, and gInd3 must be included in a critical region. */
   uint16_t *theRow;
   #pragma omp critical(updateTable)
   {
      theRow = bmalloc((1+(1<<width)+good));
      for (int row3=0; row3 < 1<<width; row3++)
         theRow[row3] = 0;
      theRow[0] = 1 + (1 << width);
      for (int row3=0; row3 < good; row3++)
         theRow[gWork[row3]]++;
      theRow[1<<width] = 0;
      for (int row3=0; row3 < (1<<width); row3++)
         theRow[row3+1] += theRow[row3];
      for (int row3=good-1; row3>=0; row3--) {
         int row4 = gWork[row3];
         theRow[--theRow[row4]] = gWork2[row3];
      }
      unsigned int h = hashRow(theRow, 1+(1<<width)+good);
      h &= (2 << (2 * width)) - 1;
      while (1) {
         if (rowHash[h] == -1) {
            rowHash[h] = (row1 << width) + row2;
            break;
         }
         /* Maybe two different row12s result in the exact same rows for the */
         /* lookup table. This prevents two different threads from trying to */
         /* build the same part of the lookup table.                         */
         if (memcmp(theRow, gInd3[rowHash[h]], 2*(1+(1<<width)+good)) == 0) {
            theRow = gInd3[rowHash[h]];
            unbmalloc(1+(1<<width)+good);
            break;
         }
         h = (h + 1) & ((2 << (2 * width)) - 1);
      }
      
      gInd3[(row1<<width)+row2] = theRow;
   }
   
   return theRow;
}

/*   We calculate the stats using a 2 * 64 << width array.  We use a
**   leading 1 to separate them.  Index 1 aaa bb cc dd represents
**   the count for a result of aaa when the last two bits of row1, row2,
**   and row3 were bb, cc, and dd, respectively.  We have to manage
**   the edge conditions appropriately.
*/
void genStatCounts() {
   int *cnt = (int*)calloc((128 * sizeof(int)), 1LL << width);
   for (int i=0; i<128<<width; i++)
      cnt[i] = 0;
   int s = 0;
   if (params[P_SYMMETRY] == SYM_ODD)
      s = 2;
   else if (params[P_SYMMETRY] == SYM_EVEN)
      s = 1;
   else
      s = width + 2;
   /* left side: never permit generation left of row4 */
   for (int row1=0; row1<2; row1++)
      for (int row2=0; row2<2; row2++)
         for (int row3=0; row3<2; row3++)
            if (evolveBit(row1, row2, row3) == 0)
               cnt[(1<<6) + (row1 << 4) + (row2 << 2) + row3]++;
   for (int nb=0; nb<width; nb++) {
      for (int row1=0; row1<8; row1++)
         for (int row2=0; row2<8; row2++)
            for (int row3=0; row3<8; row3++) {
               if (nb == width-1)
                  if ((((row1 >> s) ^ row1) & 1) ||
                      (((row2 >> s) ^ row2) & 1) ||
                      (((row3 >> s) ^ row3) & 1))
                     continue;
               int row4b = evolveBit(row1, row2, row3);
               for (int row4=0; row4<1<<nb; row4++)
                  cnt[(((((1<<nb) + row4) << 1) + row4b) << 6) +
                    ((row1 & 3) << 4) + ((row2 & 3) << 2) + (row3 & 3)] +=
                     cnt[(((1<<nb) + row4) << 6) +
                       ((row1 >> 1) << 4) + ((row2 >> 1) << 2) + (row3 >> 1)];
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
                       (row1 << 4) + (row2 << 2) + row3];
   free(cnt);
}

/* ====================================================== */
/*  Hash table for detecting equivalent partial patterns  */
/* ====================================================== */

void resetHash() {
   if (hash != 0) memset(hash,0,4*HASHSIZE);
}

int hashPhase = 0;

static inline long hashFunction(node b, row r) {
   long h = r;
   if (params[P_SYMMETRY] == SYM_ASYM) h += flip[r];
   int i;
   for (i = 0; i < nRowsInState; i++) {
      h = (h * 269) + ROW(b);
      if (params[P_SYMMETRY] == SYM_ASYM) h += flip[ROW(b)];
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

static inline int sameFlipped(node p, node q, row r) {
   int i;
   for (i = 0; i < nRowsInState; i++) {
      if (p >= QSIZE || q >= QSIZE || EMPTY(p) || EMPTY(q)) return 0;   /* sanity check */
      if (flip[ROW(p)] != r) return 0;
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
      if (hashNode == 0)
         return 0;
      else if (same(hashNode,b,r))
         return 1;
      else if (params[P_SYMMETRY] == SYM_ASYM && sameFlipped(hashNode,b,r))
         return 1;
   }
   return 0;
}

/* Add node (NOT child) to hash table */
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

/* Avoid Intel shift bug */
static inline unsigned long
safeShift(unsigned long r, int i) {
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
int bufferPattern(node b, row *pRows, int nodeRow, uint32_t lastRow, int printExpected) {
   node c;
   int nrows = 0;
   int swidth;
   int sxsNeeded;
   int p, i, j, margin;
   unsigned long *srows, *ssrows;

   uint32_t currRow = lastRow;
   int nDeepRows = 0;
   int nodeDiff;

   if (pRows != NULL){
      while (pRows[currRow] == 0){
         if (currRow == 0){
            if (!printExpected) return 0;
            printf("Success called on search root!\n");
            aborting = 1;
            return 0;
         }
         currRow--;
      }
      nDeepRows = (currRow / period) - 1;
      nodeDiff = nodeRow - period - (currRow%period);
      nodeRow -= nodeDiff;

      for (j = 0; j < nodeDiff; j++){
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
            if (!printExpected) return 0;
            printf("Success called on search root!\n");
            aborting = 1;
            return 0;
         }
      }
   }
   if (nrows < 0) nrows = 0;
   
   for (p = 0; p < period-1; p++) b = PARENT(b);
   if (b == 0) {
      if (!printExpected) return 0;
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
      if (nDeepRows > 0){
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
      switch(params[P_SYMMETRY]) {
         case SYM_ASYM:
            srows[i] = r;
            break;

         case SYM_ODD:
            srows[i] = r << (MAXWIDTH - 1);
            ssrows[i] = r >> (32 - (MAXWIDTH - 1));
            for (j = 1; j < MAXWIDTH; j++)
               if (r & (1<<j))
                  srows[i] |= 1 << (MAXWIDTH - 1 - j);
            break;

         case SYM_EVEN:
            srows[i] = r << MAXWIDTH;
            ssrows[i] = r >> (32 - MAXWIDTH);
            for (j = 0; j < MAXWIDTH; j++)
               if (r & (1<<j))
                  srows[i] |= 1 << (MAXWIDTH - 1 - j);
            break;

         case SYM_GUTTER:
            srows[i] = r << (MAXWIDTH + 1);
            ssrows[i] = r >> (32 - (MAXWIDTH + 1));
            for (j = 0; j < MAXWIDTH; j++)
               if (r & (1<<j))
                  srows[i+gutterSkew] |= 1 << (MAXWIDTH - 1 - j);
            break;

         default:
            printError("unexpected symmetry type in success()");
            return 0;
      }
   }
   
   /* normalize nrows to only include blank rows */
   nrows += MAXWIDTH;
   while (nrows>0 && srows[nrows-1] == 0 && ssrows[nrows-1] == 0) nrows--;
   while (srows[0] == 0 && ssrows[0] == 0 && nrows>0) {
      srows++;
      ssrows++;
      nrows--;
   }
   
   /* sanity check: is the pattern nonempty? */
   int allEmpty = 1;
   for (i = 0; i < nrows; i++){
      if (srows[i]){
         allEmpty = 0;
         break;
      }
   }
   if (allEmpty) return 0;
   
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

   /* make sure we didn't just output the exact same pattern */
   if (printExpected){
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
      memcpy(oldsrows, srows, nrows * sizeof(unsigned long));
      memcpy(oldssrows, ssrows, nrows * sizeof(unsigned long));
   }
   
   /* Buffer output */
   patternBuf = (char*)realloc(patternBuf, ((2 * MAXWIDTH + 4) * sxsAllocRows + 300) * sizeof(char));
   
   sprintf(patternBuf,"x = %d, y = %d, rule = %s\n", swidth - margin, nrows, baseRule);
   
   int theBufRow = -1;
   while (theBufRow++ < nrows){
      bufRow(ssrows[theBufRow], srows[theBufRow], 0);
   }
   RLEcount = 1;     /* prevents erroneous printing of '2' at end of RLE */
   RLEchar = '!';
   bufRLE('\0');
   sprintf(patternBuf+strlen(patternBuf),"\n");
   
   if (printExpected){
      numFound++;
      if (params[P_NUMSHIPS] > 0){
         if (--params[P_NUMSHIPS] == 0) aborting = 3;  /* use 3 to flag that we reached ship limit */
      }
   }
   
   return 1;
}

#ifndef QSIMPLE
void makeSubperiodTables(void);
int subperiodic(node x, row *pRows, int nodeRow, uint32_t lastRow);
#endif

void success(node b, row *pRows, int nodeRow, uint32_t lastRow) {
#ifndef QSIMPLE
   if (subperiodic(b, pRows, nodeRow, lastRow)) return;
#endif
   if (bufferPattern(b, pRows, nodeRow, lastRow, 1))
      printf("\n%s\n",patternBuf);
   fflush(stdout);
}

/* Test if this is a node at which we can stop */
int terminal(node n) {
   int p;

   for (p = 0; p < period; p++) {   /* last row in each phase must be zero */
      if (ROW(n) != 0) return 0;
      n = PARENT(n);
   }

   for (p = 0; p < period; p++) {
      if (causesBirth[ROW(n)]) return 0;
      n = PARENT(n);
   }
   return 1;
}

/* ================================================ */
/*  Queue of partial patterns still to be examined  */
/* ================================================ */

node qHead,qTail;

/* Queue dimensions required during save/restore */
node qStart; /* index of first node in queue */
node qEnd;   /* index of first unused node after end of queue */

/* Maintain phase of queue nodes.  After dequeue(), the global variable phase
** gives the phase of the dequeued item.  If the queue is compacted, this information
** needs to be reinitialized by a call to rephase(), after which phase will not be
** valid until the next call to dequeue().  Variable nextRephase points to the next
** node for which dequeue will need to increment the phase. Phase is not maintained
** when treating queue as a stack (using pop()) -- caller must do it in that case.
** It's ok to change phase since we maintain a separate copy in queuePhase.
*/
int queuePhase = 0;
node nextRephase = 0;

void rephase() {
   node x, y;
   while (qHead < qTail && EMPTY(qHead)) qHead++;   /* skip empty queue cells */
   x = qHead;   /* find next item in queue */
   queuePhase = period - 1;
   while (x != 0) {
      x = PARENT(x);
      queuePhase++;
   }
   queuePhase %= period;

   /* Now walk forward through queue finding breakpoints between each generation. */
   /* invariants: y is always the first in its generation                         */
   x = 0; y = 0;
   while (y <= qHead) {
      ++x;
      if (x >= qTail || (!EMPTY(x) && PARENT(x) >= y)) y = x;
   }
   nextRephase = y;
}

/* peekPhase() returns the phase of an element in the queue.  This only */
/* works for queue elements and should NOT be used on other nodes.      */
int peekPhase(node i) {
   return (i < nextRephase? queuePhase : (queuePhase+1)%period);
}

/* Test queue status */
static inline int qIsEmpty() {
   while (qHead < qTail && EMPTY(qHead)){
      ++qHead;
      ++deepQHead;
   }
   return (qTail == qHead);
}

void qFull() {
    if (aborting != 2) {
      printf("Exceeded %lu node limit, search aborted\n", QSIZE);
      fflush(stdout);
      aborting = 2;
   }
}

static inline void enqueue(node b, row r) {
   node tempQTail = qTail;
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
   
   /* update tail of parallel queue, but don't set value */
   deepQTail += qTail - tempQTail;
   deepRowIndices[deepQTail] = 0;
}

static inline node dequeue() {
   oldDeepQHead = deepQHead;  /* Save old parallel queue head for use in process() */
   while (qHead < qTail && EMPTY(qHead)){
      ++qHead;
      ++deepQHead;
   }
   if (qHead >= nextRephase) {
      queuePhase = (queuePhase+1)%period;
      nextRephase = qTail;
   }
   phase = queuePhase;
   ++deepQHead;
   return qHead++;
}

/* Not used, but could be useful in debugging */
static inline void pop() {
   qTail--;
   while (qTail > qHead && EMPTY(qTail-1)) qTail--;
}

void resetQ() {
   qHead = qTail = 0; deepQHead = deepQTail = 0;
}

//static inline int qTop() { return qTail - 1; }

/* =================== */
/*  Dump search state  */
/* =================== */

int dumpNum = 1;
char dumpFile[256];
const char *dumpRoot = "dump-@time-";
char loadDumpRoot[251];    /* used for loading dump root from file */
char trueDumpRoot[251];
time_t lastDumpTime;

int dumpFlag = 0;    /* Dump status flags, possible values follow */
#define DUMPRESET   (0)
#define DUMPPENDING (1)
#define DUMPFAILURE (2)
#define DUMPSUCCESS (3)

int dumpMode;  /* separate from params[P_DUMPMODE] for splitting */
#define D_DISABLED   (0)
#define D_OVERWRITE  (1)
#define D_SEQUENTIAL (2)

void parseDumpRoot() {
   char tempStr[MAXDUMPROOT + 10]  = {'\0'};
   char tempRule[151] = {'\0'};
   char *r, *t;
   const char *s;
   memset(trueDumpRoot, '\0', 251);
   
   /* replace first "@time" with hex timestamp */
   if ((s = strstr(dumpRoot, "@time"))){
      memcpy(tempStr, dumpRoot, s - dumpRoot);
      sprintf(trueDumpRoot, "%s%06lx%s", tempStr, time(NULL) & 0xffffff, s+5);
   }
   else
      memcpy(trueDumpRoot, dumpRoot, MAXDUMPROOT + 9);
   
   /* replace first "@rule" with rule string ('/' replaced with '_') */
   memcpy(tempStr, trueDumpRoot, MAXDUMPROOT + 9);
   if ((t = strstr(tempStr, "@rule"))){
      memcpy(tempRule, rule, 150);
      r = strchr(tempRule, '/');
      *r = '_';
      *t = '\0';
      sprintf(trueDumpRoot, "%s%s%s", tempStr, tempRule, t+5);
   }
   
   /* replace additional occurrences of '@' with '_' */
   r = trueDumpRoot - 1;
   while (*(++r) != '\0')
      if (*r == '@')
         *r = '_';
   
   dumpRoot = trueDumpRoot;
}

FILE * openDumpFile() {
   FILE * fp;
   
   if (dumpMode == D_OVERWRITE) {
      sprintf(dumpFile, "%s%s", dumpRoot, (++dumpNum)%2 ? "gold" : "blue");
      return fopen(dumpFile, "w");
   }
   else if (dumpMode == D_SEQUENTIAL){
      while (dumpNum < DUMPLIMIT) {
         sprintf(dumpFile, "%s%05d", dumpRoot, dumpNum++);
         if ((fp = fopen(dumpFile, "r")))
            fclose(fp);
         else
            return fopen(dumpFile, "w");
      }
      if (dumpNum == DUMPLIMIT){
         dumpMode = D_OVERWRITE;
         return openDumpFile();
      }
   }
   return (FILE *) 0;
}

void dumpState() {
   FILE * fp;
   unsigned long long i,j;
   dumpFlag = DUMPFAILURE;
   if (!(fp = openDumpFile())) return;
   fprintf(fp,"%lu\n",FILEVERSION);
   fprintf(fp,"%s\n",rule);
   fprintf(fp,"%s\n",dumpRoot);
   for (j = 0; j < NUM_PARAMS; ++j)
      fprintf(fp,"%d\n",params[j]);
   fprintf(fp,"%d\n",width);
   fprintf(fp,"%d\n",period);
   fprintf(fp,"%d\n",offset);
   fprintf(fp,"%d\n",lastDeep);
   if (params[P_DUMPMODE] == D_SEQUENTIAL)
      fprintf(fp,"1\n");
   else
      fprintf(fp,"%d\n",dumpNum%2);
   fprintf(fp,"%"PRIu32"\n",qHead-qStart);
   fprintf(fp,"%"PRIu32"\n",qEnd-qStart);
   for (i = qStart; i < qEnd; ++i)
      fprintf(fp,"%"PRIu16"\n",rows[i]);
   for (i = 0; i < QSIZE; ++i){
      if (deepRowIndices[i]){
         if (deepRowIndices[i] > 1){
            for (j = 0; j < deepRows[deepRowIndices[i]][0] + 1LU + 2LU; ++j){
               fprintf(fp,"%"PRIu16"\n",deepRows[deepRowIndices[i]][j]);
            }
         }
         else {
            fprintf(fp,"0\n");
            j = 0;
            while (deepRowIndices[i] <= 1 && i < QSIZE){
               if (deepRowIndices[i] == 1) ++j;
               ++i;
            }
            fprintf(fp,"%llu\n",j);
            if (i == QSIZE) break;
            --i;
         }
      }
   }
   fclose(fp);
   dumpFlag = DUMPSUCCESS;
}

/* ================================= */
/*  Compaction of nearly full queue  */
/* ================================= */

void putnum(long unsigned n) {
   char suffix;
   if (n >= 1000000) {
      n /= 100000;
      suffix = 'M';
   } else if (n >= 1000) {
      n /= 100;
      suffix = 'k';
   } else {
      printf("%lu", n);
      return;
   }

   if (n >= 100) printf("%lu", n/10);
   else printf("%lu.%lu", n/10, n%10);
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

/* doCompact() has two parts.  The first part compresses the
** queue.  The second part consists of the last loop which
** converts parent bits back to parent pointers.  The search
** state may be saved in between.
*/
void doCompactPart1() {
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
   
   /* make a pass forwards converting parent pointers to offset from prev parent ptr. */
   /* note that after unused nodes are eliminated, all these offsets are zero or one. */
   y = 0;
   for (x = 0; x < qTail; x++) if (!EMPTY(x)) {
      if (PARENT(x) == y) rows[x] = ROW(x);
      else {
         y = PARENT(x);
         rows[x] = (1<<width) + ROW(x);
      }
   }
   
   /* Make a pass backwards compacting gaps.
   **
   ** For most times we run this, it could be combined with the next phase, but
   ** every once in a while the process of repacking the remaining items causes them
   ** to use *MORE* space than they did before they were repacked (because of the need
   ** to leave empty space when ROFFSET gets too big) and without this phase the repacked
   ** stuff overlaps the not-yet-repacked stuff causing major badness.
   ** 
   ** For this phase, y points to the current item to be repacked, and x points
   ** to the next free place to pack an item.
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

void doCompactPart2() {
   node x,y;
   uint32_t i, j;
   int k;
   
   /* Make a pass forwards converting parent bits back to parent pointers.
   ** 
   ** For this phase, x points to the current item to be repacked, and y points
   ** to the parent of the previously repacked item.
   ** After the previous pass, x is initialized to first nonempty item,
   ** and all items after x are nonempty. 
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
   
   /* Repack nonzero depth-first extension indices to the end of the queue */
   j = QSIZE - 1;
   for (i = QSIZE; i > 0; --i){
      if (deepRowIndices[i - 1]){
         deepRowIndices[j] = deepRowIndices[i - 1];
         deepRowIndices[i - 1] = 0;
         --j;
      }
   }
   
   /* Sanity check: extension queue should not take up all available space */
   if (deepRowIndices[0]){
      fprintf(stderr,"Error: extension queue has too many elements.\n");
      exit(1);
   }
   
   /* Respace depth-first extension queue to match node queue */
   i = 0;
   j = 0;
   while (!deepRowIndices[j] && j < QSIZE) ++j;
   for (x = qHead; x < qTail && j < QSIZE; ++x){
      if (EMPTY(x)){
         ++i;
         continue;
      }
      
      deepRowIndices[i] = deepRowIndices[j];
      
      /* Sanity check: do the extension rows match the node rows? */
      if (deepRowIndices[j] > 1){
         y = x;
         for (k = 0; k < 2*period; ++k){
            uint16_t startRow = deepRows[deepRowIndices[j]][1] + 1;
            if (deepRows[deepRowIndices[j]][startRow - k] != ROW(y)){
               fprintf(stderr, "Warning: non-matching rows detected at node %u in doCompactPart2()\n",x);
               free(deepRows[deepRowIndices[j]]);
               deepRows[deepRowIndices[j]] = 0;
               deepRowIndices[i] = 0;
               break;
            }
            y = PARENT(y);
         }
      }
      if (j > i) deepRowIndices[j] = 0;
      ++i;
      ++j;
   }
   for (j = qTail - qHead; j < QSIZE; ++j){
      deepRowIndices[j] = 0;
   }
   deepQHead = 0;
   deepQTail = qTail - qHead;
}

void doCompact() {
   /* make sure we still have something left in the queue */
   if (qIsEmpty()) {
      qTail = qHead = 0;   /* nothing left, make an extremely compact queue */
      return;
   }
   /* First loop of part 1 requires qTail-1 to be non-empty.  Make it so. */
   while (EMPTY(qTail-1))
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
#ifndef QSIMPLE
   if (params[P_CACHEMEM] == 0) return 0;
#endif
   unsigned long long h = (unsigned long long)p1 +
      17 * (unsigned long long)p2 + 257 * (unsigned long long)p3 +
      513 * abn;
   h = h + (h >> 15);
   h &= (cachesize-1);
   cacheentry *ce = &(cache[omp_get_thread_num()][h]);
   if (ce->p1 == p1 && ce->p2 == p2 && ce->p3 == p3 && ce->abn == abn)
      return -2 + ce->r;
   ce->p1 = p1;
   ce->p2 = p2;
   ce->p3 = p3;
   ce->abn = abn;
   return h;
}

void setkey(int h, int v) {
#ifndef QSIMPLE
   if (params[P_CACHEMEM])
#endif
      cache[omp_get_thread_num()][h].r = v;
}
#endif

/* ========================== */
/*  Primary search functions  */
/* ========================== */

void process(node theNode);
int depthFirst(node theNode, uint16_t howDeep, uint16_t **pInd, int *pRemain, row *pRows, _Atomic int *remainingItems, _Atomic int *forceExit);

static void deepen() {
   /* compute amount to deepen, apply reduction if too deep */
   int deepeningAmount;
   int i = currentDepth();
   
   if (i >= lastDeep) deepeningAmount = MINDEEP;
   else deepeningAmount = lastDeep + MINDEEP - i;   /* go at least MINDEEP deeper */
   
   if (params[P_FIRSTDEEP]){
      deepeningAmount = params[P_FIRSTDEEP];
      params[P_FIRSTDEEP] = 0;
   }
   
   lastDeep = i + deepeningAmount;
   
   /* report what's happening */
   printf("%d, deepening %d, ", i, deepeningAmount);
   putnum(qTail - qHead);
   printf("/");
   putnum(qTail);
   fflush(stdout);
   
   _Atomic int remainingItems = 0;
   _Atomic int forceExit = 0;
   atomic_store_explicit(&remainingItems, qTail - qHead, memory_order_seq_cst);
   atomic_store_explicit(&forceExit, 0, memory_order_seq_cst);

   /* go through queue, deepening each one */
   #pragma omp parallel
   {
      uint16_t **pInd;
      int *pRemain;
      row *pRows;
      
      pInd = (uint16_t**)calloc((deepeningAmount + 4 * params[P_PERIOD]), (long long)sizeof(*pInd));
      pRemain = (int*)calloc((deepeningAmount + 4 * params[P_PERIOD]), (long long)sizeof(*pRemain));
      pRows = (row*)calloc((deepeningAmount + 4 * params[P_PERIOD]), (long long)sizeof(*pRows));
      
      long long j;
      #pragma omp for schedule(dynamic, CHUNK_SIZE)
      for (j = qHead; j < qTail; j++) {
         if (!EMPTY(j) && !depthFirst((node)j, (uint16_t)deepeningAmount, pInd, pRemain, pRows, &remainingItems, &forceExit))
            MAKEEMPTY(j);
         atomic_fetch_sub_explicit(&remainingItems, 1, memory_order_relaxed);
      }
      free(pInd);
      free(pRemain);
      free(pRows);
   }
   
   /* before reporting new queue size, shrink tree back down */
   printf(" -> ");
   fflush(stdout);
   
   /* signal time for dump */
   if (params[P_DUMPMODE] != D_DISABLED && time(NULL) - lastDumpTime > params[P_DUMPINTERVAL]){
      dumpFlag = DUMPPENDING;
      time(&lastDumpTime);
   }
   
   doCompact();
   
   /* now finish report */
   putnum(qTail - qHead);
   printf("/");
   putnum(qTail);
   printf("\n");
   
   /* Report successful/unsuccessful dump */
   if (dumpFlag == DUMPSUCCESS) {
      timeStamp();
      printf("State dumped to %s\n",dumpFile);
      if (dumpNum == DUMPLIMIT){
         timeStamp();
         printf("Sequential dump limit reached.  Changing to overwrite mode.\n");
      }
   }
   else if (dumpFlag == DUMPFAILURE) {
      timeStamp();
      printf("State dump unsuccessful\n");
   }
   dumpFlag = DUMPRESET;
   
   fflush(stdout);
}

static void breadthFirst() {
   while (!aborting && !qIsEmpty()){
      if (qTail - qHead >= (1LLU<<params[P_DEPTHLIMIT]) || qTail >= QSIZE - QSIZE/16){
         timeStamp();
         printf("Queue full, depth ");
         deepen();
      }
      else if (params[P_EVERYDEPTH] && qHead == nextRephase){
         timeStamp();
         printf("Depth ");
         deepen();
      }
      else
         process(dequeue());
   }
}

void saveDepthFirst(node theNode, uint16_t startRow, uint16_t howDeep, row *pRows) {
   uint32_t theDeepIndex;
   #pragma omp critical(findDeepIndex)
   {
      theDeepIndex = 2;
      while (theDeepIndex < 1LLU << (params[P_DEPTHLIMIT] + 1) && deepRows[theDeepIndex]) ++theDeepIndex;
      if (theDeepIndex == 1LLU << (params[P_DEPTHLIMIT] + 1)){
         fprintf(stderr,"Error: no available extension indices.\n");
         aborting = 1;
      }
      if (!aborting){
         deepRows[theDeepIndex] = (row*)calloc( startRow + howDeep + 1 + 2,
                                                sizeof(**deepRows) );
      }
   }
   if (aborting) return;
   
   memcpy( deepRows[theDeepIndex] + 2,
           pRows,
           (startRow + howDeep + 1) * (long long)sizeof(**deepRows) );
   
   deepRows[theDeepIndex][0] = startRow + howDeep;
   deepRows[theDeepIndex][1] = startRow;
   
   deepRowIndices[deepQHead + theNode - qHead] = theDeepIndex;
}

/* ========================== */
/*  Print usage instructions  */
/* ========================== */

void printHelp() {
   printf("Usage:    ./qfind "
#ifndef QSIMPLE
                             "-v <velocity> "
#endif
                                        "-w <width> -s <symmetry> [options...]\n"
          "       or\n"
          "          ./qfind -l <file> [options...]\n");
   printf("\n");
   printf("qfind is a program that searches for orthogonal spaceships and waves in Life\n"
          "and related cellular automata.  Options are read left to right, with subsequent\n"
          "occurrences of the same option overwriting the previous value.\n");
   printf("\n");
#ifdef QSIMPLE
   printf("When using qfind-s, the period and offset must be set within the code before it \n"
          "is compiled.  You have compiled with\n");
   printf("\n");
   printf("Period: %d\n",PERIOD);
   printf("Offset: %d\n",OFFSET);
   printf("\n");
#endif
   printf("Required (except when loading from a saved state):\n");
#ifndef QSIMPLE
   printf("  -v, --velocity <velocity>     written in the form <translation>c/<period>\n");
#endif
   printf("  -w, --width <number>          logical width (full width depends on symmetry)\n");
   printf("  -s, --symmetry <(asymmetric|odd|even|gutter)>  spaceship symmetry type\n");
   printf("\n");
   printf("Search options:\n");
   printf("  -r, --rule <rule>             cellular automaton rule written in Hensel\n"
          "                                notation (Default: B3/S23)\n"
          "                                '~' is used to specify a list of forbidden\n"
          "                                conditions.  For example, -r B3~6c7/S23~8\n"
          "                                searches in B3/S23 for ships that never contain\n"
          "                                the B6c, B7, or S8 neighborhoods.\n");
#ifdef _OPENMP
   printf("  -t, --threads <number>        number of threads during deepening (default: 1)\n");
#endif
   printf("  -f, --found <number>          maximum number of spaceships to output\n");
   printf("  -i, --increment <number>      minimum deepening increment (default: 3)\n");
   printf("  -g, --min-extension <number>  minimum length of saved extensions\n");
   printf("  -n, --first-depth <number>    depth of first deepening step\n");
   printf("      --fixed-depth <number>    deepen at every new depth by the given amount\n");
   printf("  -e, --extend <filename>       file containing the initial rows for a search.\n"
          "                                Use the Golly script get-rows.lua to easily\n"
          "                                generate the initial rows file.\n");
#ifdef _OPENMP
   printf("  (--enable-early-exit|--disable-early-exit)\n"
          "                                enable/disable early exit during deepening step\n"
          "                                when threads become idle (default: enabled)\n");
#endif
   printf("\n");
   printf("Memory options:\n");
   printf("  -c, --cache-mem <number>      allocate N megabytes per thread for lookahead\n"
          "                                cache (default: %d if speed is greater than c/5\n"
          "                                and disabled otherwise)\n"
          "                                Use -c 0 to disable lookahead caching.\n",DEFAULT_CACHEMEM);
   printf("  -m, --mem-limit <number>      limits lookup table memory to N megabytes\n");
   printf("  -q, --queue-bits <number>     set BFS queue size to 2^N nodes (default: %d)\n", QBITS);
   printf("  -h, --hash-bits <number>      set hash table size to 2^N nodes (default: %d)\n"
          "                                Use -h 0 to disable duplicate elimination.\n", HASHBITS);
   printf("  -b, --base-bits <number>      groups 2^N queue entries to an index node\n"
          "                                (default: 4)\n");
   printf("\n");
   printf("Save/load options:\n");
   printf("  -d, --dump-root <string>      dump filename prefix\n");
   printf("  -a, --dump-interval <number>  wait at least N seconds between dumps\n");
   printf("      --dump-mode <(overwrite|sequential|disabled)>\n"
          "                                set dump mode\n");
   printf("  -l, --load <filename>         load search state from the given dump file\n");
   printf("  -j, --split <number>          split loaded search state into at most N files\n");
   printf("  -p, --preview                 preview partial results from the loaded state\n");
   printf("\n");
   printf("Output options (enabled by default):\n");
#ifndef QSIMPLE
   printf("  (--enable-subperiod|--disable-subperiod)    enable/disable printing of\n"
          "                                              subperiodic results\n");
#endif
   printf("  (--enable-deep-print|--disable-deep-print)  enable/disable printing ships\n"
          "                                              during deepening step\n");
   printf("  (--enable-longest|--disable-longest)        enable/disable printing longest\n"
          "                                              partial result at end of search\n");
   printf("\n");
   printf("Wave options:\n");
   printf("  -o, --boundary-sym <(disabled|odd|even|gutter)>  boundary symmetry type for\n"
          "                                                   wave searches\n");
   printf("\n");
   printf("Documentation options:\n");
   printf("  --help                        print usage instructions and exit\n");
#ifndef QSIMPLE
   printf("\n");
   printf("Example search:\n"
          "    ./qfind -v c/5 -w 9 -s even -r B3/S23 -t 2\n\n"
          "  Searches Life (rule B3/S23) for c/5 orthogonal spaceships with even\n"
          "  bilateral symmetry and logical width 9 (full width 18) using two threads.\n\n");
#endif
   exit(0);
}

/* ======================== */
/*  Echo loaded parameters  */
/* ======================== */

void echoParams() {
   printf("\n");
   printf("Rule: %s\n",rule);
   printf("speed: ");
   if (params[P_OFFSET] != 1) printf("%d",params[P_OFFSET]);
   printf("c/%d\n", params[P_PERIOD]);
   printf("Width: %d\n", params[P_WIDTH]);
   if (params[P_SYMMETRY] == SYM_ASYM) printf("Symmetry: asymmetric\n");
   else if (params[P_SYMMETRY] == SYM_ODD) printf("Symmetry: odd\n");
   else if (params[P_SYMMETRY] == SYM_EVEN) printf("Symmetry: even\n");
   else if (params[P_SYMMETRY] == SYM_GUTTER) printf("Symmetry: gutter\n");
   if (params[P_BOUNDARYSYM] != SYM_UNDEF){
      printf("Wave search enabled\nBoundary symmetry: ");
      if (params[P_BOUNDARYSYM] == SYM_ODD) printf("odd\n");
      else if (params[P_BOUNDARYSYM] == SYM_EVEN) printf("even\n");
      else if (params[P_BOUNDARYSYM] == SYM_GUTTER) printf("gutter\n");
   }
#ifndef QSIMPLE
   if (params[P_FULLPERIOD] && gcd(period,offset)>1) printf("Suppress subperiodic results\n");
#endif
   if (params[P_DUMPMODE] != D_DISABLED){
      printf("Dump interval: %d second%s\n", params[P_DUMPINTERVAL], params[P_DUMPINTERVAL] == 1 ? "" : "s");
      printf("Dump mode: %s\n", params[P_DUMPMODE] == D_OVERWRITE ? "overwrite" : "sequential");
   }
   else
      printf("Dump disabled\n");
   printf("Queue size: 2^%d\n",params[P_QBITS]);
   printf("Hash table size: 2^%d\n",params[P_HASHBITS]);
   if (params[P_EVERYDEPTH])
      printf("Fixed deepening amount: %ld\n",
               params[P_FIRSTDEEP] ? (long)params[P_FIRSTDEEP] : lastDeep - currentDepth());
   else
      printf("Minimum deepening increment: %d\n",MINDEEP);
   if (params[P_PRINTDEEP] == 0) printf("Output disabled while deepening\n");
#ifndef NOCACHE
   if (params[P_CACHEMEM])
      printf("Cache memory per thread: %d megabytes\n", params[P_CACHEMEM]);
   else
      printf("Lookahead caching disabled\n");
#endif
   if (params[P_MEMLIMIT] >= 0) printf("Memory limit: %d megabytes\n",params[P_MEMLIMIT]);
#ifdef _OPENMP
   printf("Number of threads: %d\n",params[P_NUMTHREADS]);
#endif
   if (params[P_MINEXTENSION]) printf("Save depth-first extensions of length at least %d\n",params[P_MINEXTENSION]);
   if (params[P_LONGEST] == 0) printf("Printing of longest partial result disabled\n");
   printf("\n");
}

/* ========================= */
/*  Preview partial results  */
/* ========================= */

static void preview(/*int allPhases*/) {
   node i,j;
   row *pRows;
   // int ph;  /* used for allPhases option */
   // node k;  /* used for allPhases option */

   for (i = qHead; (i<qTail) && EMPTY(i); i++);
   for (j = qTail-1; (j>i) && EMPTY(j); j--);
   if (j<i) return;
   
   while (j>=i && !aborting) {
      if (!EMPTY(j)) {
         uint32_t theDeepIndex = deepRowIndices[deepQHead + j - qHead];
         
         if (theDeepIndex > 1){
            pRows = (row*) malloc((2*period + 1 + deepRows[theDeepIndex][0]
                                   - deepRows[theDeepIndex][1] + 1) * (long long)sizeof(*pRows));
            int m;
            node x = j;
            for (m = 2*period; m >= 0; --m){
               pRows[m] = ROW(x);
               x = PARENT(x);
            }
            memcpy(pRows + 2*period+1,
                   deepRows[theDeepIndex] + 2 + deepRows[theDeepIndex][1],
                   (deepRows[theDeepIndex][0] - deepRows[theDeepIndex][1] + 1) * (long long)sizeof(*pRows));
            uint32_t currRow = 2*period + 1  + deepRows[theDeepIndex][0]
                                             - deepRows[theDeepIndex][1]
                                             + 1;
            success(j, pRows, 2*period, currRow - 1);
            free(pRows);
         }
         else{
            success(j, NULL, 0, 0);
         }
         /*
         if (allPhases) {
            k=j;
            for (ph = 1; ph < period; ph++) {
               k=PARENT(k);
               success(k, NULL, 0, 0);
            }
         }
         */
      }
      j--;
   }
}

/* =============================== */
/*  Check parameters for validity  */
/* =============================== */

int loadDumpFlag = 0;
int previewFlag = 0;
int initRowsFlag = 0;

void optError(const char *errorMsg, const char *opt) {
   fprintf(stderr, "Error: %s%s\n", errorMsg, opt);
   aborting = 1;
}

void printError(const char *errorMsg) {
   optError(errorMsg,"");
}

/* Reads and reports values from nttable for a range of conditions.
** input:  a string of birth conditions or a string
**         of survival conditions (not both)
**         (e.g., "B34-w6ci" or "S0123")
**                                                   
** output: the function returns -1 if all input
**         conditions have value -1 in nttable;
**         otherwise, it returns the value in
**         nttable of the input conditions if they
**         are all the same (ignoring -1), and it
**         returns 2 if they are not the same.
*/
int checkConditions(const char *p) {
   int tempTab[256];
   int bs = 0;
   int i;
   if (*p == 's' || *p == 'S') bs = 256;
   p++;
   int val = -1;
   while (*p != '\0'){
      char dig = *p++;
      int negCount = 0;
      if (*p == '\0' || ('0' <= *p && *p <= '8')){
         for (i = 0; i < 256; i++){
            if (rulekeys[i][0] == dig && nttable[bs+i] != -1){
               if (val == -1)
                  val = nttable[bs+i];
               if (val != nttable[bs+i])
                  return 2;
            }
         }
      }
      else if (*p == '-'){
         p++;
         for (i = 0; i < 256; i++)
            tempTab[256] = 0;
         for (; *p != '\0' && !('0' <= *p && *p <= '8'); p++){
            for (i = 1; i < 256; i++)
               if (rulekeys[i][0] == dig && *p != rulekeys[i][1])
                  tempTab[i]++;
            negCount++;
         }
         for (i = 0; i < 256; i++){
            if (tempTab[i] == negCount && nttable[bs+i] != -1){
               if (val == -1)
                  val = nttable[bs+i];
               if (val != nttable[bs+i])
                  return 2;
            }
         }
      }
      else {
         for (; *p != '\0' && !('0' <= *p && *p <= '8'); p++){
            for (i = 0; i < 256; i++){
               if (rulekeys[i][0] == dig && rulekeys[i][1] == *p && nttable[bs+i] != -1){
                  if (val == -1)
                     val = nttable[bs+i];
                  if (val != nttable[bs+i])
                     return 2;
               }
            }
         }
      }
   }
   return val;
}

void checkRule() {
   /* Errors: no meaningful results possible */
   if (checkConditions("B0") == 1)
      printError("rules with B0 are not supported.");
   
   /* Included below are conditions that prevent spaceships from existing:
   ** For proofs of these conditions, see here:
   ** https://conwaylife.com/forums/viewtopic.php?f=11&t=5471
   ** https://ics.uci.edu/~eppstein/ca/glider.c
   **
   **       x <= 0: none of the conditions are satisfied in
   **               the minimum rule
   ** (x+1)%2 == 0: all of the conditions are satisfied in
   **               the maximum rule
   */
   if (checkConditions("B0") == -1){
      printError("any pattern that is not infinite in both dimensions must contain the B0\n       "
                 "neighborhood.");
   }
   if (checkConditions("B1c") == -1){
      printError("spaceships and waves must contain the B1c neighborhood.");
   }
   else if (checkConditions("B1e2a") == -1) {
      printError("spaceships and waves must contain at least one of the B1e or B2a\n       "
                 "neighborhoods.");
   }
   if (checkConditions("B1c") == 1 && checkConditions("B0") == 0){
      printError("patterns in rules with B1c and without B0 expand in all directions.");
   }
   else if (checkConditions("B1e2a") == 1 && checkConditions("B0") == 0) {
      printError("patterns in rules with B1e2a and without B0 expand in all directions.");
   }
   /* The following checks are done only for spaceship searches */
   if ( (params[P_BOUNDARYSYM] == SYM_UNDEF || params[P_SYMMETRY] == SYM_ASYM) ){ 
      if (checkConditions("B012ac3i") <= 0){
         printError("patterns in rules without any of B012ac3i cannot leave their initial\n       "
                    "bounding box.");
      }
      if (checkConditions("B012ae3a") <= 0) {
         printError("patterns in rules without any of B012ae3a cannot leave their initial\n       "
                    "bounding diamond.");
      }
      if (checkConditions("B01245") <=0 && checkConditions("S012345") <= 0) {
         printError("patterns in rules without any of B01245/S012345 cannot move a distance\n       "
                    "of more than one cell outside their initial bounding diamond.");
      }
      if ( checkConditions("B01e2a") <= 0 && 2 * params[P_OFFSET] > params[P_PERIOD]
                                          && params[P_PERIOD] > 0 ){
         printError("orthogonal spaceship speed limit in rules without any of B01e2a is c/2.");
      }
   }
   
   /* Warnings: no spaceships exist, but maybe we can get some interesting wickstretchers */
   if ( (params[P_BOUNDARYSYM] == SYM_UNDEF || params[P_SYMMETRY] == SYM_ASYM) ){
      if (checkConditions("B0") == 0 && (checkConditions("B23")+1)%2 == 0 && (checkConditions("S0")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with all of B23/S0 and without B0,\n"
                         "         because the trailing edge of a pattern cannot die.\n");
      }
      else if (checkConditions("B0") == 0 && checkConditions("B123") >= 1 && (checkConditions("S0123")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with one of B1, B2, or B3, all of S0123,\n"
                         "         and without B0, because the trailing edge of a pattern cannot die.\n");
      }
      if ((checkConditions("S012acek3aijn4a")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with all of S012acek3aijn4a and\n"
                         "         without B0, because patterns cannot shrink.\n");
      }
      if ((checkConditions("S1234-wz5-aqr6ce")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with all of S1234-wz5-aqr6ce and\n"
                         "         without B0, because connected patterns cannot shrink.\n");
      }
      if ((checkConditions("B34")+1)%2 == 0 && (checkConditions("S12345")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with all of B34/S12345 and without B0,\n"
                         "         because connected patterns cannot shrink.\n");
      }
      if ((checkConditions("B345")+1)%2 == 0 && (checkConditions("S1234")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with all of B345/S1234 and without B0,\n"
                         "         because connected patterns cannot shrink.\n");
      }
      if (checkConditions("B012") <= 0 && (checkConditions("S234567")+1)%2 == 0){
         fprintf(stderr, "Warning: no spaceships exist in rules with all of S234567 and none of B012,\n"
                         "         because patterns cannot escape their bounding diamond without an\n"
                         "         immortal triangle.\n");
      }
   }
}

void checkGutter() {
   int i = 0;
   /* if there are forbidden birth conditions, i will be less than 256 */
   while (i < 256 && nttable[i] != -1)
      i++;
   
   if (checkConditions("B2ci4ci6i") <= 0)
      gutterSkew = 0;
   else if (checkConditions("B1c2kn3ny4yz5r6i") <= 0)
      gutterSkew = 1;
   else if (checkConditions("B12aikn3cqr4cnyz5er6i") <= 0)
      gutterSkew = 2;
   else {
      printError("gutters do not work with the given birth conditions.\n       "
                 "The forbidden birth conditions for different gutter types are\n       "
                 "  Skew 0: B2ce4ci6i\n       "
                 "  Skew 1: B1c2kn3ny4yz5r6i\n       "
                 "  Skew 2: B12aikn3cqr4cnyz5er6i");
   }
   if (gutterSkew && i < 256)
      fprintf(stderr, "Warning: forbidden birth conditions cannot be checked along a skew gutter.\n");
}

void checkParams() {
   const char *ruleError;
   
   /* Errors */
   
   /* There would probably be several integer overflow bugs if sizeof(int) == 2. */
   if (sizeof(int) == 2)
      printError("This program does not work when compiled in 16-bit mode.\n       "
                 "Please recompile.");
   
   ruleError = parseRule(rule, nttable);
   
   if (ruleError != 0){
      optError("failed to parse rule ", rule);
      fprintf(stderr, "       %s\n", ruleError);
   }
   else
      checkRule();
   
   if (params[P_SYMMETRY] == SYM_GUTTER || params[P_BOUNDARYSYM] == SYM_GUTTER)
      checkGutter();
   
#ifdef QSIMPLE
   if (gcd(PERIOD,OFFSET) > 1)
      printError("qfind-s does not support gcd(PERIOD,OFFSET) > 1. Use qfind instead.");
#else
   if (params[P_PERIOD] > MAXPERIOD)
      printError("maximum allowed period (" XSTR(MAXPERIOD) ") exceeded.");
   if (params[P_OFFSET] > params[P_PERIOD] && params[P_PERIOD] > 0)
      printError("translation cannot exceed period.");
   if (params[P_OFFSET] == params[P_PERIOD] && params[P_PERIOD] > 0)
      printError("photon searches are not supported.");
#endif
   if (params[P_PERIOD] == 0)
      printError("you must specify a velocity (-v).");
   if (params[P_WIDTH] == 0)
      printError("you must specify a width (-w).");
   if (params[P_SYMMETRY] == SYM_UNDEF)
      printError("you must specify a symmetry type (-s).");
   if (params[P_BOUNDARYSYM] == SYM_ASYM)
      printError("asymmetric wave searching is not supported.");
   if (previewFlag && !loadDumpFlag)
      printError("the search state must be loaded from a file to preview partial results.\n");
   if (initRowsFlag && loadDumpFlag){
      printError("initial rows file cannot be used when the search state is loaded from a\n       "
                 "saved state.");
   }
   if (params[P_QBITS] <= 0)
      printError("queue bits (-q) must be positive.");
   if (params[P_BASEBITS] <= 0)
      printError("base bits (-b) must be positive.");
   if (params[P_BASEBITS] >= params[P_QBITS])
      printError("base bits (-b) must be less than queue bits (-q).");
   if (params[P_HASHBITS] < 0)
      printError("hash bits (-h) must be nonnegative.");
   
   /* Warnings */
   if (2 * params[P_OFFSET] > params[P_PERIOD] && params[P_PERIOD] > 0){
      fprintf(stderr, "Warning: searches for speeds exceeding c/2 may not work correctly.\n");
   }
#ifdef NOCACHE
   if (5 * params[P_OFFSET] > params[P_PERIOD] && params[P_PERIOD] > 0 && params[P_CACHEMEM] == 0){
      fprintf(stderr, "Warning: Searches for speeds exceeding c/5 may be slower without caching.\n"
                      "         It is recommended that you increase the cache memory (-c).\n");
   }
#else
   if (5 * params[P_OFFSET] <= params[P_PERIOD] && params[P_OFFSET] > 0 && params[P_CACHEMEM] > 0){
      fprintf(stderr, "Warning: Searches for speeds at or below c/5 may be slower with caching.\n"
                      "         It is recommended that you disable caching (-c 0).\n");
   }
#endif
   if (params[P_SYMMETRY] == SYM_ASYM && params[P_BOUNDARYSYM] != SYM_UNDEF){
      fprintf(stderr, "Warning: the wave symmetry settings are equivalent to a spaceship search.\n");
      params[P_SYMMETRY] = params[P_BOUNDARYSYM];
      params[P_BOUNDARYSYM] = SYM_UNDEF;
   }
   /* Reduce values to prevent integer overflow */
   if (params[P_QBITS] > 31 && !aborting){
      fprintf(stderr, "Warning: queue bits (-q) reduced to 31.\n");
      params[P_QBITS] = 31;      /* corresponds to a queue size of 2GB */
      if (params[P_BASEBITS] > params[P_QBITS]){
         fprintf(stderr, "Warning: base bits (-q) reduced to 30.\n");
         params[P_BASEBITS] = 30;
      }
   }
   if (params[P_HASHBITS] > 31 && !aborting){
      fprintf(stderr, "Warning: hash bits (-h) reduced to 31.\n");
      params[P_HASHBITS] = 31;   /* corresponds to a hash table size of 2GB */
   }
   
}

/* ============================ */
/*  Load saved state from file  */
/* ============================ */

int splitNum = 0;
char * loadFile;

void loadFail() {
   fprintf(stderr, "Load from file %s failed\n", loadFile);
   exit(1);
}

signed int loadInt(FILE *fp) {
   signed int v;
   if (fscanf(fp,"%d\n",&v) != 1) loadFail();
   return v;
}

unsigned long loadUInt(FILE *fp) {
   unsigned long v;
   if (fscanf(fp,"%lu\n",&v) != 1) loadFail();
   return v;
}

void loadParams() {
   FILE * fp;
   unsigned int i;
   
   fp = fopen(loadFile, "r");
   if (!fp) loadFail();
   if (loadUInt(fp) != FILEVERSION){
      printf("Incompatible file version\n");
      exit(1);
   }
   
   /* Load rule */
   if (fscanf(fp,"%150s\n",loadRule) != 1) loadFail();
   rule = loadRule;
   
   /* Load dump root */
   if (fscanf(fp,"%250s\n",loadDumpRoot) != 1) loadFail();
   dumpRoot = loadDumpRoot;
   
   /* Load parameters */
   for (i = 0; i < NUM_PARAMS; ++i)
      params[i] = loadInt(fp);
}

void loadState() {
   FILE * fp;
   unsigned long i, j;
   int k;
   
   fp = fopen(loadFile, "r");
   if (!fp) loadFail();
   
   /* Skip lines that are loaded in loadParams() */
   loadUInt(fp);                                   /* skip file version */
   if (fscanf(fp, "%*[^\n]\n") != 0) loadFail();   /* skip rule */
   if (fscanf(fp, "%*[^\n]\n") != 0) loadFail();   /* skip dump root */
   for (i = 0; i < NUM_PARAMS; ++i) loadInt(fp);   /* skip parameters */
   
   /* Load / initialise globals */
   width          = loadInt(fp);
   period         = loadInt(fp);
   offset         = loadInt(fp);
   lastDeep       = loadInt(fp);
   dumpNum        = loadInt(fp);
   if (params[P_DUMPMODE] == D_SEQUENTIAL)
      dumpNum = 1;
   
   
   aborting        = 0;
   nRowsInState    = period+period;   /* how many rows needed to compute successor graph? */
   
   params[P_DEPTHLIMIT] = DEFAULT_DEPTHLIMIT;
   
    /* Allocate space for the data structures */
   base = (node*)malloc((QSIZE>>BASEBITS)*sizeof(node));
   rows = (row*)malloc(QSIZE*sizeof(row));
   if (base == 0 || rows == 0) {
      printf("Unable to allocate BFS queue!\n");
      exit(1);
   }
   
   if (hashBits == 0) hash = 0;
   else {
      hash = (node*)malloc(HASHSIZE*sizeof(node));
      if (hash == 0) printf("Unable to allocate hash table, duplicate elimination disabled\n");
   }
   
   /* Load up BFS queue */
   qHead  = (node) loadUInt(fp);
   qEnd   = (node) loadUInt(fp);
   qStart = QSIZE - qEnd;
   qEnd   = QSIZE;
   qHead += qStart;
   if (qStart > QSIZE || qStart < QSIZE/16) {
      printf("BFS queue is too small for saved state\n");
      exit(1);
   }
   for (i = qStart; i < qEnd; ++i)
      rows[i] = (row) loadUInt(fp);
   
   /* Load extension rows for each queue node */
   deepRows = (row**)calloc(1LLU << (params[P_DEPTHLIMIT] + 1),sizeof(*deepRows));
   deepRowIndices = (uint32_t*)calloc(QSIZE,sizeof(deepRowIndices));
   
   uint32_t theDeepIndex = 2;
   deepQTail = 0;
   
   while ( (k = fscanf(fp,"%lu\n",&j)) != EOF){
      if (k == 0)
         loadFail();
      if (j == 0){
         j = loadUInt(fp);
         for (i = 0; i < j; ++i){
            deepRowIndices[deepQTail] = 1;
            ++deepQTail;
         }
         continue;
      }
      deepRows[theDeepIndex] = (row*)calloc( j + 1 + 2, sizeof(**deepRows));
      deepRows[theDeepIndex][0] = j;
      for (i = 1; i < j + 1 + 2; ++i){
         deepRows[theDeepIndex][i] = (row) loadUInt(fp);
      }
      deepRowIndices[deepQTail] = theDeepIndex;
      ++theDeepIndex;
      ++deepQTail;
   }
   
   fclose(fp);
   
   /* complete compaction */
   doCompactPart2();
   
   /* Let the user know that we got this far (suppress if splitting) */
   if (!splitNum) printf("State successfully loaded from file %s\n",loadFile);
   
   fflush(stdout);
}

/* ================================================= */
/*  Load initial rows for extending partial results  */
/* ================================================= */

void printRow(row theRow) {
   int i;
   for (i = width - 1; i >= 0; --i) printf("%c",(theRow & 1 << i ? 'o' : '.'));
   printf("\n");
}

void loadInitRows(char * file) {
   FILE * fp;
   int i,j;
   char rowStr[MAXWIDTH];
   row theRow = 0;
   
   loadFile = file;
   fp = fopen(loadFile, "r");
   if (!fp) loadFail();
   
   printf("Starting search from rows in %s:\n",loadFile);
   
   for (i = 0; i < 2 * period; i++){
      if (fscanf(fp,"%s",rowStr) != 1) loadFail();
      for (j = 0; j < width; j++){
         theRow |= ((rowStr[width - j - 1] == '.') ? 0:1) << j;
      }
      printRow(theRow);
      enqueue(dequeue(),theRow);
      theRow = 0;
   }
   fclose(fp);
}

/* ============== */
/*  Set Defaults  */
/* ============== */

void setDefaultParams() {
#ifdef QSIMPLE
   params[P_PERIOD] = PERIOD;
   params[P_OFFSET] = OFFSET;
#else
   params[P_PERIOD] = 0;
   params[P_OFFSET] = 0;
#endif
   params[P_WIDTH] = 0;
   params[P_SYMMETRY] = SYM_UNDEF;
   params[P_REORDER] = 1;     /* 0 and 2 are also valid values, but this cannot currently be set at runtime. */
   params[P_DUMPINTERVAL] = 1800;    /* 30 minutes */
   params[P_BASEBITS] = 4;
   params[P_QBITS] = QBITS;
   params[P_HASHBITS] = HASHBITS;
   params[P_NUMTHREADS] = 1;
   params[P_MINDEEP] = 3;
   /* A negative value for params[P_CACHEMEM] means use that amount of  */
   /* memory if speed > c/5 and turn off caching otherwise.  A positive */
   /* value forces caching even if speed <= c/5.                        */
   params[P_CACHEMEM] = -1*DEFAULT_CACHEMEM;
   params[P_MEMLIMIT] = -1;   /* negative value means no limit */
   params[P_PRINTDEEP] = 1;
   params[P_LONGEST] = 1;
   params[P_FIRSTDEEP] = 0;
   params[P_NUMSHIPS] = 0;
   params[P_MINEXTENSION] = 0;
   params[P_FULLPERIOD] = 0;
   params[P_BOUNDARYSYM] = SYM_UNDEF;
   params[P_DUMPMODE] = D_OVERWRITE;
   params[P_EVERYDEPTH] = 0;
   params[P_EARLYEXIT] = 1;
}

/* =============== */
/*  Parse options  */
/* =============== */

#define no_argument       (0)
#define required_argument (1)
#define optional_argument (2)

struct option {
   const char *name;
   int has_arg;
   int val;
};

/* Nonstandard getopt_long():
** This is a simple implementation of a getopt-like function.
** It's missing some features that ordinary getopt_long() has,
** but it has mostly the same input and output.
**
** Differences from getopt_long():
**   I/O:   char **optName and char **optArg:
**             used to save pointers to option names and arguments.
**             *optName and *optArg will point to elements of argv[].
**             If no argument is provided, *optArg will be 0.
**             If my_getopt() returns -1, then *optName and *optArg
**             will be unchanged.
*/
int my_getopt( int argc,
               char *argv[],
               const char *shortOpts,
               struct option *longOpts,
               char **optName,
               char **optArg   )
{
   static int i = 0;
   const char *charInd;
   if (++i == argc)  /* we're at the end of argv */
      return -1;
   *optName = argv[i];
   *optArg  = NULL;
   if (argv[i][0] == '-'){
      if (argv[i][1] == '\0')
         return '?';
      /* check if short option is contained in short options list */
      else if (argv[i][1] != '-' && argv[i][2] == '\0' && (charInd = strchr(shortOpts, argv[i][1]))){
         if (*(charInd+1) == ':'){
            if (argc == i+1){  /* no argument */
               if (*(charInd+2) != ':')
                  return (shortOpts[0] == ':' ? ':' : '?');
            }
            else{
               *optArg = argv[++i];
               return (int)argv[i-1][1];
            }
         }
         return (int)argv[i][1];
      }
      /* check if long option is contained in long options list */
      else if (argv[i][1] == '-') {
         while (longOpts->name != 0 && strcmp(argv[i] + 2, longOpts->name)) 
            longOpts++;
         if (longOpts->name){
            if (longOpts->has_arg){
               if (argc == i+1){  /* no argument */
                  if (longOpts->has_arg != optional_argument)
                     return (shortOpts[0] == ':' ? ':' : '?');
               }
               else
                  *optArg = argv[++i];
            }
            return longOpts->val;
         }
      }
   }
   return '?';
}

int readInt(char *opt, char *arg) {
   int i = 0;
   char c;
   if (arg == 0)
      aborting = 1;
   else if (sscanf(arg, "%d%c", &i, &c) != 1){
      fprintf(stderr, "Error: invalid argument %s in option %s.\n", arg, opt);
      aborting = 1;
   }
   return i;
}

const char *parseVelocity(char *velString, int *per, int *off) {
   int xoff = 0;
   *per = 1;
   *off = 1;
   char b = 0;
   char c = 0;
   if (!strcmp(velString,"c"))
      return 0;
   else if (sscanf(velString, "c/%d%c%c", per, &b, &c) >= 1){
      if (b == 'd' && !c)
         return "diagonal spaceship searches are not supported.";
      else if (b == '\0' || (b == 'o' && !c))
         return 0;
      else
         return "illegal characters after velocity";
   }
   else if (sscanf(velString, "%dc/%d%c%c", off, per, &b, &c) >= 2){
      if (*off == 0)
         return "oscillator searches are not supported.";
      else if (b == 'd' && !c)
         return "diagonal spaceship searches are not supported.";
      else if (*off < 0)
         return "offset must be positive.";
      else if (b == '\0' || (b == 'o' && !c))
         return 0;
      else
         return "illegal characters after velocity";
   }
   else if (sscanf(velString, "(%d,%d)c/%d%c", off, &xoff, per, &c) >= 3){
      if (c != '\0' && c != 'o')
         return "illegal characters after velocity";
      else if (xoff != 0) {
         if (*off == 0){
            *off = xoff;
            xoff = 0;
            return 0;
         }
         else if (xoff == *off || xoff * (-1) == *off)
            return "diagonal spaceship searches are not supported.";
         else
            return "oblique spaceship searches are not supported.";
      }
      else if (*off == 0)
         return "oscillator searches are not supported.";
      else if (*off < 0)
         return "offset must be positive.";
      else
         return 0;
   }
   return "Unable to read offset and period.";
}

void parseOptions(int argc, char *argv[]) {
   char *optArg = 0;
   char *optName = 0;
   int c;
   
   if (argc <= 1){
      printf("\n");
      printHelp();
   }
   printf("Input:");
   for (c = 1; c < argc; c++)
      printf(" %s", argv[c]);
   printf("\n\n");
   
   /* list of long options */
   struct option options[] = {
   // {"",                    no_argument,        -1},   /* end option parsing when "--" encountered */
      {"help",                no_argument,       256},
      {"rule",                required_argument, 'r'},
      {"width",               required_argument, 'w'},
      {"symmetry",            required_argument, 's'},
      {"boundary-sym",        required_argument, 'o'},
      {"boundary-symmetry",   required_argument, 'o'},
      {"mem-limit",           required_argument, 'm'},
      {"memory-limit",        required_argument, 'm'},
      {"cache-mem",           required_argument, 'c'},
      {"cache-memory",        required_argument, 'c'},
      {"first-depth",         required_argument, 'n'},
      {"increment",           required_argument, 'i'},
      {"queue-bits",          required_argument, 'q'},
      {"hash-bits",           required_argument, 'h'},
      {"base-bits",           required_argument, 'b'},
#ifdef _OPENMP
      {"threads",             required_argument, 't'},
#endif
      {"found",               required_argument, 'f'},
      {"min-extension",       required_argument, 'g'},
      {"minimum-extension",   required_argument, 'g'},
      {"extend",              required_argument, 'e'},
      {"dump-root",           required_argument, 'd'},
      {"load",                required_argument, 'l'},
      {"split",               required_argument, 'j'},
      {"dump-interval",       required_argument, 'a'},
      {"dump-int",            required_argument, 'a'},
      {"preview",             no_argument,       'p'},
#ifndef QSIMPLE
      {"velocity",            required_argument, 'v'},
      {"enable-subperiod",    no_argument,       257},
      {"enable-subperiodic",  no_argument,       257},
      {"disable-subperiod",   no_argument,       258},
      {"disable-subperiodic", no_argument,       258},
#endif
      {"enable-deep-print",   no_argument,       259},
      {"disable-deep-print",  no_argument,       260},
      {"enable-longest",      no_argument,       261},
      {"disable-longest",     no_argument,       262},
      {"dump-mode",           required_argument, 263},
      {"fixed-depth",         required_argument, 264},
#ifdef _OPENMP
      {"enable-early-exit",   no_argument,       265},
      {"disable-early-exit",  no_argument,       266},
#endif
      {0, 0, 0}   /* marks end of long options list */
   };
   
   while ( (c = my_getopt( argc,
                           argv,
                           ":"   /* List of short options; "<option>:" means argument required */
#ifndef QSIMPLE
                           "kKv:V:"
#endif
#ifdef _OPENMP
                           "t:T:"
#endif
                           "a:b:c:d:e:f:g:h:i:j:l:m:n:o:pq:r:s:w:z"     /* Currently unused: */
                           "A:B:C:D:E:F:G:H:I:J:L:M:N:O:PQ:R:S:W:Z",    /* u,x,y,U,X,Y       */
                           options,
                           &optName,
                           &optArg)) != -1 )
   {
      switch (c) {
         case 'r': case 'R':
            rule = optArg;
            if (strlen(rule) > 150)  /* any valid rule can be written in 141 characters or fewer */
               printError("rule string exceeds maximum allowed length (150).\n       "
                          "You must write the rule more efficiently.\n");
            break;
#ifndef QSIMPLE
         case 'v': case 'V':
         {
            const char *velError = 0;
            velError = parseVelocity(optArg, &params[P_PERIOD], &params[P_OFFSET]);
            if (velError){
               optError("invalid velocity ", optArg);
               fprintf(stderr, "       %s\n", velError);
               /* prevent additional error message in checkParams() */
               params[P_PERIOD] = 2; params[P_OFFSET] = 1;
            }
            if (params[P_PERIOD] <= 0){
               optError("invalid velocity ", optArg);
               fprintf(stderr, "       Period must be positive\n");
               /* prevent additional error message in checkParams() */
               params[P_PERIOD] = 2; params[P_OFFSET] = 1;
            }
            break;
         }
         case 'k': case 'K':
            params[P_FULLPERIOD] ^= 1;
            break;
         case 257:   /* --enable-subperiod */
            params[P_FULLPERIOD] = 0;
            break;
         case 258:   /* --disable-subperiod */
            params[P_FULLPERIOD] = 1;
            break;
#endif
         case 'w': case 'W':
            params[P_WIDTH] = readInt(optName, optArg);
            if (params[P_WIDTH] <= 0){
               printError("width must be positive");
               params[P_WIDTH] = 1;    /* prevent additional error message in checkParams() */
            }
            break;
         case 's': case 'S':
            switch(optArg[0]) {
               case 'a': case 'A':
                  params[P_SYMMETRY] = SYM_ASYM; break;
               case 'o': case 'O':
                  params[P_SYMMETRY] = SYM_ODD; break;
               case 'e': case 'E':
                  params[P_SYMMETRY] = SYM_EVEN; break;
               case 'g': case 'G':
                  params[P_SYMMETRY] = SYM_GUTTER; break;
               default:
                  optError("unrecognized symmetry type ", optArg);
                  break;
            }
            break;
         case 'o': case 'O':
            switch(optArg[0]) {
               case 'a': case 'A':
                  params[P_BOUNDARYSYM] = SYM_ASYM; break;
               case 'o': case 'O':
                  params[P_BOUNDARYSYM] = SYM_ODD; break;
               case 'e': case 'E':
                  params[P_BOUNDARYSYM] = SYM_EVEN; break;
               case 'g': case 'G':
                  params[P_BOUNDARYSYM] = SYM_GUTTER; break;
               case 'd': case 'D':
                  params[P_BOUNDARYSYM] = SYM_UNDEF; break;
               default:
                  optError("unrecognized symmetry type ", optArg);
                  break;
            }
            break;
         case 'm': case 'M':
            params[P_MEMLIMIT] = readInt(optName, optArg);
            break;
         case 'n': case 'N':
            params[P_FIRSTDEEP] = readInt(optName, optArg);
            if (params[P_FIRSTDEEP] <= 0)
               printError("first depth must be positive.");
            break;
         case 'c': case 'C':
            params[P_CACHEMEM] = readInt(optName, optArg);
            break;
         case 'i': case 'I':
            params[P_MINDEEP] = readInt(optName, optArg);
            break;
         case 'q': case 'Q':
            params[P_QBITS] = readInt(optName, optArg);
            break;
         case 'h': case 'H':
            params[P_HASHBITS] = readInt(optName, optArg);
            break;
         case 'b': case 'B':
            params[P_BASEBITS] = readInt(optName, optArg);
            break;
#ifdef _OPENMP
         case 't': case 'T':
            params[P_NUMTHREADS] = readInt(optName, optArg);
            break;
#endif
         case 'f': case 'F':
            params[P_NUMSHIPS] = readInt(optName, optArg);
            break;
         case 'g': case 'G':
            params[P_MINEXTENSION] = readInt(optName, optArg);
            break;
         case 'z': case 'Z':
            params[P_PRINTDEEP] ^= 1;
            break;
         case 'p': case 'P':
            previewFlag = 1;
            break;
         case 'd': case 'D':
            dumpRoot = optArg;
            if (strlen(dumpRoot) > MAXDUMPROOT)
               printError("dump root exceeds maximum allowed length (" XSTR(MAXDUMPROOT) ")");
            break;
         case 'j': case 'J':
            splitNum = readInt(optName, optArg);
            if (splitNum < 0) splitNum = 0;
            break;
         case 'e': case 'E':
            initRows = optArg;
            initRowsFlag = 1;
            break;
         case 'l': case 'L':
            loadFile = optArg;
            loadDumpFlag = 1;
            loadParams();
            break;
         case 'a': case 'A':
            params[P_DUMPINTERVAL] = readInt(optName, optArg);
            if (params[P_DUMPINTERVAL] < 0)
               printError("dump interval must be nonnegative");
            break;
         case 259:   /* --enable-deep-printing */
            params[P_PRINTDEEP] = 1;
            break;
         case 260:   /* --disable-deep-printing */
            params[P_PRINTDEEP] = 0;
            break;
         case 261:   /* --enable-longest-partial */
            params[P_LONGEST] = 1;
            break;
         case 262:   /* --disable-longest-partial */
            params[P_LONGEST] = 0;
            break;
         case 263:   /* --dump-mode */
            switch(optArg[0]) {
               case 'o': case 'O':
                  params[P_DUMPMODE] = D_OVERWRITE; break;
               case 's': case 'S':
                  params[P_DUMPMODE] = D_SEQUENTIAL; break;
               case 'd': case 'D':
                  params[P_DUMPMODE] = D_DISABLED; break;
               default:
                  optError("unrecognized dump mode ", optArg);
                  break;
            }
            break;
         case 264:   /* --fixed-depth */
            params[P_EVERYDEPTH] = 1;
            params[P_MINDEEP] = 1;
            params[P_FIRSTDEEP] = readInt(optName, optArg);
            if (params[P_FIRSTDEEP] <= 0)
               printError("fixed depth must be positive.");
            break;
         case 265:   /* --enable-early-exit */
            params[P_EARLYEXIT] = 1;
            break;
         case 266:   /* --disable-early-exit */
            params[P_EARLYEXIT] = 0;
            break;
         case 256:   /* --help */
            printHelp();
            break;
         case ':':
            optError("missing argument for option ", optName);
            break;
         case '?':
            optError("unrecognized option ", optName);
            break;
         default:
            printError("option parser failed");
            break;
      }
   }
}

/* ========================================= */
/*  Set up search with the given parameters  */
/* ========================================= */

void searchSetup() {
   if (params[P_CACHEMEM] < 0){
      if (5 * params[P_OFFSET] > params[P_PERIOD]) params[P_CACHEMEM] *= -1;
      else params[P_CACHEMEM] = 0;
   }
   
   checkParams();  /* Exit if parameters are invalid */
   
   if (aborting){
      fprintf(stderr, "\nUse --help for a list of available options.\n");
      exit(1);
   }
   
   if (loadDumpFlag) loadState();
   else {
      width = params[P_WIDTH];
      period = params[P_PERIOD];
      offset = params[P_OFFSET];
      hashPhase = (gcd(period,offset)>1);
      
      nRowsInState = period+period;
      
      params[P_DEPTHLIMIT] = DEFAULT_DEPTHLIMIT;
      
      base = (node*)malloc((QSIZE>>BASEBITS)*sizeof(node));
      rows = (row*)malloc(QSIZE*sizeof(row));
      if (base == 0 || rows == 0) {
         printf("Unable to allocate BFS queue!\n");
         exit(1);
      }
      
      if (hashBits == 0) hash = 0;
      else {
         hash = (node*)malloc(HASHSIZE*sizeof(node));
         if (hash == 0) printf("Unable to allocate hash table, duplicate elimination disabled\n");
      }
      
      deepRows = (row**)calloc(1LLU << (params[P_DEPTHLIMIT] + 1),sizeof(*deepRows));
      deepRowIndices = (uint32_t*)calloc(QSIZE,sizeof(deepRowIndices));
      
      resetQ();
      resetHash();
      
      enqueue(0,0);
      
      if (initRowsFlag) loadInitRows(initRows);
   }
   
#ifndef QSIMPLE
   makePhases();
   makeSubperiodTables();
#endif
   
   /* Generate proper rule string for printing patterns */
   int i;
   int j = 0;
   int k = 1;
   for (i = 0; i < 151 && rule[i] != '\0'; i++){
      if (rule[i] == '~') k = 0;
      else if (rule[i] == '/') k = 1;
      if (k) baseRule[j++] = rule[i];
   }
   baseRule[j] = '\0';
   
   if (previewFlag){
      params[P_NUMSHIPS] = 0;
      preview();
      exit(0);
   }
   
   if (params[P_MINEXTENSION] < 0) params[P_MINEXTENSION] = 0;
   if (params[P_FIRSTDEEP] < 0) params[P_FIRSTDEEP] = 0;
   
   dumpMode = params[P_DUMPMODE];   /* these may need to be different when splitting */
   
   /* split queue across multiple files */
   if (splitNum > 0){
      node x;
      uint32_t i,j;
      uint32_t deepIndex;
      int firstDumpNum = 0;
      int totalNodes = 0;
      
      dumpMode = D_SEQUENTIAL; 
      
      echoParams();
      printf("\n");
      
      if (!loadDumpFlag || qHead == 0 || splitNum == 1){
         dumpFlag = DUMPPENDING;
         if (qHead == 0){      /* can't use doCompact() here, because it tries to access rows[-1] */
            qStart = qHead;
            qEnd = qTail;
            dumpState();
         }
         else doCompact();
         if (dumpFlag == DUMPSUCCESS){
            printf("State dumped to %s\n",dumpFile);
            exit(0);
         }
         else{
            fprintf(stderr, "Error: dump failed.\n");
            exit(1);
         }
      }
      
      if (splitNum >= 100000){
         fprintf(stderr, "Warning: queue cannot be split into more than 99999 files.\n");
         splitNum = 99999;
      }
      
      /* count nodes in queue */
      for (x = qHead; x < qTail; x++){
         if (!EMPTY(x)) totalNodes++;
      }
      
      /* nodes per file is rounded up */
      unsigned long nodesPerFile = (totalNodes - 1) / splitNum + 1;
      
      printf("Splitting search state with %lu queue nodes per file\n",nodesPerFile);
      
      /* save qHead and qTail, as creating the pieces will change their values */
      node fixedQHead = qHead;
      node fixedQTail = qTail;
      
      /* delete the queue; we will reload it as needed */
      free(base);
      free(rows);
      free(hash);
      
      for (deepIndex = 0; deepIndex < 1LLU << (params[P_DEPTHLIMIT] + 1); ++deepIndex){
         if (deepRows[deepIndex]) free(deepRows[deepIndex]);
         deepRows[deepIndex] = 0;
      }
      free(deepRows);
      free(deepRowIndices);
      
      node currNode = fixedQHead;
      
      while (currNode < fixedQTail){
         
         /* load the queue */
         loadState();
         
         /* empty everything before currNode */
         j = deepQHead;
         for (x = fixedQHead; x < currNode; ++x){
            MAKEEMPTY(x);
            deepRowIndices[j] = 0;
            ++j;
         }
         
         /* skip the specified number of nonempty nodes */
         i = 0;
         while (i < nodesPerFile && x < fixedQTail){
            if (!EMPTY(x))
               ++i;
            ++x;
            ++j;
         }
         
         /* update currNode */
         currNode = x;
         
         /* empty everything after specified nonempty nodes */
         while (x < fixedQTail){
            MAKEEMPTY(x);
            deepRowIndices[j] = 0;
            ++x;
            ++j;
         }
         
         /* save the piece */
         dumpFlag = DUMPPENDING;
         doCompact();
         
         if (!firstDumpNum) firstDumpNum = dumpNum - 1;
         
         if (dumpFlag != DUMPSUCCESS){
            printf("Failed to save %s\n",dumpFile);
            exit(1);
         }
         
         if (dumpNum >= DUMPLIMIT){
            fprintf(stderr, "Error: dump file number limit (" XSTR(DUMPLIMIT) ") reached.\n");
            fprintf(stderr, "       Try splitting the queue in a new directory.\n");
            exit(1);
         }
         
         /* free memory allocated in loadState() */
         for (deepIndex = 0; deepIndex < 1LLU << (params[P_DEPTHLIMIT] + 1); ++deepIndex){
            if (deepRows[deepIndex]) free(deepRows[deepIndex]);
            deepRows[deepIndex] = 0;
         }
         
         free(base);
         free(rows);
         free(hash);
         free(deepRows);
         free(deepRowIndices);
      }
      
      printf("Saved pieces in files %s%05d to %s\n",dumpRoot,firstDumpNum,dumpFile);
      exit(0);
   }
   
   omp_set_num_threads(params[P_NUMTHREADS]);
   
   memlimit = ((long long)params[P_MEMLIMIT]) << 20;
   
   /* Allocate lookahead cache */
#ifndef NOCACHE
   cachesize = 32768;
   while (cachesize * sizeof(cacheentry) < 550000 * (unsigned long long)params[P_CACHEMEM])
      cachesize <<= 1;
   memusage += sizeof(cacheentry) * (cachesize + 5) * params[P_NUMTHREADS];
   if (params[P_MEMLIMIT] >= 0 && memusage > memlimit){
      printf("Not enough memory to allocate lookahead cache\n");
      exit(1);
   }
   totalCache = (cacheentry *)calloc(sizeof(cacheentry),
         (cachesize + 5) * params[P_NUMTHREADS]);
   cache = (cacheentry **)calloc(sizeof(**cache), params[P_NUMTHREADS]);
   
   for (int i = 0; i < params[P_NUMTHREADS]; i++)
      cache[i] = totalCache + (cachesize + 5) * i;
#endif
   
   echoParams();
   
   fasterTable();
   makeTables();
   
   rephase();
   
   parseDumpRoot();
   time(&lastDumpTime);
   
   timeStamp();
}

void finalReport() {
   timeStamp();
   printf("Search complete.\n\n");
   
   printf("%d %s%s found.\n",numFound,(params[P_BOUNDARYSYM] == SYM_UNDEF) ? "spaceship" : "wave",(numFound == 1) ? "" : "s");
   printf("Maximum depth reached: %d\n",longest);
   if (params[P_LONGEST] && aborting != 3){ /* aborting == 3 means we reached ship limit */
      if (patternBuf) printf("Longest partial result:\n\n%s",patternBuf);
      else printf("No partial results found.\n");
   }
}
