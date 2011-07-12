/* **************************frng************************** B Bunk 8/2010 
                                                            rev   12/2010
   Multiplicative lagged Fibonacci random number generator

   Implementation with 64 bits and recursion
      x(n) = x(n-13) * x(n-31) mod 2^64

   Entries: 
   - frng()    float
   - dfrng()   double
   - ifrng()   int   (31 most significant bits)
   - lfrng()   INT64 (63 most significant bits)
   - calls which fill a random vector:
      frngv(n, rvec), dfrngv(n, dvec), ifrngv(n, ivec), lfrngv(n, lvec)
   - frnget(iseed), frnset(iseed) : get/set register
            INT64 iseed[31]      (63 most significant bits)
   - frnini(num) : initialise register (num = 0 .. 714024)

   This program requires 64bit integers, named INT64, which can be either
   `long' or `long long'.
   The size of `long' integers varies e.g. between 32bit and 64bit Linux.
   Integers of type `long long' are 64bit in both cases, but require C99
   (instead of Ansi C / C90).
   If C99 is to be avoided (e.g. in C++), but `long' is 64bit anyway,
   compile with -DLONG64 .
*/
#include <stdlib.h>           /* for labs() */
#include "frng.h"             /* check prototypes */

#ifdef LONG64
#define INT64 long
#else
#define INT64 long long
#endif

/* (LEN,LAG) = (17,5), (31,13), (55,24), ... */
#define LEN 31
#define LAG 13

static unsigned INT64 ireg[LEN] =
#ifdef LONG64
  { 7569555838333416573LU, 3279508493814418161LU,14601197116916481131LU
  ,11356974270639151889LU,11533207906012338279LU,   48396622165758097LU
  , 9799506170800265737LU,10495884947063469839LU,15762222231453528583LU
  ,  592708693962015337LU,  249471268379597569LU,18293070441671520277LU
  , 1095135063400818497LU,12274476774300682419LU,12956519111941717915LU
  , 6258116919374872021LU, 7125927167643344901LU, 3558223830886234463LU
  , 6629248264238368005LU, 9006395648747532825LU, 9815540178537739879LU
  , 7397323557512132751LU,16925591859524868917LU,11024015203834196785LU
  ,11007626191176900455LU, 7857178218845272995LU,  406276102641255367LU
  , 6196339901492533131LU, 7804281187626479095LU,15539499534985592797LU
  , 4020080149601329621LU};
#else
  { 7569555838333416573LLU, 3279508493814418161LLU,14601197116916481131LLU
  ,11356974270639151889LLU,11533207906012338279LLU,   48396622165758097LLU
  , 9799506170800265737LLU,10495884947063469839LLU,15762222231453528583LLU
  ,  592708693962015337LLU,  249471268379597569LLU,18293070441671520277LLU
  , 1095135063400818497LLU,12274476774300682419LLU,12956519111941717915LLU
  , 6258116919374872021LLU, 7125927167643344901LLU, 3558223830886234463LLU
  , 6629248264238368005LLU, 9006395648747532825LLU, 9815540178537739879LLU
  , 7397323557512132751LLU,16925591859524868917LLU,11024015203834196785LLU
  ,11007626191176900455LLU, 7857178218845272995LLU,  406276102641255367LLU
  , 6196339901492533131LLU, 7804281187626479095LLU,15539499534985592797LLU
  , 4020080149601329621LLU};
#endif

static int i0 = 0, i1 = LEN - LAG;
#pragma omp threadprivate(ireg, i0, i1)

static int iup[LEN] =
  {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16
  , 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,  0};

static double dnorm63 = 1.08420217248550443400745280086994171142578125e-19;
                                                         /* 2^(-63) */
static float  rnorm31 = 4.656612873077392578125e-10;     /* 2^(-31) */

float frng(){
  float retval;
  ireg[i0] = ireg[i0] * ireg[i1];
  retval = rnorm31 * (ireg[i0] >> 33);
  i0 = iup[i0];
  i1 = iup[i1];
  return retval;
}

double dfrng(){
  double retval;
  ireg[i0] = ireg[i0] * ireg[i1];
  retval = dnorm63 * (ireg[i0] >> 1);
  i0 = iup[i0];
  i1 = iup[i1];
  return retval;
}

int ifrng(){
  int retval;
  ireg[i0] = ireg[i0] * ireg[i1];
  retval = ireg[i0] >> 33;
  i0 = iup[i0];
  i1 = iup[i1];
  return retval;
}

INT64 lfrng(){
  INT64 retval;
  ireg[i0] = ireg[i0] * ireg[i1];
  retval = ireg[i0] >> 1;
  i0 = iup[i0];
  i1 = iup[i1];
  return retval;
}

void frngv(int n, float *rvec){
  int i;
  for (i=0; i<n; i++){
    ireg[i0] = ireg[i0] * ireg[i1];
    rvec[i] = rnorm31 * (ireg[i0] >> 33);
    i0 = iup[i0];
    i1 = iup[i1];
  }
}

void dfrngv(int n, double *dvec){
  int i;
  for (i=0; i<n; i++){
    ireg[i0] = ireg[i0] * ireg[i1];
    dvec[i] = dnorm63 * (ireg[i0] >> 1);
    i0 = iup[i0];
    i1 = iup[i1];
  }
}

void ifrngv(int n, int *ivec){
  int i;
  for (i=0; i<n; i++){
    ireg[i0] = ireg[i0] * ireg[i1];
    ivec[i] = ireg[i0] >> 33;
    i0 = iup[i0];
    i1 = iup[i1];
  }
}

void lfrngv(int n, INT64 *lvec){
  int i;
  for (i=0; i<n; i++){
    ireg[i0] = ireg[i0] * ireg[i1];
    lvec[i] = ireg[i0] >> 1;
    i0 = iup[i0];
    i1 = iup[i1];
  }
}

void frnget(INT64 iseed[LEN]){
  int ipt, i;
  ipt = i0;
  for (i=0; i<LEN; i++){
    iseed[i]   = ireg[ipt] >> 1;
    ipt = iup[ipt];
  }
}

void frnset(INT64 iseed[LEN]){
  int i;
  for (i=0; i<LEN; i++){
    ireg[i] = ((unsigned INT64)iseed[i] << 1) + 1;
  }
  i0 = 0;
  i1 = LEN - LAG;
}

void frnini(int num){
  int mr=714025, ia=1366, ic=150889, mrhalf=(mr+1)/2;
  int i, ibit, irand;

  irand=labs(num) % mr;
  for (i=0; i<LEN; i++){
    ireg[i] = 0;
    for (ibit=0; ibit<=62; ibit++){
      irand = (ia*irand+ic) % mr;
      ireg[i] = 2*ireg[i] + irand/mrhalf;
    }
    ireg[i] = 2*ireg[i] + 1;   /* all odd! */
  }
  i0 = 0;
  i1 = LEN - LAG;

  /*
  for (i=0; i<LEN-1; i++) iup[i] = i+1;
  iup[LEN-1] = 0;
  */
}
