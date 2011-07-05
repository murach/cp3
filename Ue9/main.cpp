#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#include <TRandom3.h>

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;

typedef double *vektor;
typedef double **matrix;
typedef int INT;
typedef double REAL;

void geom_vec();
inline double p(double phi_re, double phi_im, double B_re, double B_im, double phi2, double lambda){ return exp(2*(B_re*phi_re+B_im*phi_im) - phi2 - lambda*(phi2-1)*(phi2-1)); }

/* stat5 static data */
static INT  nvar = 0, nbmax = 0, *nbl, *lbl, *lnew;
static REAL **blksum, *sqsum;
static const REAL zero = 0., one = 1.;
/* prototypes for stat5 functions */
void  clear5(INT nvar, INT nbmax);
void  accum5(INT ivar, REAL value);
REAL  aver5(INT ivar);
REAL  var5(INT ivar);
REAL  sigma5(INT ivar);
REAL  cov5(INT ivar, INT jvar);
REAL  covar5(INT ivar, INT jvar);
REAL  tau5(INT ivar);
REAL  rsq5(INT ivar);
REAL  tauint5(INT ivar);
void  jackout5(INT ivar1, INT *nb, REAL bj[]);
void  jackeval5(INT nb, REAL fj[], REAL *aver, REAL *sigma);
void  jack5(REAL (*fct)(INT nvar, REAL a[]), REAL *aver, REAL *sigma);
void  save5(FILE *file);
void  savef5(FILE *file);
void  get5(FILE *file);
void  getf5(FILE *file);

#define n_mc_runs 10000
#define N 10
#define ndim 3

int **nn;
int lsize[4] = {0,N,N,N};
int nvol;
int nvcell, lvec, **nnstep, **nnflag;

int main(){
  geom_vec();
  TRandom3 *ran = new TRandom3(0);

  double phi_re[nvcell][lvec], phi_im[nvcell][lvec], B_re[lvec], B_im[lvec], p_phi[lvec], phi2[lvec], phi2_neu[lvec];
  double lambda = 2;
  double kappa = 0.3;
  double h = 0.05;

  double phi_neu_re[lvec], phi_neu_im[lvec], p_phi_neu[lvec];
  double p_accept;
  int akzeptanz;
  double dummy;
  double delta = 1.;
  double phi_nachbar_phi_re[lvec], phi_nachbar_phi_im[lvec];
  int n_hits = 10;
  int n_therm = 10000;
//   double ran_vek[n_therm+n_mc_runs][nvcell][3*n_hits];

  for (int ib=0; ib<nvcell; ++ib){
//     for (int j=0; j<3*n_hits; ++j){
//       for (int k=0; k<n_therm+n_mc_runs; ++k){
//         ran_vek[k][ib][j] = ran->Uniform();
//       }
//     }
    for (int l=0; l<lvec; ++l){
      phi_re[ib][l] = 0.;//ran->Uniform();
      phi_im[ib][l] = 0.;//ran->Uniform();
    }
  }

//   accumulated data:
//     1 : Re Phi
//     2 : Im Phi
//     3 : Re (|Phi|^2*Phi)
//     4 : Re (Phi_(x+mu)*Phi_stern)
//     5 : Im (Phi_(x+mu)*Phi_stern)
//     6 : |Phi|^2
//     7 : |Phi|^4
  clear5(7, 500);

  for (int k=0; k<n_therm+n_mc_runs; ++k){
    akzeptanz = 0;
    for (int ib=0; ib<nvcell; ++ib){

      for (int l=0; l<lvec; ++l){                       // B-Berechnung Anfang
        B_re[l] = h;
        B_im[l] = 0.;
//         akzeptanz = 0;
      }
      for (int m=1; m<=ndim*2; ++m){
        int jb = nnstep[m][ib];
        if (nnflag[m][ib]){
          for (int l=0; l<lvec; ++l){
            B_re[l] += kappa*phi_re[jb][nn[m][l]];
            B_im[l] += kappa*phi_im[jb][nn[m][l]];
          }
        }
        else{
          for (int l=0; l<lvec; ++l){
            B_re[l] += kappa*phi_re[jb][l];
            B_im[l] += kappa*phi_im[jb][l];
          }
        }
      }                                                 // B-Berechnung Ende

      for (int l=0; l<lvec; ++l){
        for (int i=0; i<n_hits; ++i){

          phi2[l] = phi_re[ib][l]*phi_re[ib][l] + phi_im[ib][l]*phi_im[ib][l];
          p_phi[l] = p(phi_re[ib][l], phi_im[ib][l], B_re[l], B_im[l], phi2[l], lambda);

//           ran_vek[k][ib][3*i + 0]
          phi_neu_re[l] = phi_re[ib][l] + ran->Uniform(-1*delta, delta);
          phi_neu_im[l] = phi_im[ib][l] + ran->Uniform(-1*delta, delta);

          phi2_neu[l] = phi_neu_re[l]*phi_neu_re[l] + phi_neu_im[l]*phi_neu_im[l];
          p_phi_neu[l] = p(phi_neu_re[l], phi_neu_im[l], B_re[l], B_im[l], phi2_neu[l], lambda);

          p_accept = (p_phi_neu[l] >= p_phi[l]) ? 1 : p_phi_neu[l]/p_phi[l];

//           ran_vek[k][ib][3*i+2]
          if(p_accept==1 || p_phi_neu[l]/p_phi[l]>ran->Uniform() ) {phi_re[ib][l] = phi_neu_re[l]; phi_im[ib][l] = phi_neu_im[l]; phi2[l] = phi2_neu[l]; ++akzeptanz;}
        }

        if (k>=n_therm){
          phi_nachbar_phi_re[l]=0.;
          phi_nachbar_phi_im[l]=0.;
          for (int i=1; i<=ndim; ++i){
            int jb = nnstep[i][ib];
            if (nnflag[i][ib]){
              phi_nachbar_phi_re[l] += phi_re[nnstep[i][ib]][l]*phi_re[nnstep[i][ib]][l] + phi_im[nnstep[i][ib]][l]*phi_im[nnstep[i][ib]][l];   //TODO
              phi_nachbar_phi_im[l] += phi_im[nnstep[i][ib]][l]*phi_re[nnstep[i][ib]][l] - phi_re[nnstep[i][ib]][l]*phi_im[nnstep[i][ib]][l];
            }
            else{
              phi_nachbar_phi_re[l] += phi_re[jb][l];;
              phi_nachbar_phi_re[l] += phi_im[jb][l];;
            }
          }

          accum5(1, phi_re[ib][l]);
          accum5(2, phi_im[ib][l]);
          accum5(3, phi2[l]*phi_re[ib][l]);
          accum5(4, phi_nachbar_phi_re[l]);
          accum5(5, phi_nachbar_phi_im[l]);
          accum5(6, phi2[l]);
          accum5(7, phi2[l]*phi2[l]);
        }
      }
    }
    
    if (k<n_therm){                      // TODO: wie geht das richtig?
      dummy = (double)akzeptanz/(nvol*n_hits);
      if (dummy < 0.35) delta *= 0.95;
      else if (dummy > 0.45) delta *= 1.05;
    }
  }

  cout << "Gl. 1: " << (1-2*ndim*kappa-2*lambda)*aver5(1) + 2*lambda*aver5(3) << " = " << h << endl;
  cout << "Gl. 2: " << -2*kappa*aver5(4) + (1-2*lambda)*aver5(6) + 2*lambda*aver5(7) << " = " << 1+h*aver5(1) << endl;
  cout << "Gl. 3: " << aver5(5) << " +- " << sigma5(5) << " = " << "0" << endl;

}


/* ***********************stat5***********************************
      statistical analysis using                      B Bunk 1993
            o data blocking                           rev  3/2006
            o autocorrelation analysis
            o jackknife

         nvar  : no. of variables
         nbmax : max. no. of blocks per variable
*/


/* allocate storage and clear counters */
void clear5(INT nvar1, INT nbmax1){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   ivar, ibl;

  if(nvar1 < 1){
    printf("error in clear5: nvar is invalid\n");
    exit(501);
  }
  if(nbmax1 < 2 || nbmax1 % 2 != 0){
    printf("error in clear5: nbmax is invalid\n");
    exit(502);
  }

  if(nvar1 != nvar || nbmax1 != nbmax){
    if(nvar > 0){
      free(nbl); free(lbl); free(lnew);
      free(blksum[0]); free(blksum);
      free(sqsum);
    }
    nbl = (INT *) malloc(nvar1 * sizeof(INT));
    lbl = (INT *) malloc(nvar1 * sizeof(INT));
    lnew = (INT *) malloc(nvar1 * sizeof(INT));
    sqsum = (REAL *) malloc(nvar1 * sizeof(REAL));
    blksum = (REAL **) malloc(nvar1 * sizeof(REAL *));
    if(!nbl || !lbl || !lnew || !blksum || !sqsum){
      printf("error in clear5: allocation failed\n");
      exit(503);
    }
    blksum[0] = (REAL *) malloc(nvar1*nbmax1 * sizeof(REAL));
    if(!blksum[0]){
      printf("error in clear5: allocation failed\n");
      exit(504);
    }
    for(ivar=1; ivar<nvar1; ivar++)
      blksum[ivar] = blksum[0] + ivar*nbmax1;
  }

  nvar = nvar1;
  nbmax = nbmax1;

  for(ivar=0; ivar<nvar; ivar++){
    nbl[ivar] = 0;
    lbl[ivar] = 1;
    lnew[ivar] = 0;
    for(ibl=0; ibl<nbmax; ibl++) blksum[ivar][ibl] = zero;
    sqsum[ivar] = zero;
  }
}

/* accumulate data */

void accum5(INT ivar1, REAL value){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   ivar, ibl, iblnew;

  if(ivar1 < 1 || ivar1 > nvar){
    printf("error in accum5: ivar out of range\n");
    exit(510);
  }
  ivar = ivar1 - 1;

  if(nbl[ivar] == nbmax){
    for (ibl=0; ibl<nbmax/2; ibl++)
      blksum[ivar][ibl]=blksum[ivar][2*ibl] + blksum[ivar][2*ibl+1];
    for (ibl=nbmax/2; ibl<nbmax; ibl++)
      blksum[ivar][ibl]=zero;
    nbl[ivar] = nbmax/2;
    lbl[ivar] = 2*lbl[ivar];
  }

  iblnew = nbl[ivar];
  blksum[ivar][iblnew] += value;
  sqsum[ivar] += value*value;
  lnew[ivar]++;

  if(lnew[ivar] == lbl[ivar]){
    nbl[ivar] = iblnew + 1;
    lnew[ivar] = 0;
  }
}

/* compute averages */

REAL aver5(INT ivar1){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   ivar, ibl, nmeas;
  REAL  res;

  if(ivar1 < 1 || ivar1 > nvar){
    printf("error in aver5: ivar out of range\n");
    exit(520);
  }
  ivar = ivar1 - 1;

  nmeas = nbl[ivar]*lbl[ivar] + lnew[ivar];
  if(nmeas == 0) return zero;
  res = zero;
  for(ibl=0; ibl<=nbl[ivar] && ibl<nbmax; ibl++)
    res += blksum[ivar][ibl];
  return res /= nmeas;
}

/* compute variances */

REAL var5(INT ivar1){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   ivar, nb, nmeas, ib, it;
  REAL  av, var, gam, gamvar;
  static REAL *d = NULL;

  if(ivar1 < 1 || ivar1 > nvar){
    printf("error in var5: ivar out of range\n");
    exit(530);
  }
  ivar = ivar1 - 1;

  nb = nbl[ivar];
  if(nb<2) return zero;
  nmeas = nb*lbl[ivar] + lnew[ivar];

  d = (REAL *) realloc(d, nbmax*sizeof(REAL));
  if(!d){
    printf("error in var5: allocation failed\n");
    exit(531);
  }

  av = zero;
  for(ib=0; ib<nb; ib++) av += blksum[ivar][ib];
  av /= nb;

  var = zero;
  for(ib=0; ib<nb; ib++){
    d[ib] = blksum[ivar][ib] - av;
    var += d[ib]*d[ib];
  }
  if(var <= zero) return zero;
  var /= nb;

  gamvar = var*var;
  for(it=1; it<nb; it++){
    gam = zero;
    for(ib=0; ib<nb-it; ib++) gam += d[ib]*d[ib+it];
    gam /= nb-it;
    if( gam <= sqrt(gamvar/(nb-it)) ) break;
    var += 2*gam;
    gamvar += 2*gam*gam;
  }

  return (var/nmeas)/lbl[ivar];
}

/* compute errors */

REAL sigma5(INT ivar1){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  if(ivar1 < 1 || ivar1 > nvar){
    printf("error in sigma5: ivar out of range\n");
    exit(540);
  }

  return sqrt( var5(ivar1) );
}

/* compute elements of the covariance matrix
         - sum off-diagonal gammas
         - can violate positivity
         better use covar5 !
*/
REAL cov5(INT ivar1, INT jvar1){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  REAL  cov, avi, avj, gamij , gamji, gamii, gamjj, gamvar, sgn;
  INT   ivar, jvar, nb, nmeas, ib, it;
  static REAL  *di = NULL, *dj = NULL;

  if(ivar1 < 1 || ivar1 > nvar){
    printf("error in cov5: ivar out of range\n");
    exit(550);
  }
  if(jvar1 < 1 || jvar1 > nvar){
    printf("error in cov5: jvar out of range\n");
    exit(551);
  }

  if(ivar1 == jvar1) return var5(ivar1);

  ivar = ivar1 - 1;
  jvar = jvar1 - 1;

  nb = nbl[ivar];
  if(nb < 2) return zero;
  nmeas = nb*lbl[ivar] + lnew[ivar];
  if( nmeas != nbl[jvar]*lbl[jvar]+lnew[jvar] ){
    printf("warning in cov5: mismatch of measurements\n");
    return zero;
  }

  di = (REAL *) realloc(di, nbmax*sizeof(REAL));
  dj = (REAL *) realloc(dj, nbmax*sizeof(REAL));
  if(!di || !dj){
    printf("error in cov5: allocation failed\n");
    exit(552);
  }

  avi = zero;
  avj = zero;
  for(ib=0; ib<nb; ib++){
    avi += blksum[ivar][ib];
    avj += blksum[jvar][ib];
  }
  avi /= nb;
  avj /= nb;

  cov = zero;
  gamii = zero;
  gamjj = zero;
  for(ib=0; ib<nb; ib++){
    di[ib] = blksum[ivar][ib] - avi;
    dj[ib] = blksum[jvar][ib] - avj;
    cov += di[ib]*dj[ib];
    gamii += di[ib]*di[ib];
    gamjj += dj[ib]*dj[ib];
  }
  cov /= nb;
  sgn = (cov > zero) ? one : -one;    /* sign of cov */
  gamvar = gamii*gamjj/(nb*nb) + cov*cov;

  for(it=1; it<nb; it++){
    gamij = zero;
    gamji = zero;
    gamii = zero;
    gamjj = zero;
    for(ib=0; ib<nb-it; ib++){
      gamij += di[ib+it]*dj[ib];
      gamji += dj[ib+it]*di[ib];
      gamii += di[ib]*di[ib+it];
      gamjj += dj[ib]*dj[ib+it];
    }
    gamij /= nb-it;
    gamji /= nb-it;
    gamii /= nb-it;
    gamjj /= nb-it;
    if( sgn*(gamij+gamji) <= sqrt(2*gamvar/(nb-it)) ) break;
    cov += gamij + gamji;
    gamvar += 2*(gamii*gamjj + gamij*gamji);
  }

  return (cov/nmeas)/lbl[ivar];
}

/* compute elements of the covariance matrix
         - rescale non-diagonal elements
         with autocorrelation factors
         from the diagonal
         - positive matrix
*/
REAL covar5(INT ivar1, INT jvar1){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  REAL  covar, avi, avj, di, dj, vari, varj;
  INT   ivar, jvar, nb, nmeas, ib;

  if(ivar1 < 1 || ivar1 > nvar){
    printf("error in covar5: ivar out of range\n");
    exit(555);
  }
  if(jvar1 < 1 || jvar1 > nvar){
    printf("error in covar5: jvar out of range\n");
    exit(556);
  }

  if(ivar1 == jvar1) return var5(ivar1);

  ivar = ivar1 - 1;
  jvar = jvar1 - 1;

  nb =nbl[ivar];
  if(nb < 2) return zero;
  nmeas = nb*lbl[ivar] + lnew[ivar];
  if( nmeas != nbl[jvar]*lbl[jvar]+lnew[jvar] ){
    printf("warning in cov5: mismatch of measurements\n");
    return zero;
  }

  avi = zero;
  avj = zero;
  for(ib=0; ib<nb; ib++){
    avi += blksum[ivar][ib];
    avj += blksum[jvar][ib];
  }
  avi /= nb;
  avj /= nb;

  covar = zero;
  vari = zero;
  varj = zero;
  for(ib=0; ib<nb; ib++){
    di = blksum[ivar][ib] - avi;
    dj = blksum[jvar][ib] - avj;
    covar += di*dj;
    vari += di*di;
    varj += dj*dj;
  }
  if(vari < zero || varj < zero) return zero;

  return covar * sqrt( var5(ivar1)*var5(jvar1)/(vari*varj) );
}

/*    estimate integrated auto-correlations:
         var(average) = var(single)/nmeas * rsq
         rsq = coth( 1/(2*tau) )
         = 2 * tauint
*/
REAL tau5(INT ivar1){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   ivar, nmeas;
  REAL  av, var, rsq;
  const REAL  oneeps = 1.000001e0;

  if(ivar1 < 1 || ivar1 > nvar){
    printf("error in tau5: ivar out of range\n");
    exit(560);
  }
  ivar = ivar1 - 1;

  nmeas = nbl[ivar]*lbl[ivar] + lnew[ivar];
  if(nmeas < 2) return zero;

  av = aver5(ivar1);
  var = sqsum[ivar]/nmeas - av*av;
  if(var <= zero) return zero;

  rsq = var5(ivar1) / var * nmeas;
  if(rsq <= oneeps) return zero;
  return one/log( (rsq+one)/(rsq-one) );
}
REAL rsq5(INT ivar1){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   ivar, nmeas;
  REAL  av, var;
  const REAL  oneeps = 1.000001e0;

  if(ivar1 < 1 || ivar1 > nvar){
    printf("error in rsq5: ivar out of range\n");
    exit(561);
  }
  ivar = ivar1 - 1;

  nmeas = nbl[ivar]*lbl[ivar] + lnew[ivar];
  if(nmeas < 2) return zero;

  av = aver5(ivar1);
  var = sqsum[ivar]/nmeas - av*av;
  if(var <= zero) return zero;

  return var5(ivar1) / var * nmeas;
}
REAL tauint5(INT ivar1){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   ivar, nmeas;
  REAL  av, var, rsq;
  const REAL  oneeps = 1.000001e0;

  if(ivar1 < 1 || ivar1 > nvar){
    printf("error in tauint5: ivar out of range\n");
    exit(562);
  }
  ivar = ivar1 - 1;

  nmeas = nbl[ivar]*lbl[ivar] + lnew[ivar];
  if(nmeas < 2) return zero;

  av = aver5(ivar1);
  var = sqsum[ivar]/nmeas - av*av;
  if(var <= zero) return zero;

  return .5 * var5(ivar1) / var * nmeas;
}


/* compute jackknife blocks

         ivar   : variable (input)
         nb     : no. of blocks (output)
         bj     : vector of jackknife blocks
         bj[ib], ib=0..nb-1 (output)
*/
void jackout5(INT ivar1, INT *nb1, REAL bj[]){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   ivar, nb, ib;
  REAL  bsum;

  if(ivar1 < 1 || ivar1 > nvar){
    printf("error in jackout5: ivar out of range\n");
    exit(570);
  }
  ivar = ivar1 - 1;

  nb = nbl[ivar];
  if(nb < 2){
    printf("error in jackout5: nb < 2\n");
    exit(571);
  }
  *nb1 = nb;

  bsum = zero;
  for(ib=0; ib<nb; ib++) bsum += blksum[ivar][ib];

  for(ib=0; ib<nb; ib++)
    bj[ib] = (bsum - blksum[ivar][ib])/(lbl[ivar]*(nb-1));
}

/* jackknife analysis for a vector of (function) values

        nb     : number of jackknife values, nb > 1 (input)
        fj     : vector of jackknife values,
        fj[ib], ib=0..nb-1 (input)
        aver   : function average (pointer, output)
        sigma  : error estimate (pointer, output)
*/
void jackeval5(INT nb, REAL fj[], REAL *aver, REAL *sigma){

  INT   ib, it;
  REAL  var, gam, gamvar;
  static REAL *d = NULL;

  *aver = zero;
  *sigma = zero;
  if(nb < 2) return;

  d = (REAL *) realloc(d, nb*sizeof(REAL));
  if(!d){
    printf("error in jackeval5: allocation failed\n");
    exit(575);
  }

  for(ib=0; ib<nb; ib++) *aver += fj[ib];
  *aver /= nb;

  gam = zero;
  for(ib=0; ib<nb; ib++){
    d[ib] = fj[ib] - *aver;
    gam += d[ib]*d[ib];
  }
  if(gam < zero) return;
  gam /= nb;

  var = gam;
  gamvar = gam*gam;
  for(it=1; it<nb; it++){
    gam = zero;
    for(ib=0; ib<nb-it; ib++) gam += d[ib]*d[ib+it];
    gam /= nb-it;
    if( gam <= sqrt(gamvar/(nb-it)) ) break;
    var += 2*gam;
    gamvar += 2*gam*gam;
  }

  *sigma = (nb-1)*sqrt(var/nb);
}

/* jackknife evaluation - custom call
  fct    : function of statistical variables (user-defined)
  REAL fct(INT nvar, REAL a[])
  with a[i], i=1..nvar : variables, as in stat5
  aver   : function average (pointer, output)
  sigma  : error estimate (pointer, output)

  calls jackeval5

  note: is is assumed that all (relevant) variables were accumulated
  in sync, with ivar=1 as prototype.
*/
void jack5(REAL (*fct)(INT nvar, REAL a[]), REAL *aver, REAL *sigma){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   nb, ivar, ib;
  static REAL *bsum = NULL, *bj = NULL, *fj = NULL, *bj1;

  *aver = zero;
  *sigma = zero;
  nb = nbl[1];
  if(nb < 2) return;

  bsum = (REAL *) realloc(bsum, nvar*sizeof(REAL));
  bj = (REAL *) realloc(bj, nvar*sizeof(REAL));
  fj = (REAL *) realloc(fj, nbmax*sizeof(REAL));
  if(!bsum || !bj || !fj){
    printf("error in jack5: allocation failed\n");
    exit(576);
  }

  bj1 = bj - 1;        /* map bj[0..nvar-1] <-> bj1[1..nvar] */

  for(ivar=0; ivar<nvar; ivar++){
    bsum[ivar] = zero;
    for(ib=0; ib<nb; ib++) bsum[ivar] += blksum[ivar][ib];
  }

  for(ib=0; ib<nb; ib++){
    for(ivar=0; ivar<nvar; ivar++)
      bj[ivar] = (bsum[ivar] - blksum[ivar][ib])/(lbl[1]*(nb-1));
    fj[ib] = fct(nvar, bj1);
  }

  jackeval5(nb, fj, aver, sigma);
}

/* save and restore counters (unformatted/formatted) */

void save5(FILE *file){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   ivar;

  fwrite(&nvar, sizeof(INT), 1, file);
  fwrite(&nbmax, sizeof(INT), 1, file);
  for(ivar=0; ivar<nvar; ivar++){
    fwrite(&nbl[ivar], sizeof(INT), 1, file);
    fwrite(&lbl[ivar], sizeof(INT), 1, file);
    fwrite(&lnew[ivar], sizeof(INT), 1, file);
  }
  fwrite(blksum[0], sizeof(REAL), nvar*nbmax, file);
  fwrite(sqsum, sizeof(REAL), nvar, file);
}
void savef5(FILE *file){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   ivar, ib;

//   fprintf(file, "%i %i\n", nvar, nbmax);
//   for(ivar=0; ivar<nvar; ivar++)
//     fprintf(file, "%i %i %i\n", nbl[ivar], lbl[ivar], lnew[ivar]);
//   for(ivar=0; ivar<nvar; ivar++)
//     for(ib=0; ib<nbmax; ib++)
//       fprintf(file, "%20.16e\n", blksum[ivar][ib]);
//   for(ivar=0; ivar<nvar; ivar++)
//     fprintf(file, "%20.16e\n", sqsum[ivar]);

  for(ivar=0; ivar<nvar; ivar++)
    for(ib=0; ib<nbmax; ib++)
      fprintf(file, "%20.16e\n", blksum[ivar][ib]);
  for(ivar=0; ivar<nvar; ivar++)
    fprintf(file, "%20.16e\n", sqsum[ivar]);
}

void get5(FILE *file){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   nvar1, nbmax1, ivar;

  fread(&nvar1, sizeof(INT), 1, file);
  fread(&nbmax1, sizeof(INT), 1, file);

  if(nvar1 < 1){
    printf("error in get5: nvar is invalid\n");
    exit(580);
  }
  if(nbmax1 < 2 || nbmax1 % 2 != 0){
    printf("error in get5: nbmax is invalid\n");
    exit(581);
  }

  if(nvar1 != nvar || nbmax1 != nbmax){
    if(nvar > 0){
      free(nbl); free(lbl); free(lnew);
      free(blksum[0]); free(blksum);
      free(sqsum);
    }
    nbl = (INT *) malloc(nvar1 * sizeof(INT));
    lbl = (INT *) malloc(nvar1 * sizeof(INT));
    lnew = (INT *) malloc(nvar1 * sizeof(INT));
    sqsum = (REAL *) malloc(nvar1 * sizeof(REAL));
    blksum = (REAL **) malloc(nvar1 * sizeof(REAL *));
    if(!nbl || !lbl || !lnew || !blksum || !sqsum){
      printf("error in get5: allocation failed\n");
      exit(582);
    }
    blksum[0] = (REAL *) malloc(nvar1*nbmax1 * sizeof(REAL));
    if(!blksum[0]){
      printf("error in get5: allocation failed\n");
      exit(583);
    }
    for(ivar=1; ivar<nvar1; ivar++)
      blksum[ivar] = blksum[0] + ivar*nbmax1;
  }

  nvar = nvar1;
  nbmax = nbmax1;

  for(ivar=0; ivar<nvar; ivar++){
    fread(&nbl[ivar], sizeof(INT), 1, file);
    fread(&lbl[ivar], sizeof(INT), 1, file);
    fread(&lnew[ivar], sizeof(INT), 1, file);
  }
  fread(blksum[0], sizeof(REAL), nvar*nbmax, file);
  fread(sqsum, sizeof(REAL), nvar, file);
}
void getf5(FILE *file){

  extern INT  nvar, nbmax, *nbl, *lbl, *lnew;
  extern REAL **blksum, *sqsum;

  INT   nvar1, nbmax1, ivar, ib;

  fscanf(file, "%i %i\n", &nvar1, &nbmax1);

  if(nvar1 < 1){
    printf("error in getf5: nvar is invalid\n");
    exit(585);
  }
  if(nbmax1 < 2 || nbmax1 % 2 != 0){
    printf("error in getf5: nbmax is invalid\n");
    exit(586);
  }

  if(nvar1 != nvar || nbmax1 != nbmax){
    if(nvar > 0){
      free(nbl); free(lbl); free(lnew);
      free(blksum[0]); free(blksum);
      free(sqsum);
    }
    nbl = (INT *) malloc(nvar1 * sizeof(INT));
    lbl = (INT *) malloc(nvar1 * sizeof(INT));
    lnew = (INT *) malloc(nvar1 * sizeof(INT));
    sqsum = (REAL *) malloc(nvar1 * sizeof(REAL));
    blksum = (REAL **) malloc(nvar1 * sizeof(REAL *));
    if(!nbl || !lbl || !lnew || !blksum || !sqsum){
      printf("error in getf5: allocation failed\n");
      exit(587);
    }
    blksum[0] = (REAL *) malloc(nvar1*nbmax1 * sizeof(REAL));
    if(!blksum[0]){
      printf("error in getf5: allocation failed\n");
      exit(588);
    }
    for(ivar=1; ivar<nvar1; ivar++)
      blksum[ivar] = blksum[0] + ivar*nbmax1;
  }

  nvar = nvar1;
  nbmax = nbmax1;

  for(ivar=0; ivar<nvar; ivar++)
    fscanf(file, "%i %i %i\n", &nbl[ivar], &lbl[ivar], &lnew[ivar]);
  for(ivar=0; ivar<nvar; ivar++)
    for(ib=0; ib<nbmax; ib++)
      fscanf(file, "%le\n", &blksum[ivar][ib]);
  for(ivar=0; ivar<nvar; ivar++)
    fscanf(file, "%le\n", &sqsum[ivar]);
}

void geom_vec() {
  /*
  Angelegt und besetzt ist                            B Bunk 6/2010
      Dimension     ndim                              rev    6/2011
      Gittergroesse lsize[k], k=1..ndim

  Angelegt und berechnet wird
      Volumen                    nvol

      Volumen der Basiszelle     nvcell
      NN-Schritt fuer Basispunkt nnstep[k][ib], ib=0..(nvcell-1), k=0..2*ndim
      Permutationsflag  dazu     nnflag[k][ib], ib=0..(nvcell-1), k=0..2*ndim

      Vektorlaenge               lvec
      NN-Indexfeld im Vektor     nn[k][iv], iv=0..(lvec-1), k=0..2*ndim

  nnstep[k][ib] gibt den Nachbarn des Basispunkts ib bei einem Schritt in
      Richtung +k.
      Fuer einen Schritt in Richtung -k setze man den Index auf (ndim+k).
      nnstep[0][ib] ist reserviert.
  nnflag[k][ib] gibt die zugehoerige Permutationsnummer:
      nnflag[k][ib] = k, (k=1..2*ndim), wenn der Schritt in der Basiszelle
                           ueber den Rand springt;
      nnflag[k][ib] = 0 sonst.
      nnflag[0][ib] = 0,1 gibt die Farbe von ib (Schachbrettmuster).

  nn[k][iv] gibt, ausgehend vom Vektorindex iv, den Index im Nachbarvektor
      bei einer Verschiebung (des gesamten Vektors) in Richtung +k.
      Fuer eine Verschiebung in Richtung -k setze man den Index auf (ndim+k).
      nn[0][iv] = iv (keine Permutation).
  */
  int i, k, lsv;
  int *ibase, *ix, *lsblk, *lsvec;

  /* Hilfsfelder anlegen (ANSI C) */

  ibase = (int *) malloc((ndim+2) * sizeof(int));
  ix = (int *) malloc((ndim+1) * sizeof(int));
  lsblk = (int *) malloc((ndim+1) * sizeof(int));
  lsvec = (int *) malloc((ndim+1) * sizeof(int));

  /* Gittergroessen */

  nvol = 1;
  nvcell = 1;
  lvec = 1;
  for (k=1; k<=ndim; k++) {
    if (lsize[k] % 2) {
      printf("error in geom_vec: lsize[%i] is odd\n", lsize[k]);
      exit(10);
    }
    for (lsv=10; lsv>0; lsv--) {       /* suche Zerlegung */
      if (lsize[k] % (2*lsv) == 0) {   /* lsize[k] = lsblk[k] * lsvec[k] */
        lsvec[k] = lsv;
        lsblk[k] = lsize[k]/lsv;
        break;
      }
    }
    nvol *= lsize[k];
    nvcell *= lsblk[k];
    lvec *= lsvec[k];
  }

  /* Indexfelder neu anlegen (ANSI C) */

  if (nn) {free(nn[0]); free(nnstep[0]); free(nnflag[0]);}
  free(nn); free(nnstep); free(nnflag);

  nnstep = (int **) malloc((2*ndim+1) * sizeof(int *));
  nnstep[0] = (int *) malloc((2*ndim+1)*nvcell * sizeof(int));
  for (k=1; k<=2*ndim; k++) nnstep[k] = nnstep[0] + nvcell*k;

  nnflag = (int **) malloc((2*ndim+1) * sizeof(int *));
  nnflag[0] = (int *) malloc((2*ndim+1)*nvcell * sizeof(int));
  for (k=1; k<=2*ndim; k++) nnflag[k] = nnflag[0] + nvcell*k;

  nn = (int **) malloc((2*ndim+1) * sizeof(int *));
  nn[0] = (int *) malloc((2*ndim+1)*lvec * sizeof(int));
  for (k=1; k<=2*ndim; k++) nn[k] = nn[0] + lvec*k;

  /* Geometriefelder fuer Basiszelle */

  ibase[1] = 1;
  for (k=1; k<=ndim; k++) {
    ibase[k+1] = ibase[k]*lsblk[k];
    ix[k] = 0;
  }

  for (i=0; i<nvcell; i++) {              /* Schleife ueber Basiszelle */
    nnflag[0][i] = 0;
    for (k=1; k<=ndim; k++) nnflag[0][i] += ix[k];
    nnflag[0][i] = nnflag[0][i] % 2;      /* Kolorierung */

    for (k=1; k<=ndim; k++) {
      nnstep[k][i] = i + ibase[k];        /* Blocknachbar x + e_k */
      nnflag[k][i] = 0;
      if (ix[k] == (lsblk[k]-1)) {
        nnstep[k][i] -= ibase[k+1];
        nnflag[k][i] = k;
      }

      nnstep[ndim+k][i] = i - ibase[k];   /* Blocknachbar x - e_k */
      nnflag[ndim+k][i] = 0;
      if (ix[k] == 0) {
        nnstep[ndim+k][i] += ibase[k+1];
        nnflag[ndim+k][i] = ndim+k;
      }
    }

    for (k=1; k<=ndim; k++) {             /* naechster Punkt */
      ix[k]++;
      if (ix[k] < lsblk[k]) break;
      ix[k] = 0;
    }
  }                              /* Ende der Schleife ueber Basiszelle */

  /* Indexfeld fuer Vektor */

  ibase[1] = 1;
  for (k=1; k<=ndim; k++) {
    ibase[k+1] = ibase[k]*lsvec[k];
    ix[k] = 0;
  }

  for (i=0; i<lvec; i++) {       /* Schleife ueber Vektorindizes */
    nn[0][i] = i;                         /* keine Permutation */

    for (k=1; k<=ndim; k++) {
      nn[k][i] = i + ibase[k];            /* Nachbar (+k) */
      if (ix[k] == (lsvec[k]-1)) nn[k][i] -= ibase[k+1];

      nn[ndim+k][i] = i - ibase[k];       /* Nachbar (-k) */
      if (ix[k] == 0) nn[ndim+k][i] += ibase[k+1];
    }

    for (k=1; k<=ndim; k++) {             /* naechster Punkt */
      ix[k]++;
      if (ix[k] < lsvec[k]) break;
      ix[k] = 0;
    }
  }                              /* Ende der Schleife ueber Vektorindizes */

  /* Hilfsfelder freigeben */

  free(ibase); free(ix); free(lsblk); free(lsvec);
}