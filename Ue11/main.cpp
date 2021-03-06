#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>

#include "frng.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;

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

const int n_mc_runs = 100;
int N;
int ndim;

int **nn;
int lsize[4];
int nvol;
int nvcell, lvec, **nnstep, **nnflag;

int main(int argc, char *argv[]){
  ndim = atoi(argv[1]);
  const double lambda = atof(argv[2]);
  const double kappa = atof(argv[3]);
  const double h = atof(argv[4]);
  const int Lmax = atoi(argv[5]);

  FILE *file = fopen("output.dat", "w");
  fprintf(file, "# ndim M2 M4 U4\n");

  for (N=4; N<=Lmax; N+=2){
    lsize[0] = 0; lsize[1] = N; lsize[2] = N; lsize[3] = N;
    geom_vec();
    #pragma omp parallel
    frnini(omp_get_thread_num());

    double phi_re[nvcell][lvec], phi_im[nvcell][lvec], B_re[lvec], B_im[lvec], p_phi[lvec], phi2[lvec], phi2_neu[lvec];
    double M_re, M_im, phi2_phi_re, phi_nachbar_mal_phi_re, phi_nachbar_mal_phi_im, phi2_skalar, phi4, M2;

    double phi_neu_re[lvec], phi_neu_im[lvec], p_phi_neu[lvec];
    double p_accept[lvec];
    int akzeptanz;
    double dummy;
    double delta = 1.;
    double phi_nachbar_phi_re[lvec], phi_nachbar_phi_im[lvec];
    const int n_hits = 10;
    const int n_therm = 10000;
    double ran_vek1[lvec], ran_vek2[lvec], ran_vek3[lvec];
    int map[2][nvcell/2], icolor, ib, j;

    int l1;
    #pragma omp parallel for private(l1)
    for (ib=0; ib<nvcell; ++ib){
      for (l1=0; l1<lvec; ++l1){
        phi_re[ib][l1] = 0.;
        phi_im[ib][l1] = 0.;
      }
      map[nnflag[0][ib]][ib/2] = ib;
    }

  //   accumulated data:
  //     1 : Re Phi
  //     2 : Im Phi
  //     3 : Re (|Phi|^2*Phi)
  //     4 : Re (Phi_(x+mu)*Phi_stern)
  //     5 : Im (Phi_(x+mu)*Phi_stern)
  //     6 : |Phi|^2
  //     7 : |Phi|^4
  //     8 : |M|^2
  //     9 : |M|^4
    clear5(9, 500);

    for (int k=0; k<n_therm+n_mc_runs; ++k){
      akzeptanz=0; M_re=0; M_im=0; phi2_phi_re=0; phi_nachbar_mal_phi_re=0; phi_nachbar_mal_phi_im=0; phi2_skalar=0; phi4=0;

      for (icolor=0; icolor<2; ++icolor){
        int l, m, i, jb;
        #pragma omp parallel for private(j,i,l,m,jb,ib,B_re,B_im,phi2,p_phi,ran_vek1,ran_vek2,ran_vek3,phi_neu_re,\
        phi_neu_im,phi2_neu, p_phi_neu, p_accept, phi_nachbar_phi_re, phi_nachbar_phi_im)\
        reduction(+:akzeptanz, M_re, M_im, phi2_phi_re, phi_nachbar_mal_phi_re, phi_nachbar_mal_phi_im, phi2_skalar, phi4)
        for (j=0; j<nvcell/2; ++j){
          ib = map[icolor][j];

          for (l=0; l<lvec; ++l){                       // B-Berechnung Anfang
            B_re[l] = h;
            B_im[l] = 0.;
          }

          for (m=1; m<=ndim*2; ++m){
            jb = nnstep[m][ib];
            if (nnflag[m][ib]){
              for (l=0; l<lvec; ++l){
                B_re[l] += kappa*phi_re[jb][nn[m][l]];
                B_im[l] += kappa*phi_im[jb][nn[m][l]];
              }
            }
            else{
              for (l=0; l<lvec; ++l){
                B_re[l] += kappa*phi_re[jb][l];
                B_im[l] += kappa*phi_im[jb][l];
              }
            }
          }                                                 // B-Berechnung Ende

          for (l=0; l<lvec; ++l){
            phi2[l] = phi_re[ib][l]*phi_re[ib][l] + phi_im[ib][l]*phi_im[ib][l];
            p_phi[l] = p(phi_re[ib][l], phi_im[ib][l], B_re[l], B_im[l], phi2[l], lambda);
          }

          for (i=0; i<n_hits; ++i){
            dfrngv(lvec, ran_vek1);
            dfrngv(lvec, ran_vek2);
            dfrngv(lvec, ran_vek3);
            for (l=0; l<lvec; ++l){
              ran_vek1[l] = 2*delta*ran_vek1[l] - delta;
              ran_vek2[l] = 2*delta*ran_vek2[l] - delta;
            }

            for (l=0; l<lvec; ++l){
              phi_neu_re[l] = phi_re[ib][l] + ran_vek1[l];
              phi_neu_im[l] = phi_im[ib][l] + ran_vek2[l];

              phi2_neu[l] = phi_neu_re[l]*phi_neu_re[l] + phi_neu_im[l]*phi_neu_im[l];
              p_phi_neu[l] = p(phi_neu_re[l], phi_neu_im[l], B_re[l], B_im[l], phi2_neu[l], lambda);

              p_accept[l] = (p_phi_neu[l] >= p_phi[l]) ? 1 : p_phi_neu[l]/p_phi[l];

              if (p_accept[l]==1 || p_phi_neu[l]/p_phi[l]>ran_vek3[l]) {p_phi[l]=p_phi_neu[l]; phi_re[ib][l]=phi_neu_re[l]; phi_im[ib][l]=phi_neu_im[l]; phi2[l]=phi2_neu[l]; ++akzeptanz;}
              
              if (k>=n_therm && i+1==n_hits){
                phi_nachbar_phi_re[l]=0.;
                phi_nachbar_phi_im[l]=0.;
                for (m=1; m<=ndim; ++m){
                  jb = nnstep[m][ib];
                  if (nnflag[m][ib]){
                    phi_nachbar_phi_re[l] += phi_re[nnstep[m][ib]][l]*phi_re[nnstep[m][ib]][l] + phi_im[nnstep[m][ib]][l]*phi_im[nnstep[m][ib]][l];
                    phi_nachbar_phi_im[l] += phi_im[nnstep[m][ib]][l]*phi_re[nnstep[m][ib]][l] - phi_re[nnstep[m][ib]][l]*phi_im[nnstep[m][ib]][l];
                  }
                  else{
                    phi_nachbar_phi_re[l] += phi_re[jb][l];;
                    phi_nachbar_phi_im[l] += phi_im[jb][l];;
                  }
                }
                M_re += phi_re[ib][l];
                M_im += phi_im[ib][l];
                phi2_phi_re += phi2[l]*phi_re[ib][l];
                phi_nachbar_mal_phi_re += phi_nachbar_phi_re[l];
                phi_nachbar_mal_phi_im += phi_nachbar_phi_im[l];
                phi2_skalar += phi2[l];
                phi4 += phi2[l]*phi2[l];

              }             // end if
            }               // end for l<lvec
          }                 // end hits loop (i)
        }                   // end cell loop (j)
      }                     // end color loop (icolor)
      
      if (k>=n_therm){
        M_re /= nvol;
        M_im /= nvol;
        M2 = M_re*M_re + M_im*M_im;
        accum5(1, M_re);
        accum5(2, M_im);
        accum5(3, phi2_phi_re/nvol);
        accum5(4, phi_nachbar_mal_phi_re/nvol);
        accum5(5, phi_nachbar_mal_phi_im/nvol);
        accum5(6, phi2_skalar/nvol);
        accum5(7, phi4/nvol);
        accum5(8, M2);
        accum5(9, M2*M2);
      }
      
      if (k<n_therm){
        dummy = (double)akzeptanz/(nvol*n_hits);
        if (dummy < 0.35) delta *= 0.95;
        else if (dummy > 0.45) delta *= 1.05;
      }
    }                     // end mc-run loop (k)

    fprintf(file, "%d %f %f %f %f %f\n", N, aver5(8), sigma5(8), aver5(9), sigma5(9), 2-aver5(9)/(aver5(8)*aver5(8)));

  }                     // end L loop (L)
  fclose(file);
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