#include <string>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include <TRandom3.h>

using std::string;
using std::cout;
using std::cerr;
using std::endl;

typedef double *vektor;
typedef double **matrix;

void geom_pbc();
inline double p(double phi_re, double phi_im, double B_re, double B_im, double phi2, double lambda){ return exp(2*(B_re*phi_re+B_im*phi_im) - phi2 - lambda*(phi2-1)*(phi2-1)); }

#define N 10
#define ndim 2

// int **nn;
// int nvol;
// int lsize[ndim+1] = {0,N,N};
double delta = 0.5;

int main(){
  TRandom3 *ran = new TRandom3(0);
  double phi_re, phi_im, B_re, B_im, p_phi, lambda, phi2, kappa, r1, r2;
//   double h_re, h_im;
//   h_re = ran->Uniform();
//   h_im = ran->Uniform();
  phi_re = ran->Uniform();
  phi_im = ran->Uniform();
  kappa = ran->Uniform();
  lambda = ran->Uniform();

  B_re = ran->Uniform();
  B_im = ran->Uniform();

  phi2 = phi_re*phi_re + phi_im*phi_im;

  p_phi = p(phi_re, phi_im, B_re, B_im, phi2, lambda);

  double phi_neu_re, phi_neu_im, p_phi_neu;
  double p_accept;
  double p_grenz = 0.4;
  double wtrial, wmetro;

  for (int i=0; i<10; ++i){
    r1 = ran->Uniform(-1*delta, delta);
    r2 = ran->Uniform(-1*delta, delta);

    phi_neu_re = phi_re + r1;
    phi_neu_im = phi_im + r2;

    p_phi_neu = p(phi_neu_re, phi_neu_im, B_re, B_im, phi2, lambda);
    p_accept = (p_phi_neu >= p_phi) ? 1 : p_phi_neu/p_phi;
    cout << i << "   " << p_accept << endl;
    if(p_accept==1 || p_phi_neu/p_phi>0.5 ) p_phi = p_phi_neu;
  }
}


// void geom_pbc(){
//   /*
//   Angelegt und besetzt ist                            B Bunk 12/2005
//   Dimension     ndim                              rev     4/2011
//   Gittergroesse lsize[k], k=1..ndim
// 
//   Angelegt und berechnet wird
//   Volumen       nvol
//   NN-Indexfeld  nn[k][i], k=0..2*ndim, i=0..(nvol-1)
// 
//   nn[k][i] gibt den Index des Nachbarn von i in Richtung +k,
//   unter Beruecksichtigung periodischer Randbedingungen.
//   Fuer einen Schritt in Richtung -k setze man den Index auf (ndim+k).
//   nn[i][0] ist reserviert.
//   */
//   int   i, k;
//   int *ibase, *ix; 
// 
//   ibase = (int *) malloc((ndim+2) * sizeof(int));
//   ix = (int *) malloc((ndim+1) * sizeof(int));
// 
//   /* Basis fuer Punktindizes */
//   ibase[1] = 1;
//   for (k=1; k<=ndim; k++) ibase[k+1] = ibase[k]*lsize[k];
//   nvol = ibase[ndim+1];
// 
//   if (nn) free(nn[0]);
//   free(nn);
//   nn = (int **) malloc((2*ndim+1) * sizeof(int *));
//   nn[0] = (int *) malloc(nvol*(2*ndim+1) * sizeof(int));
//   for (k=1; k<=2*ndim; k++) nn[k] = nn[0] + nvol*k;
// 
//   for (k=1; k<=ndim; k++) ix[k] = 0;   /* Koord. des Anfangspunkts */
// 
//   for (i=0; i<nvol; i++){           /* Schleife ueber Punkte */
//     nn[0][i] = 0;
//     for (k=1; k<=ndim; k++){
//       nn[k][i] = i + ibase[k];      /* Nachbar x + e_k */
//       if (ix[k] == (lsize[k]-1)){
//         nn[k][i] -= ibase[k+1];
//         nn[0][i] = 1;               /* für Dirichlet-RB */
//       }
// 
//       nn[ndim+k][i] = i - ibase[k]; /* Nachbar x - e_k */
//       if (ix[k] == 0){
//         nn[ndim+k][i] += ibase[k+1];
//         nn[0][i] = 1;               /* für Dirichlet-RB */
//       }
// //       cout << "geomfkt: " << nn[k][i] << " " << nn[ndim+k][i] << endl;
//     }
// 
//     for (k=1; k<=ndim; k++){        /* Koord. des naechsten Punkts */
//       ix[k]++;
//       if (ix[k] < lsize[k]) break;
//       ix[k] = 0;
//     }
//   }
//   free(ibase); free(ix);
// }