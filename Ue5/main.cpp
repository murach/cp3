#include <iostream>
#include <stdlib.h>
#include <math.h>

#include <TRandom3.h>

using std::string;
using std::cout;
using std::cerr;
using std::endl;

typedef double *vektor;
typedef double **matrix;

void geom_pbc();
inline double p(double phi_re, double phi_im, double B_re, double B_im, double phi2, double lambda){ return exp(2*(B_re*phi_re+B_im*phi_im) - phi2 - lambda*(phi2-1)*(phi2-1)); }

#define n_mc_runs 1e5

int main(){
  TRandom3 *ran = new TRandom3(0);
  double phi_re, phi_im, B_re, B_im, p_phi, phi2, phi2_neu, kappa, r1, r2;
  phi_re = ran->Uniform();
  phi_im = ran->Uniform();
  kappa = ran->Uniform();
  double lambda[2]={0, ran->Uniform()};

  B_re = ran->Uniform();        // fuer uns soll B erstmal konstant sein
  B_im = ran->Uniform();
  double B2 = B_re*B_re + B_im*B_im;

  double phi_neu_re, phi_neu_im, p_phi_neu;
  double p_accept;
  double p_grenz = 0.5;
  double delta = 2.;
  double phi_mean_re=0., phi_mean_im=0., phi2_mean=0.;
  double dummy_re=0., dummy_im=0.;
  int counter;

  for (int j=0; j<2; ++j){
    counter = 0;
    for (int k=0; k<n_mc_runs; ++k){
      for (int i=0; i<10; ++i){
        r1 = ran->Uniform(-1*delta, delta);
        r2 = ran->Uniform(-1*delta, delta);

        phi2 = phi_re*phi_re + phi_im*phi_im;
        p_phi = p(phi_re, phi_im, B_re, B_im, phi2, lambda[j]);

        phi_neu_re = phi_re + r1;
        phi_neu_im = phi_im + r2;

        phi2_neu = phi_neu_re*phi_neu_re + phi_neu_im*phi_neu_im;
        p_phi_neu = p(phi_neu_re, phi_neu_im, B_re, B_im, phi2_neu, lambda[j]);

        p_accept = (p_phi_neu >= p_phi) ? 1 : p_phi_neu/p_phi;

        if(p_accept==1 || p_phi_neu/p_phi>ran->Uniform() ) {phi_re = phi_neu_re; phi_im = phi_neu_im; phi2 = phi2_neu; counter++;}
      }
      phi_mean_re += phi_re;
      phi_mean_im += phi_im;
      phi2_mean += phi2;
      dummy_re += phi_re*(phi2-1);
      dummy_im += phi_im*(phi2-1);
    }
    phi_mean_re /= n_mc_runs;
    phi_mean_im /= n_mc_runs;
    phi2_mean /= n_mc_runs;
    dummy_re /= n_mc_runs;
    dummy_im /= n_mc_runs;

    cout << endl << "Akzeptanz: " << counter*1./(10*n_mc_runs) << endl;

    if (j == 0){
      cout << endl << "######### Checks for lambda = 0: #########" << endl;
      cout << "phi_re: " << phi_mean_re << "  , B_re: " << B_re << endl;
      cout << "phi_im: " << phi_mean_im << "  , B_im: " << B_im << endl;
      cout << "phi2  : " << phi2_mean << "   , 1+|B|2: " << 1 + B2 << endl << endl;
    } else {
      cout << endl << "######### Checks for lambda > 0: #########" << endl;
      cout << "phi_re: " << phi_mean_re << ", Re(B-2*lambda*(...)): " << B_re - 2*lambda[j]*dummy_re << endl;
      cout << "phi_im: " << phi_mean_im << ", Im(B-2*lambda*(...)): " << B_im - 2*lambda[j]*dummy_im << endl << endl;
    }
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