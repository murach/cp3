#include <string>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include <TRandom3.h>

// #include "global.h"

using std::string;
using std::cout;
using std::endl;

typedef double *vektor;
typedef double **matrix;

#define N 9
#define ndim 2

matrix malloc_matrix();
inline vektor malloc_vektor(){ return new double[N*N]; }
void vec_copy(vektor a, vektor b);
void vec_addition(vektor a, vektor b, vektor c, int sign = 1);
vektor mult_matrix_vektor(matrix A, vektor x);
void mult_skalar_vektor(double alpha, vektor x, vektor c);
double skalarprodukt(vektor a, vektor b);
void print_vektor(vektor x);
void print_matrix(matrix A);
void laplace(vektor x, double m2, vektor c);

void cg(vektor x, vektor b, void (*fkt)(vektor x, double m2, vektor c), int max_it, double relerr = 1e-10, bool flag = 1);
void geom_pbc();

int **nn;
int lsize[ndim+1] = {0,N,N};
int nvol;

int main(int argc, char *argv[]) {

    TRandom3 *ran = new TRandom3(0);

    double eta[N*N];
    double phi[N*N];
    geom_pbc();

//     matrix A = malloc_matrix();
//     double A[N][N];

    for (int i=0; i<nvol; ++i){
      eta[i] = ran->Uniform();
    }

    cg(phi, eta, laplace, 1000, 1e-10, 1);

    return 0;
}

matrix malloc_matrix()
{
  matrix feld;
  feld = (matrix)malloc(N*sizeof(vektor));
  feld[0] = (vektor)malloc(N*N*sizeof(double));
  if(feld[0] == NULL) {
    fprintf(stderr, "Not enough memory for allocating matrix\n");
    exit(1);
  }
  for(int x = 1; x < N; x++)
    feld[x] = feld[0] + N * x;
  return feld;
}

void vec_copy(vektor a, vektor b){
  for (int i=0; i<nvol; ++i){
    b[i] = a[i];
  }
  return;
}

void vec_addition(vektor a, vektor b, vektor c, int sign){      // sign: defaultparameter +1
  sign = (sign>0) ? 1 : -1;
  for (int i=0; i<nvol; ++i){
    c[i] = a[i] + sign*b[i];
  }
}

vektor mult_matrix_vektor(matrix A, vektor x){
  vektor c = malloc_vektor();
//   #pragma omp parallel for
  for (int i=0; i<N; ++i){
    c[i] = 0.;
    for (int j=0; j<N; ++j){
      c[i] += A[i][j]*x[j];
    }
  }
  return c;
}

void laplace(vektor x, double m2, vektor c){
  #pragma omp parallel for
  for (int i=0; i<nvol; ++i){
    c[i] = (2*ndim + m2)*x[i];
    for (int j=1; j<=ndim; ++j){
      c[i] -= x[nn[j][i]] + x[nn[ndim+j][i]];
    }
  }
  return;
}

void mult_skalar_vektor(double alpha, vektor x, vektor c){
  #pragma omp parallel for
  for (int i=0; i<nvol; ++i){
    c[i] = alpha*x[i];
  }
  return;
}

double skalarprodukt(vektor a, vektor b){
  double erg = 0;
  #pragma omp parallel for reduction(+:erg)
  for (int i=0; i<nvol; ++i){
    erg += a[i]*b[i];
  }
  return erg;
}

void cg(vektor x, vektor b, void (*fkt)(vektor x, double m2, vektor c), int max_it, double relerr, bool flag){
  double r[N*N];
  double tol = relerr*relerr*skalarprodukt(b, b);
  double m = 0.1;
  double m2 = m*m;
  double dummyvec[N*N];
  double dummyvec2[N*N];

  if (flag){
    vec_copy(b, r);
    for (int i=0; i<nvol; ++i){
      x[i] = 0.;
    }
  }
//   else{
//     dummyvec = mult_matrix_vektor(A,x);
//     vec_addition(b, dummyvec, r, -1);
//   }

  if (skalarprodukt(r, r) < tol) exit(0);

  double p[N*N];
  vec_copy(r, p);
  double s[N*N];
  double alpha, r2_alt, r2_neu, beta;
  int counter = 0;

  for (int i=0; i<max_it; ++i){
    fkt(p, m2, s);
    r2_alt = skalarprodukt(r, r);                       // r²_k;
    alpha = r2_alt / skalarprodukt(p, s);
//     mult_skalar_vektor(alpha, p, dummyvec);
//     vec_addition(x, dummyvec, dummyvec2);		        // neues x_(k+1)
    mult_skalar_vektor(alpha, s, dummyvec);                     //checked; dummyvec = alpha * s
    vec_addition(r, dummyvec, dummyvec2, -1);	                // neues r_(k+1), checked; dummyvec2 = r - dummyvec
    vec_copy(dummyvec2, r);                                     // r = dummyvec2
    r2_neu = skalarprodukt(r, r);			// r²_(k+1)
    if (r2_neu < tol) break;
    beta = r2_neu / r2_alt;
    mult_skalar_vektor(beta, p, dummyvec);              // checked; dummyvec = beta * p
    vec_addition(r, dummyvec, p);               	// neues p_(k+1), checked; p = r + dummyvec
    ++counter;
  }

  cout << "Anzahl von Schritten: " << counter << endl;
  return;
}

void print_vektor(vektor x){
  cout << "[";
  for (int i=0; i<nvol; ++i){
    cout << " " << x[i];
  }
  cout << " ]" << endl;
  return;
}

void print_matrix(matrix A){
  for (int i=0; i<N; ++i){
    for (int j=0; j<N; ++j){
      cout << " " << A[i][j];
    }
    cout << endl;
  }
  return;
}

void geom_pbc(){
  /*
  Angelegt und besetzt ist                            B Bunk 12/2005
  Dimension     ndim                              rev     4/2011
  Gittergroesse lsize[k], k=1..ndim

  Angelegt und berechnet wird
  Volumen       nvol
  NN-Indexfeld  nn[k][i], k=0..2*ndim, i=0..(nvol-1)

  nn[k][i] gibt den Index des Nachbarn von i in Richtung +k,
  unter Beruecksichtigung periodischer Randbedingungen.
  Fuer einen Schritt in Richtung -k setze man den Index auf (ndim+k).
  nn[i][0] ist reserviert.
  */
  int   i, k;
  int *ibase, *ix; 

  ibase = (int *) malloc((ndim+2) * sizeof(int));
  ix = (int *) malloc((ndim+1) * sizeof(int));

  /* Basis fuer Punktindizes */
  ibase[1] = 1;
  for (k=1; k<=ndim; k++) ibase[k+1] = ibase[k]*lsize[k];
  nvol = ibase[ndim+1];

  if (nn) free(nn[0]);
  free(nn);
  nn = (int **) malloc((2*ndim+1) * sizeof(int *));
  nn[0] = (int *) malloc(nvol*(2*ndim+1) * sizeof(int));
  for (k=1; k<=2*ndim; k++) nn[k] = nn[0] + nvol*k;

  for (k=1; k<=ndim; k++) ix[k] = 0;   /* Koord. des Anfangspunkts */

  for (i=0; i<nvol; i++){           /* Schleife ueber Punkte */
    nn[0][i] = 0;
    for (k=1; k<=ndim; k++){
      nn[k][i] = i + ibase[k];      /* Nachbar x + e_k */
      if (ix[k] == (lsize[k]-1)){
	nn[k][i] -= ibase[k+1];
	nn[0][i] = 1;		    /* für Dirichlet-RB */
      }

      nn[ndim+k][i] = i - ibase[k]; /* Nachbar x - e_k */
      if (ix[k] == 0){
	nn[ndim+k][i] += ibase[k+1];
	nn[0][i] = 1;		    /* für Dirichlet-RB */
      }
//       cout << "geomfkt: " << nn[k][i] << " " << nn[ndim+k][i] << endl;
    }

    for (k=1; k<=ndim; k++){        /* Koord. des naechsten Punkts */
      ix[k]++;
      if (ix[k] < lsize[k]) break;
      ix[k] = 0;
    }
  }
  free(ibase); free(ix);
}