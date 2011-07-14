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

#define N 10
#define ndim 2

void init_matrix();
inline vektor malloc_vektor(){ return new double[N*N]; }
inline void vec_copy(vektor a, vektor b);
inline void vec_addition(vektor a, vektor b, vektor c, int sign = 1);
void mult_matrix_vektor(matrix A, vektor x, vektor c);
inline void mult_skalar_vektor(double alpha, vektor x, vektor c);
inline double skalarprodukt(vektor a, vektor b);
void print_vektor(vektor x);
void print_matrix(matrix A);
void laplace(vektor x, double m2, vektor c);
double vektor_norm(vektor x);

int cg(vektor b, void (*fkt)(vektor x, double m2, vektor c), int max_it, double relerr = 1e-10);
void geom_pbc();

int **nn;
int lsize[ndim+1] = {0,N,N};
int nvol;
matrix A;

int main(int argc, char *argv[]) {

    TRandom3 *ran = new TRandom3(0);

    geom_pbc();
    init_matrix();
    double s[nvol];

    double dummy;
    for (int i=0; i<nvol; ++i){
      for (int j=0; j<nvol; ++j){
	    if (abs(i-j)<=1){
			dummy = ran->Uniform();
			A[i][j] = dummy;
			A[j][i] = dummy;
		}
		else A[i][j] = 0.;
      }
      s[i] = ran->Uniform();
    }

    int steps = cg(s, laplace, N+2, 1e-10);
    cout << "Anzahl der Schritte: " << steps << endl;

    return 0;
}

void init_matrix()
{
  A = (matrix)malloc(nvol*sizeof(vektor));
  A[0] = (vektor)malloc(nvol*nvol*sizeof(double));
  if(A[0] == NULL) {
    fprintf(stderr, "Not enough memory for allocating matrix\n");
    exit(1);
  }
  for(int i=1; i<nvol; ++i)
    A[i] = A[0] + nvol*i;
  return;
}

void vec_copy(vektor a, vektor b){
  #pragma omp parallel for
  for (int i=0; i<nvol; ++i){
    b[i] = a[i];
  }
  return;
}

void vec_addition(vektor a, vektor b, vektor c, int sign){      // sign: defaultparameter +1
  sign = (sign>0) ? 1 : -1;
  #pragma omp parallel for
  for (int i=0; i<nvol; ++i){
    c[i] = a[i] + sign*b[i];
  }
}

void mult_matrix_vektor(matrix A, vektor x, vektor c){
//   #pragma omp parallel for
  for (int i=0; i<nvol; ++i){
    c[i] = 0.;
    for (int j=0; j<nvol; ++j){
      c[i] += A[i][j]*x[j];
    }
  }
  return;
}

void laplace(vektor x, double m2, vektor c){
  double dummy = 0;
  int j = 0;
  #pragma omp parallel for private(j) reduction(-:dummy)
//   #pragma omp parallel for reduction() private() oder so
  for (int i=0; i<nvol; ++i){
    dummy = (2*ndim + m2)*x[i];
    for (j=1; j<=ndim; ++j){
      dummy -= x[nn[j][i]] + x[nn[ndim+j][i]];
    }
    c[i] = dummy;
  }
  // anzahl rechenschritte: N*N*(1 + ndim*2)
  return;
}

double vektor_norm(vektor x){
  return sqrt(skalarprodukt(x,x));
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

int cg(vektor s, void (*fkt)(vektor x, double m2, vektor c), int max_it, double relerr){
  double v[nvol];
  double q[nvol];
  double v_neu[nvol];
  double tol = relerr*relerr*skalarprodukt(s, s);
  double dummyvec[nvol];
  double dummyvec2[nvol];
  
  print_matrix(A);

  double s_norm = vektor_norm(s);
  if (vektor_norm(s) < sqrt(tol)) exit(0);

  mult_skalar_vektor(1/s_norm, s, dummyvec);
  vec_copy(dummyvec, v);
//   print_vektor(v);

  mult_matrix_vektor(A, v, q);
//   print_vektor(q);

  double alpha, beta;
  int counter = 0;

  for (int i=0; i<max_it; ++i){
    alpha = skalarprodukt(v, q);
//     cout << alpha << endl;

    mult_skalar_vektor(alpha, v, dummyvec);
    vec_addition(q, dummyvec, dummyvec2, -1);
    vec_copy(dummyvec2, q);

    beta = vektor_norm(q);
    cout << beta << endl;
    if (beta < tol) break;
    mult_skalar_vektor(1/beta, q, v_neu);

    mult_matrix_vektor(A, v_neu, dummyvec);
    mult_skalar_vektor(beta, v, dummyvec2);
    vec_addition(dummyvec, dummyvec2, q, -1);

    vec_copy(v_neu, v);

    ++counter;
  }

  return counter;
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
