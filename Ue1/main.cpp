#include <string>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <omp.h>

#include <TRandom3.h>

using std::string;
using std::cout;
using std::endl;

typedef double *vektor;
typedef double **matrix;

const int N = 10;

matrix malloc_matrix(int N);
inline vektor malloc_vektor(int N){ return (vektor)malloc(N*sizeof(double)); }
vektor vec_copy(vektor a);
vektor vec_addition(vektor a, vektor b, int sign = 1);
vektor mult_matrix_vektor(matrix A, vektor x);
vektor mult_skalar_vektor(double alpha, vektor x);
double skalarprodukt(vektor a, vektor b);
void print_vektor(vektor x);
void print_matrix(matrix A);

vektor cg(int N, matrix A, vektor x, vektor b, int max_it, double relerr = 1e-10, bool flag = 1);

int main(int argc, char *argv[]) {

    TRandom3 *ran = new TRandom3(0);

    vektor b = new double[N];
    vektor x = new double[N];
//     vektor b = malloc_vektor(N);
//     vektor x = malloc_vektor(N);

//     matrix A = new double[N][N];
    matrix A = malloc_matrix(N);
    double dummy = 0.;

    for (int i=0; i<N; ++i){
      for (int j=0; j<N; ++j){
        dummy = ran->Uniform();
        A[i][j] = dummy;
        A[j][i] = dummy;
      }
      b[i] = ran->Uniform();
      x[i] = ran->Uniform();
    }

//     print_matrix(A);
    vektor x2 = malloc_vektor(N);
    x2 = vec_copy(cg(N, A, x, b, 1000, 1e-10, 0));
//     print_vektor(x2);

    delete[] b;
    delete[] x;
//     delete[] A;
    return 0;
}

matrix malloc_matrix(int N)
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

vektor vec_copy(vektor a){
  vektor b = malloc_vektor(N);
  for (int i=0; i<N; ++i){
    b[i] = a[i];
  }
  return b;
}

vektor vec_addition(vektor a, vektor b, int sign){	// sign: defaultparameter +1
  vektor c = malloc_vektor(N);
  (sign>0) ? (sign = 1) : (sign = -1);
  for (int i=0; i<N; ++i){
    c[i] = a[i] + sign*b[i];
  }
  return c;
}

vektor mult_matrix_vektor(matrix A, vektor x){
  vektor c = malloc_vektor(N);
  #pragma omp parallel for
  for (int i=0; i<N; ++i){
    c[i] = 0.;
    for (int j=0; j<N; ++j){
      c[i] += A[i][j]*x[j];
    }
  }
  return c;
}

vektor mult_skalar_vektor(double alpha, vektor x){
  vektor c = malloc_vektor(N);
  #pragma omp parallel for
  for (int i=0; i<N; ++i){
    c[i] = alpha*x[i];
  }
  return c;
}

double skalarprodukt(vektor a, vektor b){
  double erg = 0;
  int i;
  #pragma omp parallel for reduction(+:erg)
  for (i=0; i<N; ++i){
    erg += a[i]*b[i];
  }
  return erg;
}

vektor cg(int N, matrix A, vektor x, vektor b, int max_it, double relerr, bool flag){
  vektor r;
  double tol = relerr*relerr*skalarprodukt(b, b);

  if (flag){
    r = vec_copy(b);
    for (int i=0; i<N; ++i){
      x[i] = 0.;
    }
  }
  else r = vec_addition(b, mult_matrix_vektor(A,x), -1);

  if (skalarprodukt(r, r) < tol) exit(0);

  vektor p = vec_copy(r);
  vektor s = malloc_vektor(N);
  double alpha, r2_alt, r2_neu, beta;
  int counter = 0;

  for (int i=0; i<max_it; ++i){
    s = mult_matrix_vektor(A, p);
    alpha = skalarprodukt(p, r) / skalarprodukt(p, s);
    x = vec_addition(x, mult_skalar_vektor(alpha, p));		// neues x_(k+1)
    r2_alt = skalarprodukt(r, r);				// r²_k; TODO kann man effizienter aus r2_neu holen
    r = vec_addition(r, mult_skalar_vektor(alpha, s), -1);	// neues r_(k+1)
    r2_neu = skalarprodukt(r, r);				// r²_(k+1)
    if (r2_neu < tol) break;
    beta = r2_neu / r2_alt;
    p = vec_addition(r, mult_skalar_vektor(beta, p));		// neues p_(k+1)
    ++counter;
  }

  cout << "Anzahl von Schritten: " << counter << endl;
  return x;
}

void print_vektor(vektor x){
  cout << "[";
  for (int i=0; i<N; ++i){
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