/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
       int n = *la;
    int m = *lab;   

    for (int j = 0; j < n; j++)
    {
        // Colonne 0 = 0 partout
        AB[j*m + 0] = 0.0;

        // Colonne 1 = -1 sauf première ligne = 0
        AB[j*m + 1] = (j == 0) ? 0.0 : -1.0;

        // Colonne 2 = 2 partout
        AB[j*m + 2] = 2.0;

        // Colonne 3 = -1 sauf dernière ligne = 0
        AB[j*m + 3] = (j == n-1) ? 0.0 : -1.0;
    }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the identity matrix
  // Only the main diagonal should have 1, all other entries are 0
    int i, j;
    int n = *la;      // nombre de points
    int kd = *kv;     // demi-largeur de bande
    int ldab = *lab;  // dimension principale du tableau AB

    // Initialisation à zéro
    for (i = 0; i < ldab * n; i++) {
        AB[i] = 0.0;
    }

    // Remplissage de la diagonale principale avec des 1
    for (j = 0; j < n; j++) {
        AB[kd + j*ldab] = 1.0;  // seule la diagonale principale vaut 1
    }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // TODO: Compute RHS vector
   int n = *la;

    // Initialiser le RHS à zéro
    for (int i = 0; i < n; i++) {
        RHS[i] = 0.0;
    }

    // Ajouter les contributions des conditions aux limites
    if (n > 0) {
        RHS[0]   += *BC0;      // première ligne (bord gauche)
        RHS[n-1] += *BC1;      // dernière ligne (bord droit)
    }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D
  int n = *la;

    for (int i = 0; i < n; i++) {
        // solution linéaire entre BC0 et BC1
        EX_SOL[i] = *BC0 + ( *BC1 - *BC0 ) * X[i];
    }
}  

void set_grid_points_1D(double* x, int* la){
  // TODO: Generate uniformly spaced grid points in [0,1]
   int n = *la;           // nombre de points intérieurs
    double h = 1.0 / (n + 1);  // pas de discrétisation

    for (int i = 0; i < n; i++) {
        x[i] = (i + 1) * h;  // points intérieurs : x_1 = h, x_2 = 2h, ..., x_n = n*h
    }
}

double relative_forward_error(double* x, double* y, int* la){
  // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)
  //return 0.0;
  int n = *la;
    double num = 0.0;   // norme de x - y
    double denom = 0.0; // norme de y

    for (int i = 0; i < n; i++) {
        double diff = x[i] - y[i];
        num   += diff * diff;
        denom += y[i] * y[i];
    }

    num = sqrt(num);
    denom = sqrt(denom);

    if (denom == 0.0) return 0.0; // éviter division par zéro
    return num / denom;
}

int indexABCol(int i, int j, int *lab){
  // TODO: Return the correct index formula for column-major band storage
 // return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices
}


//docker build -f docker/Dockerfile --progress plain -t tp-cn:latest .
//docker run --rm -it tp-cn:latest
//cat AB.dat

