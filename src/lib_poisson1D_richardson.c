/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "lib_poisson1D.h"
#include <math.h>


void eig_poisson1D(double* eigval, int *la){
  // Compute all eigenvalues for the 1D Poisson operator
  int n = *la;
    for (int k = 1; k <= n; k++) {
        eigval[k-1] = 2.0 * (1.0 - cos(M_PI * k / (n + 1)));
    }
}

double eigmax_poisson1D(int *la){
  // Compute and return the maximum eigenvalue for the 1D Poisson operator
  int n = *la;
  return 2.0 * (1.0 - cos(M_PI * n / (n + 1)));
}

double eigmin_poisson1D(int *la){
  // Compute and return the minimum eigenvalue for the 1D Poisson operator
  int n = *la;
  return 2.0 * (1.0 - cos(M_PI / (n + 1)));
  
}

double richardson_alpha_opt(int *la){
  // Compute alpha_opt
  int n = *la;
  double lambda_max = eigmax_poisson1D(la);
  double lambda_min = eigmin_poisson1D(la);
  return 2.0 / (lambda_max + lambda_min);
}

/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iteration
  // 1. Compute residual r = b - A*x (use dgbmv for matrix-vector product)
  // 2. Update x = x + alpha*r (use daxpy)
  // 3. Check convergence: ||r||_2 < tol (use dnrm2)
  // 4. Store residual norm in resvec and repeat
    
    int n = *la;
    double alpha = *alpha_rich;

    // Variables pour BLAS/LAPACK
    int incx = 1;
    double one = 1.0;
    double minus_one = -1.0;

    double *r = (double *) malloc(sizeof(double) * n); // vecteur résiduel
    if (!r) {
        *nbite = 0;
        return;
    }

    // Calcul de la norme de b
    double norm_b = dnrm2_(&n, RHS, &incx);
    double norm_r;

    int iter;
    for (iter = 0; iter < *maxit; iter++) {

        // r = b
        for (int i = 0; i < n; i++) r[i] = RHS[i];

        // r = r - A*x  => r = RHS - A*X
        dgbmv_("N", la, la, kl, ku, &one, AB, lab, X, &incx, &minus_one, r, &incx);

        // Calcul de la norme relative du résidu
        norm_r = dnrm2_(&n, r, &incx);
        resvec[iter] = norm_r / norm_b;

        // Vérification convergence
        if (norm_r / norm_b < *tol) {
            iter++; // On compte l’itération finale
            break;
        }

        // Mise à jour x = x + alpha * r
        daxpy_(&n, &alpha, r, &incx, X, &incx);
    }

    *nbite = iter; // Nombre d'itérations effectuées
    free(r);

}

/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  //  Extract diagonal elements from AB and store in MB
  // MB should contain only the diagonal of A
    int n = *la;        // nombre de points
    int kd = *kv;       // demi-largeur de bande
    int ldab = *lab;    // dimension principale de AB

    for (int i = 0; i < n; i++){
        MB[i] = AB[kd + i*ldab];  // diagonale principale se trouve à ligne kv
    }
}
/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // Extract diagonal and lower diagonal from AB
  // MB should contain the lower triangular part (including diagonal) of A
    int n = *la;
    int ldab = *lab;
    int kd = *kv;

    for (int i = 0; i < n; i++){
        // Diagonale principale
        MB[i*2] = AB[kd + i*ldab];

        // Sous-diagonale (0 pour la première ligne)
        MB[i*2 + 1] = (i == 0) ? 0.0 : AB[kd+1 + (i-1)*ldab];
    }
}

/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  //Implement Richardson iterative method

    int n = *la;
    int ldab = *lab;
    int kd = *ku;  // demi-largeur superdiagonale pour col-major
    int incx = 1;
    double one = 1.0, minus_one = -1.0;

    // Vecteur résiduel
    double *r = (double *) malloc(sizeof(double)*n);
    if(!r) { *nbite=0; return; }

    // Norme initiale du RHS
    double norm_b = dnrm2_(&n, RHS, &incx);
    double norm_r;
    int iter;

    for(iter=0; iter<*maxit; iter++){
        // r = b - A*x
        for(int i=0;i<n;i++) r[i] = RHS[i];
        dgbmv_("N", la, la, kl, ku, &one, AB, lab, X, &incx, &minus_one, r, &incx);

        // Calcul de la norme relative
        norm_r = dnrm2_(&n, r, &incx);
        resvec[iter] = norm_r / norm_b;

        // Vérification convergence
        if(norm_r / norm_b < *tol){ iter++; break; }

        // Mise à jour x selon MB
        // Détecter si MB est Jacobi ou Gauss-Seidel
        // Jacobi : MB[i] contient seulement D[i] → sous-diagonale = 0
        // Gauss-Seidel : MB[i*2] = D[i], MB[i*2+1] = E[i]
        if(MB[1] == 0.0){ 
            // Jacobi
            for(int i=0;i<n;i++) X[i] += r[i] / MB[i];
        } else { 
            // Gauss-Seidel
            for(int i=0;i<n;i++){
                double sum = r[i];
                if(i>0) sum -= MB[(i-1)*2+1]*X[i-1]; // contribution sous-diagonale
                X[i] += sum / MB[i*2];                // diviser par la diagonale
            }
        }
    }

    *nbite = iter;  // nombre d'itérations effectuées
    free(r);
}

