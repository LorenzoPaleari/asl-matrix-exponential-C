#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "mkl.h"

#ifdef FLOPS
double flops_opt = 0;
double flops_opt_N = 0;
double flops_opt_N2 = 0;
double flops_opt_N3 = 0;
#define FLOPS_RESET flops_opt = 0; flops_opt_N = 0; flops_opt_N2 = 0; flops_opt_N3 = 0;
#define FLOPS_ADD(x) flops_opt += x;
#define FLOPS_ADD_N(x) flops_opt_N += x;
#define FLOPS_ADD_N2(x) flops_opt_N2 += x;
#define FLOPS_ADD_N3(x) flops_opt_N3 += x;
#define FLOPS_PRINT printf("MAIN \nFlops: %f\n", flops_opt); printf("Flops N: %f\n", flops_opt_N); printf("Flops N^2: %f\n", flops_opt_N2); printf("Flops N^3: %f\n", flops_opt_N3);
#else
#define FLOPS_RESET
#define FLOPS_ADD(x)
#define FLOPS_ADD_N(x)
#define FLOPS_ADD_N2(x)
#define FLOPS_ADD_N3(x)
#define FLOPS_PRINT
#endif

double * I;

void replace_diag_exp(double *A, double *X, int n, int s) {
    double temp = pow(2, s);
    for (int i = 0; i < n; i++)
        X[i * n + i] = exp(A[i * n + i] * temp);
}

double sinch(double x, double y) {
    double temp = y - x;
    return (exp(temp * 0.5) - exp((x - y) * 0.5)) / temp;
}

void replace_superdiag_exp(double *A, double *X, int n, int s) {
    double lambda1;
    double lambda2;
    double t12;
    double temp = pow(2, s);

    FLOPS_ADD_N(17);
    for (int i = 0; i < n - 1; i++){
        lambda1 = A[i * n + i] * temp;
        lambda2 = A[(i + 1) * n + i + 1] * temp;
        t12 = A[i * n + i + 1] * temp;

        X[i * n + i + 1] = t12 * exp((lambda1 + lambda2) * 0.5) * sinch(lambda1, lambda2);
    }
}

void fragment_2_1(double *A, double **X, int n, int s){
    replace_diag_exp(A, *X, n, 0);
    FLOPS_ADD_N(2)
    FLOPS_ADD(1)

    for (int i = 1; i <= s; i++){
        // matrix_multiply_blas(n, n, n, *X, *X, X);
        replace_diag_exp(A, *X, n, i);
        // FLOPS_ADD_N3(2)
        FLOPS_ADD_N(2)
        FLOPS_ADD(1)
        replace_superdiag_exp(A, *X, n, i);
    }
}

int ell(double *A, int n, int m, double u){
    double norm_A, alpha;
    double *abs_A = allocate_matrix(n, n);
    FLOPS_ADD(5)
    int res;

    //1-norm of A
    norm_A = norm_1(A, n);
    FLOPS_ADD_N2(2)
    FLOPS_ADD_N(1)
    FLOPS_ADD(1)

    //absolut value of matrix A: all elements are positive
    m_abs(A, abs_A, n);
    FLOPS_ADD_N2(1)

    //Computation of alpha using the 1-norm estimation of abs_A and the coefficient c(m)
    alpha = _norm1est(abs_A, n, 2*m + 1) /norm_A;
    res = fmax(ceil(log2(alpha/(u * c[m])) / (2.0*m)), 0);
    FLOPS_ADD(8);

    free(abs_A);
    return res;
}

void pq_m_13(double *A, double *A2, double *A4, double *A6, int n, int m, double *P_m, double *Q_m){
    double *U1 = allocate_matrix(n, n);
    double *V1 = allocate_matrix(n, n);
    double *U2 = allocate_matrix(n, n);
    double *V2 = allocate_matrix(n, n);
    double *neg_U = allocate_matrix(n, n);
    double *temp = allocate_matrix(n, n);
    FLOPS_ADD(30)

    mcm(I, U2, b13[1], n);
    mcm(I, V2, b13[0], n);
    FLOPS_ADD_N2(2)

    // First part that has to be multiplied by A6
    mcm(A2, temp, b13[9], n);
    mcm(A4, U1, b13[11], n);
    m_sum(U1, temp, U1, n);
    FLOPS_ADD_N2(3)

    mcm(A2, temp, b13[8], n);
    mcm(A4, V1, b13[10], n);
    m_sum(V1, temp, V1, n);
    FLOPS_ADD_N2(3)

    mcm(A6, temp, b13[13], n);
    m_sum(U1, temp, U1, n);
    mcm(A6, temp, b13[12], n);
    m_sum(V1, temp, V1, n);
    FLOPS_ADD_N2(4)

    // Second part
    mcm(A2, temp, b13[3], n);
    m_sum(U2, temp, U2, n);
    mcm(A2, temp, b13[2], n);
    m_sum(V2, temp, V2, n);
    FLOPS_ADD_N2(4)

    mcm(A4, temp, b13[5], n);
    m_sum(U2, temp, U2, n);
    mcm(A4, temp, b13[4], n);
    m_sum(V2, temp, V2, n);
    FLOPS_ADD_N2(4)

    mcm(A6, temp, b13[7], n);
    m_sum(U2, temp, U2, n);
    mcm(A6, temp, b13[6], n);
    m_sum(V2, temp, V2, n);
    FLOPS_ADD_N2(4)

    // matrix_multiply_blas(n, n, n, U1, A6, &U1);
    // matrix_multiply_blas(n, n, n, V1, A6, &V1);
    // FLOPS_ADD_N3(4)
    // FLOPS_ADD(10)

    m_sum(U1, U2, U1, n);
    // matrix_multiply_blas(n, n, n, A, U1, &U1);
    m_sum(V1, V2, V1, n);
    FLOPS_ADD_N2(2)
    // FLOPS_ADD(5)
    // FLOPS_ADD_N3(2)

    m_sum(U1, V1, P_m, n);

    mcm(U1, neg_U, -1.0, n);
    m_sum(neg_U, V1, Q_m, n);
    FLOPS_ADD_N2(3)

    free(U1);
    free(V1);
    free(U2);
    free(V2);
    free(neg_U);
    free(temp);
}


void pq_m_3579(double *A, double *A2, double *A4, double *A6, double*A8, int n, int m, double *P_m, double *Q_m){
    double *U = allocate_matrix(n, n);
    double *V = allocate_matrix(n, n);
    double *neg_U = allocate_matrix(n, n);
    double *temp = allocate_matrix(n, n);
    FLOPS_ADD(20)

    if (m == 3){
        FLOPS_ADD(1)
        mcm(I, U, b3[1], n);
        mcm(I, V, b3[0], n);
        FLOPS_ADD_N2(2)

        mcm(A2, temp, b3[2+1], n);
        m_sum(U, temp, U, n);
        mcm(A2, temp, b3[2], n);
        m_sum(V, temp, V, n);
        FLOPS_ADD_N2(4)
        
        // matrix_multiply_blas(n, n, n, A, U, &U);
        // FLOPS_ADD_N3(2)
        // FLOPS_ADD(5)

        m_sum(U, V, P_m, n);

        mcm(U, neg_U, -1.0, n);
        m_sum(neg_U, V, Q_m, n);
        FLOPS_ADD_N2(3)
    }
    else if (m == 5){
        FLOPS_ADD(1)
        mcm(I, U, b5[1], n);
        mcm(I, V, b5[0], n);
        FLOPS_ADD_N2(2)

        mcm(A2, temp, b5[2+1], n);
        m_sum(U, temp, U, n);
        mcm(A2, temp, b5[2], n);
        m_sum(V, temp, V, n);
        FLOPS_ADD_N2(4)

        mcm(A4, temp, b5[4+1], n);
        m_sum(U, temp, U, n);
        mcm(A4, temp, b5[4], n);
        m_sum(V, temp, V, n);
        FLOPS_ADD_N2(4)

        //matrix_multiply_blas(n, n, n, A, U, &U);
        //FLOPS_ADD_N3(2)
        //FLOPS_ADD(5)

        m_sum(U, V, P_m, n);

        mcm(U, neg_U, -1.0, n);
        m_sum(neg_U, V, Q_m, n);
        FLOPS_ADD_N2(3)
    }
    else if (m == 7){
        FLOPS_ADD(1)
        mcm(I, U, b7[1], n);
        mcm(I, V, b7[0], n);
        FLOPS_ADD_N2(2)

        mcm(A2, temp, b7[2+1], n);
        m_sum(U, temp, U, n);
        mcm(A2, temp, b7[2], n);
        m_sum(V, temp, V, n);
        FLOPS_ADD_N2(4)

        mcm(A4, temp, b7[4+1], n);
        m_sum(U, temp, U, n);
        mcm(A4, temp, b7[4], n);
        m_sum(V, temp, V, n);
        FLOPS_ADD_N2(4)

        mcm(A6, temp, b7[6+1], n);
        m_sum(U, temp, U, n);
        mcm(A6, temp, b7[6], n);
        m_sum(V, temp, V, n);
        FLOPS_ADD_N2(4)

        //matrix_multiply_blas(n, n, n, A, U, &U);
        //FLOPS_ADD_N3(2)
        //FLOPS_ADD(5)

        m_sum(U, V, P_m, n);

        mcm(U, neg_U, -1.0, n);
        m_sum(neg_U, V, Q_m, n);
        FLOPS_ADD_N2(3)
    }
    else if (m == 9){
        FLOPS_ADD(1)
        mcm(I, U, b9[1], n);
        mcm(I, V, b9[0], n);
        FLOPS_ADD_N2(2)

        mcm(A2, temp, b9[2+1], n);
        m_sum(U, temp, U, n);
        mcm(A2, temp, b9[2], n);
        m_sum(V, temp, V, n);
        FLOPS_ADD_N2(4)

        mcm(A4, temp, b9[4+1], n);
        m_sum(U, temp, U, n);
        mcm(A4, temp, b9[4], n);
        m_sum(V, temp, V, n);
        FLOPS_ADD_N2(4)

        mcm(A6, temp, b9[6+1], n);
        m_sum(U, temp, U, n);
        mcm(A6, temp, b9[6], n);
        m_sum(V, temp, V, n);
        FLOPS_ADD_N2(4)

        mcm(A8, temp, b9[8+1], n);
        m_sum(U, temp, U, n);
        mcm(A8, temp, b9[8], n);
        m_sum(V, temp, V, n);
        FLOPS_ADD_N2(4)

        //matrix_multiply_blas(n, n, n, A, U, &U);
        //FLOPS_ADD_N3(2)
        //FLOPS_ADD(5)

        m_sum(U, V, P_m, n);

        mcm(U, neg_U, -1.0, n);
        m_sum(neg_U, V, Q_m, n);
        FLOPS_ADD_N2(3)
    } else {
        exit(EXIT_FAILURE);
    }
    free(U);
    free(V);
    free(neg_U);
    free(temp);
}

int mexp_no_blas(double *A, double *A2, double *A4, double *A6, double*A8, double *A10, double **X, int n){
    FLOPS_RESET
    double d4, d6, d8, d10;
    double eta1, eta3, eta4, eta5;
    double d6_factor = 1.0/6.0;
    FLOPS_ADD(2)
    //It is the "unit roundoff" of IEEE double precision arithmetic.
    double u = pow(2, -53);
    double s;
    double *temp1 = allocate_matrix(n, n);
    double *temp2 = allocate_matrix(n, n);
    FLOPS_ADD(15)
    I = identity_matrix(n);
    // A2 = NULL;
    // A4 = NULL;
    // A6 = NULL;
    // A8 = NULL;

    //matrix_multiply_blas(n, n, n, A, A, &A2);
    //FLOPS_ADD_N3(2)
    d6 = pow(_norm1est(A2, n, 3), d6_factor);
    d4 = pow(_norm1est(A2, n, 2), 0.25);
    eta1 = fmax(d4, d6);
    FLOPS_ADD(3)

    FLOPS_ADD(2)
    if (eta1 <= theta[3] && ell(A, n, 3, u) == 0) {
        pq_m_3579(A, A2, A4, A6, A8, n, 3, temp1, temp2);

        solve_system(temp1, temp2, X, n);

        free(temp1);
        free(temp2);
        
        free(I);
        FLOPS_PRINT
        return 1;
    }
    
    // matrix_multiply_blas(n, n, n, A2, A2, &A4);
    // FLOPS_ADD_N3(2)
    
    FLOPS_ADD(2)
    if (eta1 <= theta[5] && ell(A, n, 5, u) == 0){
        pq_m_3579(A, A2, A4, A6, A8, n, 5, temp1, temp2);

        solve_system(temp1, temp2, X, n);

        free(temp1);
        free(temp2);
        
        
        free(I);
        FLOPS_PRINT
        return 2;

    }

    // matrix_multiply_blas(n, n, n, A2, A4, &A6);
    // matrix_multiply_blas(n, n, n, A2, A6, &A8);
    // FLOPS_ADD_N3(4)
    d8 = pow(_norm1est(A4, n, 2), 0.125);
    eta3 = fmax(d6, d8);
    FLOPS_ADD(2)
    for (int m = 7; m < 10; m+=2) {
        if (eta3 <= theta[m] && ell(A, n, m, u) == 0){
            FLOPS_ADD(2)
            pq_m_3579(A, A2, A4, A6, A8, n, m, temp1, temp2);

            solve_system(temp1, temp2, X, n);
            free(temp1);
            free(temp2);
            
            
            
            
            free(I);
            FLOPS_PRINT
            return 3;

        }
    }
   
   // Up to line 27
   // Continues below
    // matrix_multiply_blas(n, n, n, A4, A6, &temp1);
    // FLOPS_ADD_N3(2)
    // FLOPS_ADD(5)
    d10 = pow(_norm1est(A10, n, 1), 0.1);
    eta4 = fmax(d8, d10);
    eta5 = fmin(eta3, eta4);
    FLOPS_ADD(3)

    s = fmax(ceil(log2(eta5 / theta[13])), 0);
    FLOPS_ADD(4)

    mcm(A, temp1, pow(2, -s), n);
    FLOPS_ADD_N2(1)
    FLOPS_ADD(2)
    s += ell(temp1, n, 13, u);

    //line 31 - end
    double *A_copy = allocate_matrix(n, n);
    mcm(A, A_copy, pow(2, -s), n);
    mcm(A2, A2, pow(4, -s), n);
    mcm(A4, A4, pow(16, -s), n);
    mcm(A6, A6, pow(64, -s), n);
    FLOPS_ADD_N2(4)
    FLOPS_ADD(9)

    pq_m_13(A_copy, A2, A4, A6, n, 13, temp1, temp2);
    
    solve_system(temp1, temp2, X, n);

    if (s==0){ 
        free(temp1);
        free(temp2);
        
        
        
        
        free(I);
        free(A_copy);
        FLOPS_PRINT
        return 4;
    }


    free(temp1);
    free(temp2);

    FLOPS_ADD_N2(0.5)
    if (is_triangular(A_copy, n)) {
        fragment_2_1(A_copy, X, n, s);
        free(I);
        free(A_copy);
        FLOPS_PRINT
        return 5;
    } else {
        // for (int i = 0; i < s; i++) {
        //     matrix_multiply_blas(n, n, n, *X, *X, X);
        //     FLOPS_ADD_N3(2)
        //     FLOPS_ADD(5)
        // }
        free(I);
        free(A_copy);
        FLOPS_PRINT
        return 6;
    }
}
