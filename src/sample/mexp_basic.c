#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "mkl.h"

double * I;

void replace_diag_exp(double *A, double *X, int n, int s) {
    for (int i = 0; i < n; i++)
        X[i * n + i] = exp(A[i * n + i] * pow(2, s));
}

double sinch(double x, double y) {
    return (exp((y - x) / 2.0) - exp((x - y) / 2.0)) / (y - x);
}

void replace_superdiag_exp(double *A, double *X, int n, int s) {
    double lambda1;
    double lambda2;
    double t12;

    for (int i = 0; i < n - 1; i++){
        lambda1 = A[i * n + i] * pow(2, s);
        lambda2 = A[(i + 1) * n + i + 1] * pow(2, s);
        t12 = A[i * n + i + 1] * pow(2, s);

        X[i * n + i + 1] = t12 * exp((lambda1 + lambda2) / 2.0) * sinch(lambda1, lambda2);
    }
}

void fragment_2_1(double *A, double **X, int n, int s){
    replace_diag_exp(A, *X, n, 0);

    for (int i = 1; i <= s; i++){
        matrix_multiply_blas(n, n, n, *X, *X, X);
        replace_diag_exp(A, *X, n, i);
        replace_superdiag_exp(A, *X, n, i);
    }
}

int ell(double *A, int n, int m){
    double u, norm_A, alpha;
    double *abs_A = allocate_matrix(n, n);
    int res;

    //It is the "unit roundoff" of IEEE double precision arithmetic.
    u = pow(2, -53);

    //1-norm of A
    norm_A = norm_1(A, n);

    //absolut value of matrix A: all elements are positive
    m_abs(A, abs_A, n);

    //Computation of alpha using the 1-norm estimation of abs_A and the coefficient c(m)
    alpha = _norm1est(abs_A, n, 2*m + 1) /norm_A;
    res = fmax(ceil(log2(alpha/(u * c[m])) / (2.0*m)), 0);

    free(abs_A);
    return res;
}

void pq_m_13(double *A, int n, int m, double *P_m, double *Q_m){
    double *U1 = allocate_matrix(n, n);
    double *V1 = allocate_matrix(n, n);
    double *U2 = allocate_matrix(n, n);
    double *V2 = allocate_matrix(n, n);
    double *neg_U = allocate_matrix(n, n);
    double *temp = allocate_matrix(n, n);

    mcm(I, U2, b13[1], n);
    mcm(I, V2, b13[0], n);

    // First part that has to be multiplied by A6
    mcm(A2, temp, b13[9], n);
    mcm(A4, U1, b13[11], n);
    m_sum(U1, temp, U1, n);

    mcm(A2, temp, b13[8], n);
    mcm(A4, V1, b13[10], n);
    m_sum(V1, temp, V1, n);

    mcm(A6, temp, b13[13], n);
    m_sum(U1, temp, U1, n);
    mcm(A6, temp, b13[12], n);
    m_sum(V1, temp, V1, n);

    // Second part
    mcm(A2, temp, b13[3], n);
    m_sum(U2, temp, U2, n);
    mcm(A2, temp, b13[2], n);
    m_sum(V2, temp, V2, n);

    mcm(A4, temp, b13[5], n);
    m_sum(U2, temp, U2, n);
    mcm(A4, temp, b13[4], n);
    m_sum(V2, temp, V2, n);

    mcm(A6, temp, b13[7], n);
    m_sum(U2, temp, U2, n);
    mcm(A6, temp, b13[6], n);
    m_sum(V2, temp, V2, n);

    matrix_multiply_blas(n, n, n, U1, A6, &U1);
    matrix_multiply_blas(n, n, n, V1, A6, &V1);

    m_sum(U1, U2, U1, n);
    matrix_multiply_blas(n, n, n, A, U1, &U1);
    m_sum(V1, V2, V1, n);

    m_sum(U1, V1, P_m, n);

    mcm(U1, neg_U, -1.0, n);
    m_sum(neg_U, V1, Q_m, n);

    free(U1);
    free(V1);
    free(U2);
    free(V2);
    free(neg_U);
    free(temp);
}


void pq_m_3579(double *A, int n, int m, double *P_m, double *Q_m){
    double *U = allocate_matrix(n, n);
    double *V = allocate_matrix(n, n);
    double *neg_U = allocate_matrix(n, n);
    double *temp = allocate_matrix(n, n);

    if (m == 3){
        mcm(I, U, b3[1], n);
        mcm(I, V, b3[0], n);

        mcm(A2, temp, b3[2+1], n);
        m_sum(U, temp, U, n);
        mcm(A2, temp, b3[2], n);
        m_sum(V, temp, V, n);
        
        matrix_multiply_blas(n, n, n, A, U, &U);

        m_sum(U, V, P_m, n);

        mcm(U, neg_U, -1.0, n);
        m_sum(neg_U, V, Q_m, n);
    }
    else if (m == 5){
        mcm(I, U, b5[1], n);
        mcm(I, V, b5[0], n);

        mcm(A2, temp, b5[2+1], n);
        m_sum(U, temp, U, n);
        mcm(A2, temp, b5[2], n);
        m_sum(V, temp, V, n);

        mcm(A4, temp, b5[4+1], n);
        m_sum(U, temp, U, n);
        mcm(A4, temp, b5[4], n);
        m_sum(V, temp, V, n);

        matrix_multiply_blas(n, n, n, A, U, &U);

        m_sum(U, V, P_m, n);

        mcm(U, neg_U, -1.0, n);
        m_sum(neg_U, V, Q_m, n);
    }
    else if (m == 7){
        mcm(I, U, b7[1], n);
        mcm(I, V, b7[0], n);

        mcm(A2, temp, b7[2+1], n);
        m_sum(U, temp, U, n);
        mcm(A2, temp, b7[2], n);
        m_sum(V, temp, V, n);

        mcm(A4, temp, b7[4+1], n);
        m_sum(U, temp, U, n);
        mcm(A4, temp, b7[4], n);
        m_sum(V, temp, V, n);

        mcm(A6, temp, b7[6+1], n);
        m_sum(U, temp, U, n);
        mcm(A6, temp, b7[6], n);
        m_sum(V, temp, V, n);

        matrix_multiply_blas(n, n, n, A, U, &U);

        m_sum(U, V, P_m, n);

        mcm(U, neg_U, -1.0, n);
        m_sum(neg_U, V, Q_m, n);
    }
    else if (m == 9){
        mcm(I, U, b9[1], n);
        mcm(I, V, b9[0], n);

        mcm(A2, temp, b9[2+1], n);
        m_sum(U, temp, U, n);
        mcm(A2, temp, b9[2], n);
        m_sum(V, temp, V, n);

        mcm(A4, temp, b9[4+1], n);
        m_sum(U, temp, U, n);
        mcm(A4, temp, b9[4], n);
        m_sum(V, temp, V, n);

        mcm(A6, temp, b9[6+1], n);
        m_sum(U, temp, U, n);
        mcm(A6, temp, b9[6], n);
        m_sum(V, temp, V, n);

        mcm(A8, temp, b9[8+1], n);
        m_sum(U, temp, U, n);
        mcm(A8, temp, b9[8], n);
        m_sum(V, temp, V, n);

        matrix_multiply_blas(n, n, n, A, U, &U);

        m_sum(U, V, P_m, n);

        mcm(U, neg_U, -1.0, n);
        m_sum(neg_U, V, Q_m, n);
    } else {
        exit(EXIT_FAILURE);
    }
    free(U);
    free(V);
    free(neg_U);
    free(temp);
}

int mexp(double *A, double **X, int n){
    I = identity_matrix(n);
    A2 = NULL;
    A4 = NULL;
    A6 = NULL;
    A8 = NULL;

    matrix_multiply_blas(n, n, n, A, A, &A2);
    double d6 = pow(_norm1est(A2, n, 3), 1.0/6.0);
    double eta1 = fmax(pow(_norm1est(A2, n, 2), 1.0/4.0), d6);

    if (eta1 <= theta[3] && ell(A, n, 3) == 0) {
        double* p3 = allocate_matrix(n, n);;
        double* q3 = allocate_matrix(n, n);

        pq_m_3579(A, n, 3, p3, q3);

        solve_system(p3, q3, X, n);

        free(p3);
        free(q3);
        free(A2);
        free(I);
        return 1;
    }
    

    matrix_multiply_blas(n, n, n, A2, A2, &A4);
    double d4 = pow(norm_1(A4, n), 1.0 / 4.0);
    double eta2 = fmax(d4, d6);

    if (eta2 <= theta[5] && ell(A, n, 5) == 0){
        double* p5 = allocate_matrix(n, n);
        double* q5 = allocate_matrix(n, n);

        pq_m_3579(A, n, 5, p5, q5);

        solve_system(p5, q5, X, n);

        free(p5);
        free(q5);
        free(A2);
        free(A4);
        free(I);
        return 2;

    }

    matrix_multiply_blas(n, n, n, A2, A4, &A6);
    matrix_multiply_blas(n, n, n, A2, A6, &A8);
    d6 = pow(norm_1(A6, n), 1.0 / 6.0);
    double d8 = pow(_norm1est(A4, n, 2), 1.0 / 8.0);
    double eta3 = fmax(d6, d8);

    for (int m = 7; m < 10; m+=2) {
        if (eta3 <= theta[m] && ell(A, n, m) == 0){
            double* p_m = allocate_matrix(n, n);
            double* q_m = allocate_matrix(n, n);

            pq_m_3579(A, n, m, p_m, q_m);

            solve_system(p_m, q_m, X, n);
            free(p_m);
            free(q_m);
            free(A2);
            free(A4);
            free(A6);
            free(A8);
            free(I);
            return 3;

        }
    }
   
   // Up to line 27
   // Continues below
    double *A_temp = allocate_matrix(n, n);
    int s;

    matrix_multiply_blas(n, n, n, A4, A6, &A_temp);
    double eta4 = fmax(d8, pow(_norm1est(A_temp, n, 1), 1.0 / 10.0));
    double eta5 = fmin(eta3, eta4);

    s = fmax(ceil(log2(eta5 / theta[13])), 0);

    mcm(A, A_temp, pow(2.0, -s), n);
    s += ell(A_temp, n, 13);

    //line 31 - end
    double *A_copy = allocate_matrix(n, n);
    mcm(A, A_copy, pow(2.0, -s), n);
    mcm(A2, A2, pow(2.0, -2*s), n);
    mcm(A4, A4, pow(2.0, -4*s), n);
    mcm(A6, A6, pow(2.0, -6*s), n);

    double* p13 = allocate_matrix(n, n);
    double* q13 = allocate_matrix(n, n);

    pq_m_13(A_copy, n, 13, p13, q13);
    
    if (s==0){
        solve_system(p13, q13, X, n);
        free(p13);
        free(q13);
        free(A2);
        free(A4);
        free(A6);
        free(A8);
        free(I);
        free(A_temp);
        free(A_copy);
        return 4;
    }

    solve_system(p13, q13, X, n);


    free(p13);
    free(q13);
    free(A_temp);


    if (is_triangular(A_copy, n)) {
        fragment_2_1(A_copy, X, n, s);
        free(I);
        free(A2);
        free(A4);
        free(A6);
        free(A8);
        free(A_copy);
        return 5;
    } else {
        for (int i = 0; i < s; i++) {
            matrix_multiply_blas(n, n, n, *X, *X, X);
        }
        free(I);
        free(A2);
        free(A4);
        free(A6);
        free(A8);
        free(A_copy);
        return 6;
    }
}
