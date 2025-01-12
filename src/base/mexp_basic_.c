#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "utils_.h"
#include "mkl.h"

#ifdef FLOPS
double flops_opt___ = 0;
double flops_opt_N___ = 0;
double flops_opt_N2___ = 0;
double flops_opt_N3___ = 0;
#define FLOPS_RESET___ flops_opt___ = 0; flops_opt_N___ = 0; flops_opt_N2___ = 0; flops_opt_N3___ = 0;
#define FLOPS_ADD___(x) flops_opt___ += x;
#define FLOPS_ADD_N___(x) flops_opt_N___ += x;
#define FLOPS_ADD_N2___(x) flops_opt_N2___ += x;
#define FLOPS_ADD_N3___(x) flops_opt_N3___ += x;
#define FLOPS_PRINT___ printf("BASE\n"); printf("MAIN \nFlops: %f\n", flops_opt___); printf("Flops N: %f\n", flops_opt_N___); printf("Flops N^2: %f\n", flops_opt_N2___); printf("Flops N^3: %f\n", flops_opt_N3___);
#else
#define FLOPS_RESET___
#define FLOPS_ADD___(x)
#define FLOPS_ADD_N___(x)
#define FLOPS_ADD_N2___(x)
#define FLOPS_ADD_N3___(x)
#define FLOPS_PRINT___
#endif

double * I_;

void replace_diag_exp_(double *A, double *X, int n, int s) {
    for (int i = 0; i < n; i++)
        X[i * n + i] = exp(A[i * n + i] * pow(2, s));
}

double sinch_(double x, double y) {
    return (exp((y - x) / 2.0) - exp((x - y) / 2.0)) / (y - x);
}

void replace_superdiag_exp_(double *A, double *X, int n, int s) {
    double lambda1;
    double lambda2;
    double t12;

    FLOPS_ADD_N___(20);
    for (int i = 0; i < n - 1; i++){
        lambda1 = A[i * n + i] * pow(2, s);
        lambda2 = A[(i + 1) * n + i + 1] * pow(2, s);
        t12 = A[i * n + i + 1] * pow(2, s);

        X[i * n + i + 1] = t12 * exp((lambda1 + lambda2) / 2.0) * sinch_(lambda1, lambda2);
    }
}

void fragment_2_1_(double *A, double **X, int n, int s){
    replace_diag_exp_(A, *X, n, 0);
    FLOPS_ADD_N___(3);

    for (int i = 1; i <= s; i++){
        matrix_multiply_blas_(n, n, n, *X, *X, X);
        FLOPS_ADD_N3___(2)
        replace_diag_exp_(A, *X, n, i);
        FLOPS_ADD_N___(3);
        replace_superdiag_exp_(A, *X, n, i);
    }
}

int ell_(double *A, int n, int m){
    double u, norm_A, alpha;
    double *abs_A = allocate_matrix_(n, n);
    int res;

    //It is the "unit roundoff" of IEEE double precision arithmetic.
    u = pow(2, -53);
    FLOPS_ADD___(1);
    //1-norm of A
    norm_A = norm_1_(A, n);
    FLOPS_ADD_N2___(2)
    FLOPS_ADD_N___(1)

    //absolut value of matrix A: all elements are positive
    m_abs_(A, abs_A, n);
    FLOPS_ADD_N2___(1)

    //Computation of alpha using the 1-norm estimation of abs_A and the coefficient c(m)
    alpha = _norm1est_(abs_A, n, 2*m + 1) /norm_A;
    res = fmax(ceil(log2(alpha/(u * c_[m])) / (2.0*m)), 0);
    FLOPS_ADD___(8);

    free(abs_A);
    return res;
}

void pq_m_13_(double *A, int n, int m, double *P_m, double *Q_m){
    double *U1 = allocate_matrix_(n, n);
    double *V1 = allocate_matrix_(n, n);
    double *U2 = allocate_matrix_(n, n);
    double *V2 = allocate_matrix_(n, n);
    double *neg_U = allocate_matrix_(n, n);
    double *temp = allocate_matrix_(n, n);

    mcm_(I_, U2, b13_[1], n);
    mcm_(I_, V2, b13_[0], n);
    FLOPS_ADD_N2___(2)

    // First part that has to be multiplied by A6_
    mcm_(A2_, temp, b13_[9], n);
    mcm_(A4_, U1, b13_[11], n);
    m_sum_(U1, temp, U1, n);
    FLOPS_ADD_N2___(3)

    mcm_(A2_, temp, b13_[8], n);
    mcm_(A4_, V1, b13_[10], n);
    m_sum_(V1, temp, V1, n);
    FLOPS_ADD_N2___(3)

    mcm_(A6_, temp, b13_[13], n);
    m_sum_(U1, temp, U1, n);
    mcm_(A6_, temp, b13_[12], n);
    m_sum_(V1, temp, V1, n);
    FLOPS_ADD_N2___(4)

    // Second part
    mcm_(A2_, temp, b13_[3], n);
    m_sum_(U2, temp, U2, n);
    mcm_(A2_, temp, b13_[2], n);
    m_sum_(V2, temp, V2, n);
    FLOPS_ADD_N2___(4)

    mcm_(A4_, temp, b13_[5], n);
    m_sum_(U2, temp, U2, n);
    mcm_(A4_, temp, b13_[4], n);
    m_sum_(V2, temp, V2, n);
    FLOPS_ADD_N2___(4)

    mcm_(A6_, temp, b13_[7], n);
    m_sum_(U2, temp, U2, n);
    mcm_(A6_, temp, b13_[6], n);
    m_sum_(V2, temp, V2, n);
    FLOPS_ADD_N2___(4)

    matrix_multiply_blas_(n, n, n, U1, A6_, &U1);
    matrix_multiply_blas_(n, n, n, V1, A6_, &V1);
    FLOPS_ADD_N3___(4)

    m_sum_(U1, U2, U1, n);
    matrix_multiply_blas_(n, n, n, A, U1, &U1);
    m_sum_(V1, V2, V1, n);
    FLOPS_ADD_N3___(2)
    FLOPS_ADD_N2___(2)

    m_sum_(U1, V1, P_m, n);

    mcm_(U1, neg_U, -1.0, n);
    m_sum_(neg_U, V1, Q_m, n);
    FLOPS_ADD_N2___(3)

    free(U1);
    free(V1);
    free(U2);
    free(V2);
    free(neg_U);
    free(temp);
}


void pq_m_3579_(double *A, int n, int m, double *P_m, double *Q_m){
    double *U = allocate_matrix_(n, n);
    double *V = allocate_matrix_(n, n);
    double *neg_U = allocate_matrix_(n, n);
    double *temp = allocate_matrix_(n, n);

    if (m == 3){
        FLOPS_ADD___(1)
        mcm_(I_, U, b3_[1], n);
        mcm_(I_, V, b3_[0], n);
        FLOPS_ADD_N2___(2)

        mcm_(A2_, temp, b3_[2+1], n);
        m_sum_(U, temp, U, n);
        mcm_(A2_, temp, b3_[2], n);
        m_sum_(V, temp, V, n);
        FLOPS_ADD_N2___(4)
        
        matrix_multiply_blas_(n, n, n, A, U, &U);
        FLOPS_ADD_N3___(2)

        m_sum_(U, V, P_m, n);

        mcm_(U, neg_U, -1.0, n);
        m_sum_(neg_U, V, Q_m, n);
        FLOPS_ADD_N2___(3)
    }
    else if (m == 5){
        FLOPS_ADD___(1)
        mcm_(I_, U, b5_[1], n);
        mcm_(I_, V, b5_[0], n);
        FLOPS_ADD_N2___(2)

        mcm_(A2_, temp, b5_[2+1], n);
        m_sum_(U, temp, U, n);
        mcm_(A2_, temp, b5_[2], n);
        m_sum_(V, temp, V, n);
        FLOPS_ADD_N2___(4)

        mcm_(A4_, temp, b5_[4+1], n);
        m_sum_(U, temp, U, n);
        mcm_(A4_, temp, b5_[4], n);
        m_sum_(V, temp, V, n);
        FLOPS_ADD_N2___(4)

        matrix_multiply_blas_(n, n, n, A, U, &U);
        FLOPS_ADD_N3___(2)

        m_sum_(U, V, P_m, n);

        mcm_(U, neg_U, -1.0, n);
        m_sum_(neg_U, V, Q_m, n);
        FLOPS_ADD_N2___(3)
    }
    else if (m == 7){
        mcm_(I_, U, b7_[1], n);
        mcm_(I_, V, b7_[0], n);
        FLOPS_ADD_N2___(2)

        mcm_(A2_, temp, b7_[2+1], n);
        m_sum_(U, temp, U, n);
        mcm_(A2_, temp, b7_[2], n);
        m_sum_(V, temp, V, n);
        FLOPS_ADD_N2___(4)

        mcm_(A4_, temp, b7_[4+1], n);
        m_sum_(U, temp, U, n);
        mcm_(A4_, temp, b7_[4], n);
        m_sum_(V, temp, V, n);
        FLOPS_ADD_N2___(4)

        mcm_(A6_, temp, b7_[6+1], n);
        m_sum_(U, temp, U, n);
        mcm_(A6_, temp, b7_[6], n);
        m_sum_(V, temp, V, n);
        FLOPS_ADD_N2___(4)

        matrix_multiply_blas_(n, n, n, A, U, &U);
        FLOPS_ADD_N3___(2)

        m_sum_(U, V, P_m, n);

        mcm_(U, neg_U, -1.0, n);
        m_sum_(neg_U, V, Q_m, n);
        FLOPS_ADD_N2___(3)
    }
    else if (m == 9){
        mcm_(I_, U, b9_[1], n);
        mcm_(I_, V, b9_[0], n);
        FLOPS_ADD_N2___(2)

        mcm_(A2_, temp, b9_[2+1], n);
        m_sum_(U, temp, U, n);
        mcm_(A2_, temp, b9_[2], n);
        m_sum_(V, temp, V, n);
        FLOPS_ADD_N2___(4)

        mcm_(A4_, temp, b9_[4+1], n);
        m_sum_(U, temp, U, n);
        mcm_(A4_, temp, b9_[4], n);
        m_sum_(V, temp, V, n);
        FLOPS_ADD_N2___(4)

        mcm_(A6_, temp, b9_[6+1], n);
        m_sum_(U, temp, U, n);
        mcm_(A6_, temp, b9_[6], n);
        m_sum_(V, temp, V, n);
        FLOPS_ADD_N2___(4)

        mcm_(A8_, temp, b9_[8+1], n);
        m_sum_(U, temp, U, n);
        mcm_(A8_, temp, b9_[8], n);
        m_sum_(V, temp, V, n);
        FLOPS_ADD_N2___(4)

        matrix_multiply_blas_(n, n, n, A, U, &U);
        FLOPS_ADD_N3___(2)

        m_sum_(U, V, P_m, n);

        mcm_(U, neg_U, -1.0, n);
        m_sum_(neg_U, V, Q_m, n);
        FLOPS_ADD_N2___(3)
    } else {
        exit(EXIT_FAILURE);
    }
    free(U);
    free(V);
    free(neg_U);
    free(temp);
}

int mexp_model(double *A, double **X, int n){
    FLOPS_RESET___
    I_ = identity_matrix_(n);
    A2_ = NULL;
    A4_ = NULL;
    A6_ = NULL;
    A8_ = NULL;

    matrix_multiply_blas_(n, n, n, A, A, &A2_);
    FLOPS_ADD_N3___(2)
    double d6 = pow(_norm1est_(A2_, n, 3), 1.0/6.0);
    double eta1 = fmax(pow(_norm1est_(A2_, n, 2), 1.0/4.0), d6);
    FLOPS_ADD___(5)

    FLOPS_ADD___(2)
    if (eta1 <= theta_[3] && ell_(A, n, 3) == 0) {
        double* p3 = allocate_matrix_(n, n);;
        double* q3 = allocate_matrix_(n, n);

        pq_m_3579_(A, n, 3, p3, q3);

        solve_system_(p3, q3, X, n);

        free(p3);
        free(q3);
        free(A2_);
        free(I_);
        FLOPS_PRINT___
        return 1;
    }
    

    matrix_multiply_blas_(n, n, n, A2_, A2_, &A4_);
    FLOPS_ADD_N3___(2)
    double d4 = pow(norm_1_(A4_, n), 1.0 / 4.0);
    double etA2_ = fmax(d4, d6);
    FLOPS_ADD___(4)

    if (etA2_ <= theta_[5] && ell_(A, n, 5) == 0){
        double* p5 = allocate_matrix_(n, n);
        double* q5 = allocate_matrix_(n, n);

        pq_m_3579_(A, n, 5, p5, q5);

        solve_system_(p5, q5, X, n);

        free(p5);
        free(q5);
        free(A2_);
        free(A4_);
        free(I_);
        FLOPS_PRINT___
        return 2;

    }

    matrix_multiply_blas_(n, n, n, A2_, A4_, &A6_);
    matrix_multiply_blas_(n, n, n, A2_, A6_, &A8_);
    FLOPS_ADD_N3___(4)
    d6 = pow(norm_1_(A6_, n), 1.0 / 6.0);
    double d8 = pow(_norm1est_(A4_, n, 2), 1.0 / 8.0);
    double eta3 = fmax(d6, d8);
    FLOPS_ADD___(5)

    for (int m = 7; m < 10; m+=2) {
        if (eta3 <= theta_[m] && ell_(A, n, m) == 0){
            FLOPS_ADD___(2)
            double* p_m = allocate_matrix_(n, n);
            double* q_m = allocate_matrix_(n, n);

            pq_m_3579_(A, n, m, p_m, q_m);

            solve_system_(p_m, q_m, X, n);
            free(p_m);
            free(q_m);
            free(A2_);
            free(A4_);
            free(A6_);
            free(A8_);
            free(I_);
            FLOPS_PRINT___
            return 3;

        }
    }
   
   // Up to line 27
   // Continues below
    double *A_temp = allocate_matrix_(n, n);
    int s;

    matrix_multiply_blas_(n, n, n, A4_, A6_, &A_temp);
    FLOPS_ADD_N3___(2)
    double eta4 = fmax(d8, pow(_norm1est_(A_temp, n, 1), 1.0 / 10.0));
    double eta5 = fmin(eta3, eta4);
    FLOPS_ADD___(4)

    s = fmax(ceil(log2(eta5 / theta_[13])), 0);
    FLOPS_ADD___(4)

    mcm_(A, A_temp, pow(2.0, -s), n);
    FLOPS_ADD___(1)
    FLOPS_ADD_N2___(1)
    s += ell_(A_temp, n, 13);
    FLOPS_ADD___(1)

    //line 31 - end
    double *A_copy = allocate_matrix_(n, n);
    mcm_(A, A_copy, pow(2.0, -s), n);
    mcm_(A2_, A2_, pow(2.0, -2*s), n);
    mcm_(A4_, A4_, pow(2.0, -4*s), n);
    mcm_(A6_, A6_, pow(2.0, -6*s), n);
    FLOPS_ADD_N2___(4)
    FLOPS_ADD___(4)

    double* p13 = allocate_matrix_(n, n);
    double* q13 = allocate_matrix_(n, n);

    pq_m_13_(A_copy, n, 13, p13, q13);
    
    if (s==0){
        solve_system_(p13, q13, X, n);
        free(p13);
        free(q13);
        free(A2_);
        free(A4_);
        free(A6_);
        free(A8_);
        free(I_);
        free(A_temp);
        free(A_copy);
        FLOPS_PRINT___
        return 4;
    }

    solve_system_(p13, q13, X, n);


    free(p13);
    free(q13);
    free(A_temp);

    FLOPS_ADD_N2___(0.5)
    if (is_triangular_(A_copy, n)) {
        fragment_2_1_(A_copy, X, n, s);
        free(I_);
        free(A2_);
        free(A4_);
        free(A6_);
        free(A8_);
        free(A_copy);
        FLOPS_PRINT___
        return 5;
    } else {
        for (int i = 0; i < s; i++) {
            matrix_multiply_blas_(n, n, n, *X, *X, X);
        }
        FLOPS_ADD_N3___(s*2)
        free(I_);
        free(A2_);
        free(A4_);
        free(A6_);
        free(A8_);
        free(A_copy);
        FLOPS_PRINT___
        return 6;
    }
}
