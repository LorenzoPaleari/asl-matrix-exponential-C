#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

int mexp(double *A, double **X, int n){
    FLOPS_RESET
    double d4, d6, d8, d10;
    double eta1, eta3, eta4, eta5;
    int s;
    double s_temp;
    double normA, alpha;
    double ell;
    double d6_factor = 1.0/6.0;
    //It is the "unit roundoff" of IEEE double precision arithmetic.
    double u = pow(2, -53);
    FLOPS_ADD(22)
    double *temp1 = allocate_matrix(n, n);
    double *temp2 = allocate_matrix(n, n);
    //double *temp3 = allocate_matrix(n, n);
    double *temp4 = allocate_matrix(n, n);
    double *tempA = allocate_matrix(n, n);
    A2 = NULL;
    A4 = NULL;
    A6 = NULL;
    A8 = NULL;

    normA = norm_1(A, n);   //normA used in ell
    m_abs(A, tempA, n);   //tempA = abs(A)
    FLOPS_ADD(1)
    FLOPS_ADD_N(1)
    FLOPS_ADD_N2(3)

    matrix_multiply_blas(n, n, n, A, A, &A2);
    FLOPS_ADD_N3(2)
    FLOPS_ADD(5)
    d6 = pow(_norm1est(A2, n, 3), d6_factor);
    d4 = pow(_norm1est(A2, n, 2), 0.25);
    eta1 = d4 > d6 ? d4 : d6;
    FLOPS_ADD(3)

    //matrix_multiply_blas(n, n, n, temp1, temp1, &temp3);  //2
    //matrix_multiply_blas(n, n, n, temp3, temp3, &tempA);  //4
    //matrix_multiply_blas(n, n, n, tempA, temp3, &temp3); //6
    //matrix_multiply_blas(n, n, n, temp3, temp1, &temp3);  //7
    FLOPS_ADD(1)
    if (eta1 <= 0.01495585217958292){
        alpha = _norm1est(tempA, n, 7);
        ell = normA * u * 100800.0;
        FLOPS_ADD(3)

        if (alpha <= ell) {
            mcm(A2, temp1, 1.0, n);
            m_sum_diag(temp1, 60.0, n);
            mcm(A2, temp2, 12.0, n);
            m_sum_diag(temp2, 120.0, n);
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(2)
        
            matrix_multiply_blas(n, n, n, A, temp1, &temp1);
            FLOPS_ADD_N3(2)
            FLOPS_ADD(5)

            m_sum(temp1, temp2, temp4, n);
            m_sub(temp2, temp1, temp2, n);
            FLOPS_ADD_N2(2)

            //DISCOVER THE NUMBER TO USE BLOCKING
            solve_system(temp4, temp2, X, n);

            free(temp1);
            free(temp2);
            free(temp4);
            free(tempA);
            free(A2);
            FLOPS_PRINT
            return 1;
        }
    }
    
    matrix_multiply_blas(n, n, n, A2, A2, &A4);
    FLOPS_ADD_N3(2)
    FLOPS_ADD(6)
    //matrix_multiply_blas(n, n, n, temp3, tempA, &temp3);  //11

    if (eta1 <= 0.2539398330063230){
        alpha = _norm1est(tempA, n, 11);
        ell = normA * u * 10059033600.0;
        FLOPS_ADD(3)

        if (alpha <= ell){
            //4 matrixes
            mcm(A2, temp1, 420.0, n);
            m_sum_diag(temp1, 15120.0, n);
            mcm(A2, temp2, 3360.0, n);
            m_sum_diag(temp2, 30240.0, n);
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(2)

            m_sum(A4, temp1, temp1, n);
            mcm_sum(A4, temp2, 30.0, n);
            FLOPS_ADD_N2(3)

            matrix_multiply_blas(n, n, n, A, temp1, &temp1);
            FLOPS_ADD_N3(2)
            FLOPS_ADD(5)

            //3 matrixes
            m_sum(temp1, temp2, temp4, n);
            m_sub(temp2, temp1, temp2, n);
            FLOPS_ADD_N2(2)

            solve_system(temp4, temp2, X, n);

            free(temp1);
            free(temp2);
            free(temp4);
            free(tempA);
            free(A2);
            free(A4);
            FLOPS_PRINT
            return 2;
        }
    }

    matrix_multiply_blas(n, n, n, A2, A4, &A6);
    FLOPS_ADD_N3(2)
    d8 = pow(_norm1est(A4, n, 2), 0.125);
    eta3 = d6 > d8 ? d6 : d8;
    FLOPS_ADD(8)

    //matrix_multiply_blas(n, n, n, temp3, tempA, &temp3);  //15

    if (eta3 <= 0.9504178996162932){
        alpha = _norm1est(tempA, n, 15);
        ell = normA * u * 4487938430976000.0;
        FLOPS_ADD(3)

        if (alpha <= ell){
            //5 matrixes
            mcm(A2, temp1, 277200.0, n);
            m_sum_diag(temp1, 8648640.0, n);
            mcm(A2, temp2, 1995840.0, n);
            m_sum_diag(temp2, 17297280.0, n);
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(2)

            mcm_sum(A4, temp1, 1512.0, n);
            mcm_sum(A4, temp2, 25200.0, n);
            FLOPS_ADD_N2(4)

            m_sum(A6, temp1, temp1, n);
            mcm_sum(A6, temp2, 56.0, n);
            FLOPS_ADD_N2(3)

            matrix_multiply_blas(n, n, n, A, temp1, &temp1);
            FLOPS_ADD_N3(2)
            FLOPS_ADD(5)

            //3 matrixes
            m_sum(temp1, temp2, temp4, n);
            m_sub(temp2, temp1, temp2, n);
            FLOPS_ADD_N2(2)

            solve_system(temp4, temp2, X, n);
            free(temp1);
            free(temp2);
            free(temp4);
            free(tempA);
            free(A2);
            free(A4);
            free(A6);
            free(A8);
            FLOPS_PRINT
            return 3;
        }
    }

    //matrix_multiply_blas(n, n, n, temp3, tempA, &temp3);  //19
    FLOPS_ADD(1)
    if (eta3 <= 2.097847961257068){
        alpha = _norm1est(tempA, n, 19);
        ell = normA * u * 5914384781877411840000.0;
        FLOPS_ADD(3)
        if (alpha <= ell){
            matrix_multiply_blas(n, n, n, A2, A6, &A8);
            FLOPS_ADD_N3(2)
            FLOPS_ADD(5)

            //6 matrixes
            mcm(A2, temp1, 302702400.0, n);
            m_sum_diag(temp1, 8821612800.0, n);
            mcm(A2, temp2, 2075673600.0, n);
            m_sum_diag(temp2, 17643225600.0, n);
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(2)

            mcm_sum(A4, temp1, 2162160.0, n);
            mcm_sum(A4, temp2, 30270240.0, n);
            FLOPS_ADD_N2(4)

            mcm_sum(A6, temp1, 3960.0, n);
            mcm_sum(A6, temp2, 110880.0, n);
            FLOPS_ADD_N2(4)

            m_sum(A8, temp1, temp1, n);
            mcm_sum(A8, temp2, 90.0, n);
            FLOPS_ADD_N2(3)

            matrix_multiply_blas(n, n, n, A, temp1, &temp1);
            FLOPS_ADD_N3(2)
            FLOPS_ADD(5)

            //3 matrixes
            m_sum(temp1, temp2, temp4, n);
            m_sub(temp2, temp1, temp2, n);
            FLOPS_ADD_N2(2)

            solve_system(temp4, temp2, X, n);
            free(temp1);
            free(temp2);
            free(temp4);
            free(tempA);
            free(A2);
            free(A4);
            free(A6);
            free(A8);
            FLOPS_PRINT
            return 3;
        }
    }
   
   // Up to line 27
   // Continues below
    matrix_multiply_blas(n, n, n, A4, A6, &temp1);
    FLOPS_ADD_N3(2)
    FLOPS_ADD(5)
    d10 = pow(_norm1est(temp1, n, 1), 0.1);
    eta4 = d8 > d10 ? d8 : d10;
    eta5 = eta3 > eta4 ? eta4 : eta3;
    FLOPS_ADD(3)

    s = 0;
    s_temp = eta5 * 0.2352941176;
    FLOPS_ADD(2)
    if (s_temp > 1.0){
        FLOPS_ADD(2)
        s = ceil(log2(s_temp));
    } 

    //mcm(tempA, tempA, pow(2, -s*4), n);
    //mcm(temp3, temp3, pow(2, -s*19), n);
    eta1 = pow(2, -s);
    FLOPS_ADD(2)
    normA = normA * eta1;
    mcm(tempA, tempA, eta1, n);
    FLOPS_ADD_N2(1)

    //matrix_multiply_blas(n, n, n, tempA, tempA, &tempA);  //8
    //matrix_multiply_blas(n, n, n, tempA, temp3, &temp3);  //27

    alpha = ceil(log2(_norm1est(tempA, n, 27.0)/(u * 113250775606021113483283660800000000.0 * normA)) * 0.03846153846);
    ell = alpha > 0 ? alpha : 0;
    s += ell;
    FLOPS_ADD(8)

    //line 31 - end

    //7-8 matrixes -- THIS CAN BE A LONG BLOCKING OPTIMIZATION (MAYBE ENDING UP ADDING SOME COMPUTATION)
    mcm(A, tempA, pow(2, -s), n);
    mcm(A2, A2, pow(4, -s), n);
    mcm(A4, A4, pow(16, -s), n);
    mcm(A6, A6, pow(64, -s), n);
    FLOPS_ADD_N2(4)
    FLOPS_ADD(4)

    //P_M & Q_M
    mcm(A2, temp1, 40840800.0, n);
    mcm_sum(A4, temp1, 16380.0, n);
    FLOPS_ADD_N2(3)

    mcm(A2, temp2, 1323241920.0, n);
    mcm_sum(A4, temp2, 960960.0, n);
    FLOPS_ADD_N2(3)

    m_sum(A6, temp1, temp1, n);
    mcm_sum(A6, temp2, 182.0, n);
    FLOPS_ADD_N2(3)

    matrix_multiply_blas(n, n, n, temp1, A6, &temp1);
    matrix_multiply_blas(n, n, n, temp2, A6, &temp2);
    FLOPS_ADD_N3(4)
    FLOPS_ADD(10)

    //Second part
    mcm_sum(A2, temp1, 1187353796428800.0, n);
    mcm_sum(A2, temp2, 7771770303897600.0, n);
    FLOPS_ADD_N2(4)

    mcm_sum(A4, temp1, 10559470521600.0, n);
    mcm_sum(A4, temp2, 129060195264000.0, n);
    FLOPS_ADD_N2(4)

    mcm_sum(A6, temp1, 33522128640.0, n);
    mcm_sum(A6, temp2, 670442572800.0, n);
    FLOPS_ADD_N2(4)

    m_sum_diag(temp1, 32382376266240000.0, n);
    m_sum_diag(temp2, 64764752532480000.0, n);
    FLOPS_ADD_N(2)

    //PROBABBLY 4 MATRIXES HERE

    matrix_multiply_blas(n, n, n, tempA, temp1, &temp1);
    FLOPS_ADD_N3(2)
    FLOPS_ADD(5)

    //3 matrixes

    m_sum(temp1, temp2, temp4, n);
    m_sub(temp2, temp1, temp2, n);
    FLOPS_ADD_N2(2)

    solve_system(temp4, temp2, X, n);

    if (s==0){ 
        free(temp1);
        free(temp2);
        free(temp4);
        free(tempA);
        free(A2);
        free(A4);
        free(A6);
        FLOPS_PRINT
        return 4;
    }

    if (is_triangular(tempA, n)) {
        double lambda1, lambda2, t12, temp5, temp6, temp7, temp8, lambda3, lambda4, lambda5, t13, t14, t15;
        double temp = 1.0;
        int index;
        for (int i = 0; i < n; i++)
            (*X)[i * n + i] = exp(tempA[i * n + i]);
        FLOPS_ADD_N(1)

        for (int _ = 1; _ <= s; _++){
            temp *= 2;
            FLOPS_ADD(1)

            matrix_multiply_blas(n, n, n, *X, *X, X);
            FLOPS_ADD_N3(2)
            FLOPS_ADD(5)

            lambda1 = tempA[0] * temp;
            FLOPS_ADD_N(15)
            int i;
            for (i = 0; i < n - 4; i+=4){
                index = n + 1;  // 9
                lambda2 = tempA[(i + 1) * index] * temp;
                lambda3 = tempA[(i + 2) * index] * temp;
                lambda4 = tempA[(i + 3) * index] * temp;
                lambda5 = tempA[(i + 4) * index] * temp;
                (*X)[i * index] = exp(lambda1);
                (*X)[(i + 1) * index] = exp(lambda2);
                (*X)[(i + 2) * index] = exp(lambda3);
                (*X)[(i + 3) * index] = exp(lambda4);
                t12 = tempA[i*index + 1] * temp;
                t13 = tempA[(i + 1) * index + 1] * temp;
                t14 = tempA[(i + 2) * index + 1] * temp;
                t15 = tempA[(i + 3) * index + 1] * temp;
                temp5 = lambda2 - lambda1;
                temp6 = lambda3 - lambda2;
                temp7 = lambda4 - lambda3;
                temp8 = lambda5 - lambda4;

                (*X)[i*index + 1] = t12 * exp((lambda1 + lambda2) * 0.5) * ((exp(temp5 * 0.5) - exp(temp5 * -0.5)) / temp5);
                (*X)[(i + 1) * index + 1] = t13 * exp((lambda2 + lambda3) * 0.5) * ((exp(temp6 * 0.5) - exp(temp6 * -0.5)) / temp6);
                (*X)[(i + 2) * index + 1] = t14 * exp((lambda3 + lambda4) * 0.5) * ((exp(temp7 * 0.5) - exp(temp7 * -0.5)) / temp7);
                (*X)[(i + 3) * index + 1] = t15 * exp((lambda4 + lambda5) * 0.5) * ((exp(temp8 * 0.5) - exp(temp8 * -0.5)) / temp8);
                lambda1 = lambda5;
            }
            for (; i < n - 1; i++){
                index = n + 1; //9
                (*X)[i * index] = exp(lambda1);
                lambda2 = tempA[(i + 1) * index] * temp;
                t12 = tempA[i*index + 1] * temp;
                temp5 = lambda2 - lambda1;

                (*X)[i*index + 1] = t12 * exp((lambda1 + lambda2) * 0.5) * ((exp(temp5 * 0.5) - exp(temp5 * -0.5)) / temp5);
                lambda1 = lambda2;
            }
            (*X)[n*n - 1] = exp(lambda1);

        }
        free(temp1);
        free(temp2);
        free(temp4);
        free(A2);
        free(A4);
        free(A6);
        free(tempA);
        FLOPS_PRINT
        return 5;
    } else {
        for (int i = 0; i < s; i++) {
            matrix_multiply_blas(n, n, n, *X, *X, X);
            FLOPS_ADD_N3(2)
            FLOPS_ADD(5)
        }
        free(temp1);
        free(temp2);
        free(temp4);
        free(A2);
        free(A4);
        free(A6);
        free(tempA);
        FLOPS_PRINT
        return 6;
    }
}
