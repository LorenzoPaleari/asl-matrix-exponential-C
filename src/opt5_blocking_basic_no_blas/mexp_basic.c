#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

// Function to allocate memory for a matrix
__attribute__((always_inline)) double *allocate_matrix(int n, int m) {
    // Align the memory to a 32-byte boundary to make the most out of AVX2 instructions
    size_t alignment = 32;

    // Calculate the size of the memory to allocate
    size_t size = n * m * sizeof(double);

    // Ensure the size is a multiple of the alignment
    if (size % alignment != 0) {
        size += alignment - (size % alignment);
    }

    // Allocate the memory
    double *A = (double *)aligned_alloc(alignment, size);
    if (A == NULL) {
        printf("Error: memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }
    return A;
}

__attribute__((always_inline)) void matrix_multiply_blas(int m, int n, int k, double *A, double *B, double **X) {
    double alpha = 1.0;
    double beta = 0.0;
    double *C =  allocate_matrix(m, n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, alpha, A, k, B, n, beta, C, n);
    
    if(*X != NULL){
        free(*X);
    }
    *X = C;
}

//NORMEST
__attribute__((always_inline)) double _norm1est(double *A, int n, int m) {
    int max_iter = 5;
    srand(24);
    double est_old = 0;
    int *ind_hist = (int *)malloc((n + 1) * sizeof(int));
    int ind_hist_len = 0;
    double *Y = (double *)malloc(n * 2 * sizeof(double));

    // Initialize Y
    double div_n = 1.0/n;
    FLOPS_ADD(2);
    double row_sum, row_sum2, row_sum3, row_sum4, row_sum5, row_sum6, row_sum7, row_sum8;
    if (n < 8)
        for (int i = 0; i < n; i++) {
            row_sum = 0.0;
            for (int j = 0; j < n; j++) {
                row_sum += fabs(A[i * n + j]);
            }
            long long n = (long long)(row_sum * div_n);
            Y[i * 2] = (double)n;
        }
    else
        for (int i = 0; i < n; i++) {
            row_sum = 0.0;
            row_sum2 = 0.0;
            row_sum3 = 0.0;
            row_sum4 = 0.0;
            row_sum5 = 0.0;
            row_sum6 = 0.0;
            row_sum7 = 0.0;
            row_sum8 = 0.0;
            int j;
            for (j = 0; j < n - 7; j+=8) {
                row_sum += fabs(A[i*n + j]);
                row_sum2 += fabs(A[i*n + j + 1]);
                row_sum3 += fabs(A[i*n + j + 2]);
                row_sum4 += fabs(A[i*n + j + 3]);
                row_sum5 += fabs(A[i*n + j + 4]);
                row_sum6 += fabs(A[i*n + j + 5]);
                row_sum7 += fabs(A[i*n + j + 6]);
                row_sum8 += fabs(A[i*n + j + 7]);
            }
            for (; j < n; j++) {
                row_sum += fabs(A[i*n + j]);
            }
            row_sum += row_sum2 + row_sum3 + row_sum4 + row_sum5 + row_sum6 + row_sum7 + row_sum8;
        long long n = (long long)(row_sum * div_n);
        Y[i * 2] = (double)n;
    }
    FLOPS_ADD_N2(1)
    FLOPS_ADD_N(1)

    // Random columns
    for (int i = 1; i < 2; i++) {
        int col_idx = rand() % n;
        for (int j = 0; j < n; j++) {
            Y[2 * j + i] = A[j * n + col_idx];
        }
        ind_hist[ind_hist_len++] = col_idx;
    }
    FLOPS_ADD_N(2)


    for (int k = 0; k < max_iter; k++) {
        
        for (int _ = 0; _ < m - 1; _++) {
                FLOPS_ADD_N2(4);
                FLOPS_ADD(5);
                matrix_multiply_blas(n, 2, n, A, Y, &Y);
        }
        
        
        double y_sums_0 = 0;
        double y_sums_1 = 0;
        double y_sums_2 = 0;
        double y_sums_3 = 0;
        double y_sums_4 = 0;
        double y_sums_5 = 0;
        double y_sums_6 = 0;
        double y_sums_7 = 0;
        FLOPS_ADD_N(2);
        if (n < 8)
            for (int j = 0; j < 2*n; j+=2) {
                y_sums_0 += Y[j] > 0 ? Y[j] : -Y[j];
                y_sums_1 += Y[j + 1] > 0 ? Y[j + 1] : -Y[j + 1];
            }
        else{
            int j;
            for (j = 0; j < 2*n - 7; j+=8) {
                y_sums_0 += Y[j] > 0 ? Y[j] : -Y[j];
                y_sums_1 += Y[j + 1] > 0 ? Y[j + 1] : -Y[j + 1];
                y_sums_2 += Y[j + 2] > 0 ? Y[j + 2] : -Y[j + 2];
                y_sums_3 += Y[j + 3] > 0 ? Y[j + 3] : -Y[j + 3];
                y_sums_4 += Y[j + 4] > 0 ? Y[j + 4] : -Y[j + 4];
                y_sums_5 += Y[j + 5] > 0 ? Y[j + 5] : -Y[j + 5];
                y_sums_6 += Y[j + 6] > 0 ? Y[j + 6] : -Y[j + 6];
                y_sums_7 += Y[j + 7] > 0 ? Y[j + 7] : -Y[j + 7];
            }
            for (; j < 2*n; j+=2) {
                y_sums_0 += Y[j] > 0 ? Y[j] : -Y[j];
                y_sums_1 += Y[j + 1] > 0 ? Y[j + 1] : -Y[j + 1];
            }
            y_sums_0 = y_sums_0 + y_sums_2 + y_sums_4 + y_sums_6;
            y_sums_1 = y_sums_1 + y_sums_3 + y_sums_5 + y_sums_7;
        }
        FLOPS_ADD_N(4);
        
        double est; 
        FLOPS_ADD(1);
        if (y_sums_0 < y_sums_1)
        {
            est = y_sums_1;
        }else{
            est = y_sums_0;
        }
        

        FLOPS_ADD(1);
        if (est <= est_old) {
            est = est_old;
            break;
        }
        est_old = est;

    bool all_in_hist = true;
    for (int i = 0; i < 2; i++) {
        bool in_hist = false;
        for (int j = 0; j < ind_hist_len; j++) {
            FLOPS_ADD(1);
            if (ind_hist[j] == i) {
                in_hist = true;
                break;
            }
        }
        FLOPS_ADD(1);
        if (!in_hist) {
            all_in_hist = false;
            break;
        }
    }
    FLOPS_ADD(1);
    if (all_in_hist) {
        break;
    } else {
        int pick;
        for (int i = 0; i < 2; i++) {
            int is_in_history;
            int limit = 0;
            do {
                is_in_history = 0;
                FLOPS_ADD_N(1);
                pick = rand() % n;
                FLOPS_ADD(1);
                limit++;
                for (int j = 0; j < ind_hist_len; j++) {
                    FLOPS_ADD(1);
                    if (ind_hist[j] == pick) {
                        is_in_history = 1;
                        break;
                    }
                }
            } while (is_in_history && limit < max_iter);
            FLOPS_ADD(2);
            ind_hist[ind_hist_len++] = pick;
            for (int j = 0; j < n; j++) {
                Y[2 * j + i] = A[j * n + pick];
            }
        }
    }
}
free(ind_hist);
free(Y);
return est_old;
}

// calculates the norm-1 of a matrix
__attribute__((always_inline)) double norm_1(double *A, int n) {
    double result = 0.0;
    double column_sum, column_sum2, column_sum3, column_sum4, column_sum5, column_sum6, column_sum7, column_sum8;

    if (n < 8) {
        for (int j = 0; j < n; j++) {
            column_sum = 0.0;
            for (int i = 0; i < n; i++) {
                column_sum += fabs(A[i * n + j]);
            }
            result = fmax(column_sum, result);
        }
    }else{
        for (int j = 0; j < n; j++) {
            column_sum = 0.0;
            column_sum2 = 0.0;
            column_sum3 = 0.0;
            column_sum4 = 0.0;
            column_sum5 = 0.0;
            column_sum6 = 0.0;
            column_sum7 = 0.0;
            column_sum8 = 0.0;
            int i;
            for (i = 0; i < n - 7; i+=8) {
                column_sum += fabs(A[i * n + j]);
                column_sum2 += fabs(A[(i+1) * n + j]);
                column_sum3 += fabs(A[(i+2) * n + j]);
                column_sum4 += fabs(A[(i+3) * n + j]);
                column_sum5 += fabs(A[(i+4) * n + j]);
                column_sum6 += fabs(A[(i+5) * n + j]);
                column_sum7 += fabs(A[(i+6) * n + j]);
                column_sum8 += fabs(A[(i+7) * n + j]);
            }
            for (; i < n; i++) {
                column_sum += fabs(A[i * n + j]);
            }
            column_sum = column_sum + column_sum2 + column_sum3 + column_sum4 + column_sum5 + column_sum6 + column_sum7 + column_sum8;
            result = fmax(column_sum, result);
        }

    }
    return result;
}

// calculates the absolute value of a matrix whihch is the same matrix but all elements are positive
__attribute__((always_inline)) void m_abs(double *A, double *result, int n){
    for(int i = 0; i < n*n; i++)
        result[i] = fabs(A[i]);
}

// calculates the product between a matrix and a scalar (matrix-constant-multiplcation)
__attribute__((always_inline)) void mcm(double *A, double *result, double scalar, int n){
    for(int i = 0; i < n*n; i++)
        result[i] = A[i] * scalar;
}

__attribute__((always_inline)) void mcm_sum(double *A, double *result, double scalar, int n){
    for(int i = 0; i < n*n; i++)
        result[i] += A[i] * scalar;
}

__attribute__((always_inline)) void m_sum_diag(double *A, double scalar, int n){
    for(int i = 0; i < n; i++)
        A[i*n + i] += scalar;
}

__attribute__((always_inline)) void m_sub(double *A, double *B, double *result, int n){
    for(int i = 0; i < n*n; i++)
        result[i] = A[i] - B[i];
}

// calculates the sum between matrixes
__attribute__((always_inline)) void m_sum(double *A, double *B, double *result, int n){
    for(int i = 0; i < n*n; i++)
        result[i] = A[i] + B[i];
}

__attribute__((always_inline)) int is_triangular(double *A, int n) {
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++)
            if (A[i * n + j] != 0.0)
                return 0;
    }
    return 1;
}

__attribute__((always_inline)) bool invert_matrix(const double *input, double *output, int n) {
    double *A = allocate_matrix(n, n);
    FLOPS_ADD(5)
    memcpy(A, input, n * n * sizeof(double));

    // Initialize output as an identity matrix
    for (int row = 0; row < n; row++) {
        int row_idx = row * n;
        for (int col = 0; col < n; col++) {
            output[row_idx + col] =  0.0;
        }
        output[row_idx + row] = 1.0;
    }

    for (int pivot = 0; pivot < n; pivot++) {
        int maxrow = pivot;

        double max_value = fabs(A[pivot * n + pivot]);
        double value = 0.0;
        for (int row = pivot + 1; row < n; row++) {
            value = fabs(A[row * n + pivot]);
            if (value > max_value) {
                maxrow = row;
                max_value = value;
            }
        }

        // Inline swap_rows and avoid repeated calculations
        if (maxrow != pivot) {
            int pivot_idx = pivot * n;
            int maxrow_idx = maxrow * n;
            for (int i = 0; i < n; i++) {
                double temp_a = A[pivot_idx + i];
                A[pivot_idx + i] = A[maxrow_idx + i];
                A[maxrow_idx + i] = temp_a;

                double temp_output = output[pivot_idx + i];
                output[pivot_idx + i] = output[maxrow_idx + i];
                output[maxrow_idx + i] = temp_output;
            }
        }
        double pivot_value = 1.0/A[pivot * n + pivot];

        for (int col = 0; col < n; col++) {
            A[pivot * n + col] *= pivot_value;
            output[pivot * n + col] *= pivot_value;
        }

        for (int row = 0; row < n; row++) {
            if (row != pivot) {
                double factor = A[row * n + pivot];

                for (int col = 0; col < n; col++) {
                    A[row * n + col] -= factor * A[pivot * n + col];
                    output[row * n + col] -= factor * output[pivot * n + col];
                }
            }
        }
    }
    FLOPS_ADD_N3(4)
    FLOPS_ADD_N2(3)
    FLOPS_ADD_N(3)
    FLOPS_ADD_N2(1)

    free(A);
    return true;
}

__attribute__((always_inline)) void solve_system(double *P, double *Q, double **X, int n) {
    double *Q_inv = allocate_matrix(n, n);
    if (*X == NULL) {
        *X = allocate_matrix(n, n);
    }
    FLOPS_ADD(6)
    if (!invert_matrix(Q, Q_inv, n)) {
        printf("Error: matrix is singular and cannot be inverted.\n");
        return;
    }
    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    //             n, n, n, 1.0, Q_inv, n, P, n, 0.0, *X, n);
    // FLOPS_ADD_N3(2)
    // FLOPS_ADD(5)
    *X = Q_inv;
}

int mexp_no_blas(double *A, double *A2, double *A4, double *A6, double *A8, double *A10, double **X, int n){
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

    //ALL work space we need for the algorithm
    FLOPS_ADD(40)

    //ALLOCATION OF ALLIGNED MEMORY
    size_t alignment = 32;
    size_t size = n * n * sizeof(double);
    if (size % alignment != 0) {
        size += alignment - (size % alignment);
    }
    double *temp1 = (double *)aligned_alloc(alignment, size);
    double *temp2 = (double *)aligned_alloc(alignment, size);
    double *temp4 = (double *)aligned_alloc(alignment, size);
    double *tempA = (double *)aligned_alloc(alignment, size);

    size = n*sizeof(double);
    if (size % alignment != 0) {
        size += alignment - (size % alignment);
    }
    double *temp3 = (double *)aligned_alloc(alignment, size);
    double *temp5 = (double *)aligned_alloc(alignment, size);

    normA = norm_1(A, n);   //normA used in ell
    m_abs(A, tempA, n);   //tempA = abs(A)
    FLOPS_ADD(1)
    FLOPS_ADD_N(1)
    FLOPS_ADD_N2(3)

    // // Calculate the size of the memory to allocate
    // size = n * n * sizeof(double);
    // // Ensure the size is a multiple of the alignment
    // if (size % alignment != 0) {
    //     size += alignment - (size % alignment);
    // }
    // // Allocate the memory
    // A2 = (double *)aligned_alloc(alignment, size);
    // if (A2 == NULL) {
    //     printf("Error: memory allocation failed.\n");
    //     exit(EXIT_FAILURE);
    // }
    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    //             n, n, n, 1.0, A, n, A, n, 0.0, A2, n);
    // FLOPS_ADD_N3(2)
    // FLOPS_ADD(5)
    d6 = pow(_norm1est(A2, n, 3), d6_factor);
    d4 = pow(_norm1est(A2, n, 2), 0.25);
    eta1 = d4 > d6 ? d4 : d6;
    FLOPS_ADD(3)

    FLOPS_ADD(1)
    if (eta1 <= 0.01495585217958292){
        alpha = _norm1est(tempA, n, 7);
        ell = normA * u * 100800.0;
        FLOPS_ADD(3)

        if (alpha <= ell) {
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    temp1[i*n + j] = A2[i*n + j];
                    temp2[i*n + j] = A2[i*n + j] * 12.0;
                }
                temp1[i*n + i] += 60.0;
                temp2[i*n + i] += 120.0;
            }
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(2)
        
            // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            //     n, n, n, 1.0, A, n, temp1, n, 0.0, temp4, n);
            // FLOPS_ADD_N3(2)

            for(int i = 0; i < n*n; i++) {
                temp4[i] = temp2[i] + temp1[i];
                temp2[i] = temp2[i] - temp1[i];
            }
            FLOPS_ADD_N2(2)

            //DISCOVER THE NUMBER TO USE BLOCKING
            solve_system(temp4, temp2, X, n);

            free(temp1);
            free(temp2);
            free(temp4);
            free(tempA);
            FLOPS_PRINT
            return 1;
        }
    }
    
    // size = n * n * sizeof(double);
    // if (size % alignment != 0) {
    //     size += alignment - (size % alignment);
    // }
    // // Allocate the memory
    // A4 = (double *)aligned_alloc(alignment, size);
    // if (A4 == NULL) {
    //     printf("Error: memory allocation failed.\n");
    //     exit(EXIT_FAILURE);
    // }
    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    //             n, n, n, 1.0, A2, n, A2, n, 0.0, A4, n);
    // FLOPS_ADD_N3(2)
    FLOPS_ADD(1)
    //matrix_multiply_blas(n, n, n, temp3, tempA, &temp3);  //11

    if (eta1 <= 0.2539398330063230){
        alpha = _norm1est(tempA, n, 11);
        ell = normA * u * 10059033600.0;
        FLOPS_ADD(3)

        if (alpha <= ell){
            //4 matrixes
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(2)
            FLOPS_ADD_N2(3)

            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    temp1[i*n + j] = A2[i*n + j] * 420.0 + A4[i*n + j];
                    temp2[i*n + j] = A2[i*n + j] * 3360.0 + A4[i*n + j] * 30.0;
                }
                temp1[i*n + i] += 15120.0;
                temp2[i*n + i] += 30240.0;
            }

            //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            //    n, n, n, 1.0, A, n, temp1, n, 0.0, temp4, n);
            //FLOPS_ADD_N3(2)

            for(int i = 0; i < n*n; i++) {
                temp4[i] = temp2[i] + temp1[i];
                temp2[i] = temp2[i] - temp1[i];
            }
            FLOPS_ADD_N2(2)

            //DISCOVER THE NUMBER TO USE BLOCKING
            solve_system(temp4, temp2, X, n);

            free(temp1);
            free(temp2);
            free(temp3);
            free(temp4);
            free(temp5);
            free(tempA);
            FLOPS_PRINT;
            return 2;
        }
    }

    // size = n * n * sizeof(double);
    // if (size % alignment != 0) {
    //     size += alignment - (size % alignment);
    // }
    // // Allocate the memory
    // A6 = (double *)aligned_alloc(alignment, size);
    // if (A6 == NULL) {
    //     printf("Error: memory allocation failed.\n");
    //     exit(EXIT_FAILURE);
    // }
    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    //             n, n, n, 1.0, A4, n, A2, n, 0.0, A6, n);
    // FLOPS_ADD_N3(2)
    FLOPS_ADD(1)
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
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    temp1[i*n + j] = A2[i*n + j] * 277200.0 + A4[i*n + j] * 1512.0 + A6[i*n + j];
                    temp2[i*n + j] = A2[i*n + j] * 1995840.0 + A4[i*n + j] * 25200.0 + A6[i*n + j]*56.0;
                }
                temp1[i*n + i] += 8648640.0;
                temp2[i*n + i] += 17297280.0;
            }

            //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            //    n, n, n, 1.0, A, n, temp1, n, 0.0, temp4, n);
            //FLOPS_ADD_N3(2)

            for(int i = 0; i < n*n; i++) {
                temp4[i] = temp2[i] + temp1[i];
                temp2[i] = temp2[i] - temp1[i];
            }
            FLOPS_ADD_N2(2)

            //DISCOVER THE NUMBER TO USE BLOCKING
            solve_system(temp4, temp2, X, n);
            free(temp1);
            free(temp2);
            free(temp3);
            free(temp4);
            free(tempA);
            free(temp5);
            FLOPS_PRINT;
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
            // size = n * n * sizeof(double);
            // if (size % alignment != 0) {
            //     size += alignment - (size % alignment);
            // }
            // // Allocate the memory
            // A8 = (double *)aligned_alloc(alignment, size);
            // if (A8 == NULL) {
            //     printf("Error: memory allocation failed.\n");
            //     exit(EXIT_FAILURE);
            // }
            // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            //             n, n, n, 1.0, A4, n, A4, n, 0.0, A8, n);
            // FLOPS_ADD_N3(2)
            // FLOPS_ADD(5)

            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    temp1[i*n + j] = A2[i*n + j] * 302702400.0 + A4[i*n + j] * 2162160.0 + A6[i*n + j]*3960.0 + A8[i*n + j];
                    temp2[i*n + j] = A2[i*n + j] * 2075673600.0 + A4[i*n + j] * 30270240.0 + A6[i*n + j]*110880.0 + A8[i*n + j]*90.0;
                }
                temp1[i*n + i] += 8821612800.0;
                temp2[i*n + i] += 17643225600.0;
            }

            //6 matrixes
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(2)
            FLOPS_ADD_N2(4)
            FLOPS_ADD_N2(4)
            FLOPS_ADD_N2(3)

            //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            //    n, n, n, 1.0, A, n, temp1, n, 0.0, temp4, n);
            //FLOPS_ADD_N3(2)

            for(int i = 0; i < n*n; i++) {
                temp4[i] = temp2[i] + temp1[i];
                temp2[i] = temp2[i] - temp1[i];
            }
            FLOPS_ADD_N2(2)

            //DISCOVER THE NUMBER TO USE BLOCKING
            solve_system(temp4, temp2, X, n);
            free(temp1);
            free(temp2);
            free(temp3);
            free(temp4);
            free(tempA);
            free(temp5);
            FLOPS_PRINT;
            return 3;
        }
    }
   
   // Up to line 27
   // Continues below
    // matrix_multiply_blas(n, n, n, A4, A6, &temp1);
    // FLOPS_ADD_N3(2)
    // FLOPS_ADD(5)
    d10 = pow(_norm1est(A10, n, 1), 0.1);
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
    double value = pow(2, -s);
    double value2 = pow(4, -s);
    double value3 = pow(16, -s);
    double value4 = pow(64, -s);
    for (int i = 0; i < n*n; i++){
        tempA[i] = A[i] * value;
        A2[i] = A2[i] * value2;
        A4[i] = A4[i] * value3;
        A6[i] = A6[i] * value4;
    }
    FLOPS_ADD_N2(4)
    FLOPS_ADD(4)

    //P_M & Q_M
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            temp1[i*n + j] = A2[i*n + j] * 40840800.0 + A4[i*n + j] * 16380.0 + A6[i*n + j];
            temp2[i*n + j] = A2[i*n + j] * 1323241920.0 + A4[i*n + j] * 960960.0 + A6[i*n + j]*182.0;
        }
    }
    FLOPS_ADD_N2(3)
    FLOPS_ADD_N2(3)
    FLOPS_ADD_N2(3)

    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    //             n, n, n, 1.0, A6, n, temp1, n, 0.0, temp4, n);  //Temp1 --> Temp4
    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    //             n, n, n, 1.0, A6, n, temp2, n, 0.0, temp1, n);  //Temp2 --> Temp1
    // FLOPS_ADD_N3(4)

    //Second part
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            temp1[i*n + j] += A2[i*n + j] * 1187353796428800.0 + A4[i*n + j] * 10559470521600.0 + A6[i*n + j]*33522128640.0;
            temp2[i*n + j] += A2[i*n + j] * 7771770303897600.0 + A4[i*n + j] * 129060195264000.0 + A6[i*n + j]*670442572800.0;
        }
        temp1[i*n + i] += 32382376266240000.0;
        temp2[i*n + i] += 64764752532480000.0;
    }
    FLOPS_ADD_N2(4)
    FLOPS_ADD_N2(4)
    FLOPS_ADD_N2(4)
    FLOPS_ADD_N(2)

    //PROBABBLY 4 MATRIXES HERE

    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,      //Temp2 --> Temp1
    //             n, n, n, 1.0, tempA, n, temp4, n, 0.0, temp2, n);  //Temp1 --> Temp2
    // FLOPS_ADD_N3(2)

    for(int i = 0; i < n*n; i++) {
        temp4[i] = temp2[i] + temp1[i];
        temp2[i] = temp2[i] - temp1[i];
    }
    FLOPS_ADD_N2(2)

    //DISCOVER THE NUMBER TO USE BLOCKING
    solve_system(temp4, temp2, X, n);

    if (s==0){ 
        free(temp1);
        free(temp2);
        free(temp4);
        free(tempA);
        FLOPS_PRINT
        return 4;
    }

    if (is_triangular(tempA, n)) {
        double lambda1, lambda2, t12, temp5_2, temp6, temp7, temp8, lambda3, lambda4, lambda5, t13, t14, t15;
        double temp = 1.0;
        int index;
        for (int i = 0; i < n; i++)
            (*X)[i * n + i] = exp(tempA[i * n + i]);
        FLOPS_ADD_N(1)

        int _;
        for (_ = 1; _ < s; _+=2){
            temp *= 2;
            FLOPS_ADD(1)

            // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            //     n, n, n, 1.0, *X, n, *X, n, 0.0, temp1, n);
            // FLOPS_ADD_N3(2)

            lambda1 = tempA[0] * temp;
            FLOPS_ADD_N(15)
            int i;
            for (i = 0; i < n - 4; i+=4){
                index = n + 1;  // 9
                lambda2 = tempA[(i + 1) * index] * temp;
                lambda3 = tempA[(i + 2) * index] * temp;
                lambda4 = tempA[(i + 3) * index] * temp;
                lambda5 = tempA[(i + 4) * index] * temp;
                temp1[i * index] = exp(lambda1);
                temp1[(i + 1) * index] = exp(lambda2);
                temp1[(i + 2) * index] = exp(lambda3);
                temp1[(i + 3) * index] = exp(lambda4);
                t12 = tempA[i*index + 1] * temp;
                t13 = tempA[(i + 1) * index + 1] * temp;
                t14 = tempA[(i + 2) * index + 1] * temp;
                t15 = tempA[(i + 3) * index + 1] * temp;
                temp5_2 = lambda2 - lambda1;
                temp6 = lambda3 - lambda2;
                temp7 = lambda4 - lambda3;
                temp8 = lambda5 - lambda4;

                temp1[i*index + 1] = t12 * exp((lambda1 + lambda2) * 0.5) * ((exp(temp5_2 * 0.5) - exp(temp5_2 * -0.5)) / temp5_2);
                temp1[(i + 1) * index + 1] = t13 * exp((lambda2 + lambda3) * 0.5) * ((exp(temp6 * 0.5) - exp(temp6 * -0.5)) / temp6);
                temp1[(i + 2) * index + 1] = t14 * exp((lambda3 + lambda4) * 0.5) * ((exp(temp7 * 0.5) - exp(temp7 * -0.5)) / temp7);
                temp1[(i + 3) * index + 1] = t15 * exp((lambda4 + lambda5) * 0.5) * ((exp(temp8 * 0.5) - exp(temp8 * -0.5)) / temp8);
                lambda1 = lambda5;
            }
            for (; i < n - 1; i++){
                index = n + 1; //9
                temp1[i * index] = exp(lambda1);
                lambda2 = tempA[(i + 1) * index] * temp;
                t12 = tempA[i*index + 1] * temp;
                temp5_2 = lambda2 - lambda1;

                temp1[i*index + 1] = t12 * exp((lambda1 + lambda2) * 0.5) * ((exp(temp5_2 * 0.5) - exp(temp5_2 * -0.5)) / temp5_2);
                lambda1 = lambda2;
            }
            temp1[n*n - 1] = exp(lambda1);

            //SECOND ITERATION
            temp *= 2;
            FLOPS_ADD(1)

            // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            //     n, n, n, 1.0, temp1, n, temp1, n, 0.0, *X, n);
            // FLOPS_ADD_N3(2)

            lambda1 = tempA[0] * temp;
            FLOPS_ADD_N(15)
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
                temp5_2 = lambda2 - lambda1;
                temp6 = lambda3 - lambda2;
                temp7 = lambda4 - lambda3;
                temp8 = lambda5 - lambda4;

                (*X)[i*index + 1] = t12 * exp((lambda1 + lambda2) * 0.5) * ((exp(temp5_2 * 0.5) - exp(temp5_2 * -0.5)) / temp5_2);
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
                temp5_2 = lambda2 - lambda1;

                (*X)[i*index + 1] = t12 * exp((lambda1 + lambda2) * 0.5) * ((exp(temp5_2 * 0.5) - exp(temp5_2 * -0.5)) / temp5_2);
                lambda1 = lambda2;
            }
            (*X)[n*n - 1] = exp(lambda1);
        }
        //CLEANUP IF S IS ODD
        for (; _ <= s; _++){
            temp *= 2;
            FLOPS_ADD(1)

            // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            //     n, n, n, 1.0, *X, n, *X, n, 0.0, temp1, n);
            // FLOPS_ADD_N3(2)

            lambda1 = tempA[0] * temp;
            FLOPS_ADD_N(15)
            int i;
            for (i = 0; i < n - 4; i+=4){
                index = n + 1;  // 9
                lambda2 = tempA[(i + 1) * index] * temp;
                lambda3 = tempA[(i + 2) * index] * temp;
                lambda4 = tempA[(i + 3) * index] * temp;
                lambda5 = tempA[(i + 4) * index] * temp;
                temp1[i * index] = exp(lambda1);
                temp1[(i + 1) * index] = exp(lambda2);
                temp1[(i + 2) * index] = exp(lambda3);
                temp1[(i + 3) * index] = exp(lambda4);
                t12 = tempA[i*index + 1] * temp;
                t13 = tempA[(i + 1) * index + 1] * temp;
                t14 = tempA[(i + 2) * index + 1] * temp;
                t15 = tempA[(i + 3) * index + 1] * temp;
                temp5_2 = lambda2 - lambda1;
                temp6 = lambda3 - lambda2;
                temp7 = lambda4 - lambda3;
                temp8 = lambda5 - lambda4;

                temp1[i*index + 1] = t12 * exp((lambda1 + lambda2) * 0.5) * ((exp(temp5_2 * 0.5) - exp(temp5_2 * -0.5)) / temp5_2);
                temp1[(i + 1) * index + 1] = t13 * exp((lambda2 + lambda3) * 0.5) * ((exp(temp6 * 0.5) - exp(temp6 * -0.5)) / temp6);
                temp1[(i + 2) * index + 1] = t14 * exp((lambda3 + lambda4) * 0.5) * ((exp(temp7 * 0.5) - exp(temp7 * -0.5)) / temp7);
                temp1[(i + 3) * index + 1] = t15 * exp((lambda4 + lambda5) * 0.5) * ((exp(temp8 * 0.5) - exp(temp8 * -0.5)) / temp8);
                lambda1 = lambda5;
            }
            for (; i < n - 1; i++){
                index = n + 1; //9
                temp1[i * index] = exp(lambda1);
                lambda2 = tempA[(i + 1) * index] * temp;
                t12 = tempA[i*index + 1] * temp;
                temp5_2 = lambda2 - lambda1;

                temp1[i*index + 1] = t12 * exp((lambda1 + lambda2) * 0.5) * ((exp(temp5_2 * 0.5) - exp(temp5_2 * -0.5)) / temp5_2);
                lambda1 = lambda2;
            }
            temp1[n*n - 1] = exp(lambda1);

            *X = temp1;
        }

        if (s % 2 == 0)
            free(temp1);
        free(temp2);
        free(temp4);
        free(tempA);
        FLOPS_PRINT
        return 5;
    } else {
        free(temp1);
        free(temp2);
        free(temp4);
        free(tempA);
        FLOPS_PRINT
        return 6;
    }
}
