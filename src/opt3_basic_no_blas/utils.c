#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include "mkl.h"

#ifdef FLOPS
double flops_opt_ = 0;
double flops_opt_N_ = 0;
double flops_opt_N2_ = 0;
double flops_opt_N3_ = 0;
#define FLOPS_RESET_ flops_opt_ = 0; flops_opt_N_ = 0; flops_opt_N2_ = 0; flops_opt_N3_ = 0;
#define FLOPS_ADD_(x) flops_opt_ += x;
#define FLOPS_ADD_N_(x) flops_opt_N_ += x;
#define FLOPS_ADD_N2_(x) flops_opt_N2_ += x;
#define FLOPS_ADD_N3_(x) flops_opt_N3_ += x;
#define FLOPS_PRINT_ printf("UTILS \nFlops: %f\n", flops_opt_); printf("Flops N: %f\n", flops_opt_N_); printf("Flops N^2: %f\n", flops_opt_N2_); printf("Flops N^3: %f\n", flops_opt_N3_);
#else
#define FLOPS_RESET_
#define FLOPS_ADD_(x)
#define FLOPS_ADD_N_(x)
#define FLOPS_ADD_N2_(x)
#define FLOPS_ADD_N3_(x)
#define FLOPS_PRINT_
#endif

// calculates the norm-1 of a matrix
double norm_1(double *A, int n) {
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

// calculates the infinite-norm of a matrix whihch is the maximum absolute row sum
double norm_inf(double *A, int n){
    double result = 0.0;
    double column_sum, column_sum2, column_sum3, column_sum4, column_sum5, column_sum6, column_sum7, column_sum8;

    if (n < 8){
        for (int i = 0; i < n; i++) {
            column_sum = 0.0;
            for (int j = 0; j < n; j++) {
                column_sum += fabs(A[i * n + j]);
            }
            result = fmax(column_sum, result);
        }
    }else{
        for (int i = 0; i < n; i++) {
            column_sum = 0.0;
            column_sum2 = 0.0;
            column_sum3 = 0.0;
            column_sum4 = 0.0;
            column_sum5 = 0.0;
            column_sum6 = 0.0;
            column_sum7 = 0.0;
            column_sum8 = 0.0;
            int j;
            for (j = 0; j < n - 7; j+=8) {
                column_sum += fabs(A[i*n + j]);
                column_sum2 += fabs(A[i*n + (j+1)]);
                column_sum3 += fabs(A[i*n + (j+2)]);
                column_sum4 += fabs(A[i*n + (j+3)]);
                column_sum5 += fabs(A[i*n + (j+4)]);
                column_sum6 += fabs(A[i*n + (j+5)]);
                column_sum7 += fabs(A[i*n + (j+6)]);
                column_sum8 += fabs(A[i*n + (j+7)]);
            }
            for (; j < n; j++) {
                column_sum += fabs(A[i*n + j]);
            }
            column_sum = column_sum + column_sum2 + column_sum3 + column_sum4 + column_sum5 + column_sum6 + column_sum7 + column_sum8;
            result = fmax(column_sum, result);
        }
    }
    return result;
}

// calculates the absolute value of a matrix whihch is the same matrix but all elements are positive
void m_abs(double *A, double *result, int n){
    for(int i = 0; i < n*n; i++)
        result[i] = fabs(A[i]);
}

// calculates the product between a matrix and a scalar (matrix-constant-multiplcation)
void mcm(double *A, double *result, double scalar, int n){
    for(int i = 0; i < n*n; i++)
        result[i] = A[i] * scalar;
}

void mcm_sum(double *A, double *result, double scalar, int n){
    for(int i = 0; i < n*n; i++)
        result[i] += A[i] * scalar;
}

void m_sum_diag(double *A, double scalar, int n){
    for(int i = 0; i < n; i++)
        A[i*n + i] += scalar;
}

void m_sub(double *A, double *B, double *result, int n){
    for(int i = 0; i < n*n; i++)
        result[i] = A[i] - B[i];
}

// calculates the sum between matrixes
void m_sum(double *A, double *B, double *result, int n){
    for(int i = 0; i < n*n; i++)
        result[i] = A[i] + B[i];
}

// calculates the transpose of a matrix
void m_transpose(double *A, double *result, int n){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++)
            result[i*n + j] = A[j*n + i];
    }
}

int is_triangular(double *A, int n) {
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++)
            if (A[i * n + j] != 0.0)
                return 0;
    }
    return 1;
}

// Function to allocate memory for a matrix
double *allocate_matrix(int n, int m) {
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

double *identity_matrix(int n){
    double *In = allocate_matrix(n, n);
    for (int i = 0; i < n*n; i++)
        In[i] = 0.0;
    for (int i = 0; i < n; i++){
        In[i*n + i] = 1.0;
    }
    return In;
}

bool invert_matrix(const double *input, double *output, int n) {
    double *A = allocate_matrix(n, n);
    FLOPS_ADD_(5)
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
    FLOPS_ADD_N3_(4)
    FLOPS_ADD_N2_(3)
    FLOPS_ADD_N_(3)
    FLOPS_ADD_N2_(1)

    free(A);
    return true;
}

bool invert_matrix2(const double *input, double *output, int n) {
    double *A = (double *)malloc(n * n * sizeof(double));
    FLOPS_ADD_(5)
    memcpy(A, input, n * n * sizeof(double));

    // Create row pointers for A matrix
    double **row_ptrs = (double **)malloc(n * sizeof(double *));
    FLOPS_ADD_(5)
    for (int row = 0; row < n; row++) {
        row_ptrs[row] = &A[row * n];
    }

    // Initialize output as an identity matrix
    for (int row = 0; row < n; row++) {
        int row_idx = row * n;
        for (int col = 0; col < n; col++) {
            output[row_idx + col] = (row == col) ? 1.0 : 0.0;
        }
    }
    FLOPS_ADD_N2_(1)

    for (int pivot = 0; pivot < n; pivot++) {
        int maxrow = pivot;

        for (int row = pivot + 1; row < n; row++) {
            if (fabs(row_ptrs[row][pivot]) > fabs(row_ptrs[maxrow][pivot])) {
                maxrow = row;
            }
        }

        if (fabs(row_ptrs[maxrow][pivot]) < 1e-8) {
            free(A);
            free(row_ptrs);
            return false;
        }

        // Swap rows by swapping pointers in row_ptrs
        if (maxrow != pivot) {
            double *temp_ptr = row_ptrs[pivot];
            row_ptrs[pivot] = row_ptrs[maxrow];
            row_ptrs[maxrow] = temp_ptr;
        }

        double pivot_value = 1.0 / row_ptrs[pivot][pivot];

        for (int col = 0; col < n; col++) {
            row_ptrs[pivot][col] *= pivot_value;
            output[pivot * n + col] *= pivot_value;
        }

        for (int row = 0; row < n; row++) {
            if (row != pivot) {
                double factor = row_ptrs[row][pivot];

                for (int col = 0; col < n; col++) {
                    row_ptrs[row][col] -= factor * row_ptrs[pivot][col];
                    output[row * n + col] -= factor * output[pivot * n + col];
                }
            }
        }
    }
    FLOPS_ADD_N3_(4)
    FLOPS_ADD_N2_(5.5)
    FLOPS_ADD_N_(3)

    free(A);
    free(row_ptrs);
    return true;
}

// Matrix multiplication with BLAS
// C = A * B
// A: m x k
// B: k x n
void matrix_multiply_blas(int m, int n, int k, double *A, double *B, double **X) {
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

void solve_system(double *P, double *Q, double **X, int n) {
    FLOPS_RESET_
    double *Q_inv = allocate_matrix(n, n);
    FLOPS_ADD_(6)
    if (!invert_matrix(Q, Q_inv, n)) {
        printf("Error: matrix is singular and cannot be inverted.\n");
        return;
    }
    // matrix_multiply_blas(n, n, n, Q_inv, P, X);
    // FLOPS_ADD_N3_(2)
    // FLOPS_ADD_(5)
    *X = Q_inv;
    FLOPS_PRINT_
}
