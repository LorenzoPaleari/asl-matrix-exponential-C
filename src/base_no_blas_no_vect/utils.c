#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
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

//calculate the 1-norm of a matrix whihch is the maximum absolute column sum
double norm_1(double *A, int n){
    double result = 0.0;
    double temp[n];

    for(int j = 0; j < n; j++){
        temp[j] = 0.0;
        for(int i = 0; i < n; i++)
            temp[j] += fabs(A[i*n + j]);
    }

    for(int i = 0; i < n; i++)
        result = fmax(temp[i], result);

    return result;
}

//calculate the absolute value of a matrix whihch is the same matrix but all elements are positive
void m_abs(double *A, double *result, int n){
    for(int i = 0; i < n*n; i++)
        result[i] = fabs(A[i]);
}

//calculate the product between a matrix and a scalar (matrix-constant-multiplcation)
void mcm(double *A, double *result, double scalar, int n){
    for(int i = 0; i < n*n; i++)
        result[i] = A[i] * scalar;
}

//calculate the sum between matrixes
void m_sum(double *A, double *B, double *result, int n){
    for(int i = 0; i < n*n; i++)
        result[i] = A[i] + B[i];
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
    double *A = malloc(n * m * sizeof(double));
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

void swap_rows(double *matrix, int row1, int row2, int n) {
    for (int i = 0; i < n; i++) {
        double temp = matrix[row1 * n + i];
        matrix[row1 * n + i] = matrix[row2 * n + i];
        matrix[row2 * n + i] = temp;
    }
}

bool invert_matrix(const double *input, double *output, int n) {
    double *augmented = (double *)malloc(n * 2 * n * sizeof(double));

    // Create an augmented matrix [A | I]
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            augmented[row * (2 * n) + col] = input[row * n + col];
            augmented[row * (2 * n) + col + n] = (row == col) ? 1.0 : 0.0;
        }
    }
    FLOPS_ADD_N2_(1)

    // Gaussian elimination with partial pivoting
    for (int pivot = 0; pivot < n; pivot++) {
        int maxrow = pivot;

        // Find the row with the maximum value in the current column
        for (int row = pivot + 1; row < n; row++) {
            if (fabs(augmented[row * (2 * n) + pivot]) > fabs(augmented[maxrow * (2 * n) + pivot])) {
                maxrow = row;
            }
        }

        // If the maximum value is 0, the matrix is singular and non-invertible
        if (fabs(augmented[maxrow * (2 * n) + pivot]) < 1e-8) {
            free(augmented);
            return false;
        }

        // Swap the rows to move the pivot to the diagonal
        swap_rows(augmented, pivot, maxrow, 2 * n);

        // Normalize the pivot row
        double pivot_value = augmented[pivot * (2 * n) + pivot];
        for (int col = 0; col < 2 * n; col++) {
            augmented[pivot * (2 * n) + col] /= pivot_value;
        }

        // Eliminate the pivot column in all other rows
        for (int row = 0; row < n; row++) {
            if (row != pivot) {
                double factor = augmented[row * (2 * n) + pivot];
                for (int col = 0; col < 2 * n; col++) {
                    augmented[row * (2 * n) + col] -= factor * augmented[pivot * (2 * n) + col];
                }
            }
        }
    }
    FLOPS_ADD_N3_(4)
    FLOPS_ADD_N2_(5.5)
    FLOPS_ADD_N_(2)

    // Copy the inverted matrix to the output
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            output[row * n + col] = augmented[row * (2 * n) + col + n];
        }
    }

    free(augmented);
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
    FLOPS_ADD_(1)
    if (!invert_matrix(Q, Q_inv, n)) {
        printf("Error: matrix is singular and cannot be inverted.\n");
        return;
    }
    // matrix_multiply_blas_(n, n, n, Q_inv, P, X);
    // FLOPS_ADD_N3_(2)
    *X = Q_inv;
    FLOPS_PRINT_
}