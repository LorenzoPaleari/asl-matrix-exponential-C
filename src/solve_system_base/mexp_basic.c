#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
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
#define FLOPS_PRINT printf("BASE\n"); printf("MAIN \nFlops: %f\n", flops_opt); printf("Flops N: %f\n", flops_opt_N); printf("Flops N^2: %f\n", flops_opt_N2); printf("Flops N^3: %f\n", flops_opt_N3);
#else
#define FLOPS_RESET
#define FLOPS_ADD(x)
#define FLOPS_ADD_N(x)
#define FLOPS_ADD_N2(x)
#define FLOPS_ADD_N3(x)
#define FLOPS_PRINT
#endif

// Function to allocate memory for a matrix
double *allocate_matrix(int n, int m) {
    double *A = malloc(n * m * sizeof(double));
    if (A == NULL) {
        printf("Error: memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }
    return A;
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
    FLOPS_ADD_N2(1)

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
    FLOPS_ADD_N3(4)
    FLOPS_ADD_N2(5.5)
    FLOPS_ADD_N(2)

    // Copy the inverted matrix to the output
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            output[row * n + col] = augmented[row * (2 * n) + col + n];
        }
    }

    free(augmented);
    return true;
}

void solve_system(double *P, double *Q, double **X, int n) {
    double *Q_inv = allocate_matrix(n, n);
    FLOPS_ADD(1)
    if (!invert_matrix(Q, Q_inv, n)) {
        printf("Error: matrix is singular and cannot be inverted.\n");
        return;
    }
    *X = Q_inv;
}

int mexp_no_blas(double *A, double *A2, double *A4, double *A6, double *A8, double *A10, double **X, int n){
    FLOPS_RESET
    solve_system(A2, A, X, n);
    FLOPS_PRINT
}
