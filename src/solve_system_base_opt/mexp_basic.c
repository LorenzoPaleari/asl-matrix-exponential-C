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
    FLOPS_ADD(6)
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
