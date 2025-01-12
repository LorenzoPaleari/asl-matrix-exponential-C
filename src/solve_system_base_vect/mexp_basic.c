#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mkl.h"
#include <immintrin.h>

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
    double temp;
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

        double max_value = A[pivot*n + pivot] > 0 ? A[pivot*n + pivot] : -A[pivot*n + pivot];
        for (int j = pivot + 1; j < n; j++){
            temp = A[j*n + pivot] > 0 ? A[j*n + pivot] : -A[j*n + pivot];
            //temp = Q[pivot*n + j] > 0 ? Q[pivot*n + j] : -Q[pivot*n + j];
            if (temp > max_value){
                max_value = temp;
                maxrow = j;
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

        int col;
        __m256d pivot_vec = _mm256_set1_pd(pivot_value);
        for (col = 0; col < n - 7; col+=8) {
            __m256d A_vec1 = _mm256_load_pd(&A[pivot * n + col]);
            __m256d A_vec2 = _mm256_load_pd(&A[pivot * n + col + 4]);
            __m256d output_vec1 = _mm256_load_pd(&output[pivot * n + col]);
            __m256d output_vec2 = _mm256_load_pd(&output[pivot * n + col + 4]);

            _mm256_store_pd(&A[pivot * n + col], _mm256_mul_pd(A_vec1, pivot_vec));
            _mm256_store_pd(&A[pivot * n + col + 4], _mm256_mul_pd(A_vec2, pivot_vec));
            _mm256_store_pd(&output[pivot * n + col], _mm256_mul_pd(output_vec1, pivot_vec));
            _mm256_store_pd(&output[pivot * n + col + 4], _mm256_mul_pd(output_vec2, pivot_vec));
        }
        for (;col < n; col++){
            A[pivot * n + col] *= pivot_value;
            output[pivot * n + col] *= pivot_value;
        }

        int row;
        for (row = 0; row < pivot - 3; row+=4) {
            __m256d factor_vec = _mm256_set1_pd(A[row * n + pivot]);
            __m256d factor_vec2 = _mm256_set1_pd(A[(row + 1) * n + pivot]);
            __m256d factor_vec3 = _mm256_set1_pd(A[(row + 2) * n + pivot]);
            __m256d factor_vec4 = _mm256_set1_pd(A[(row + 3) * n + pivot]);

            for (int col = 0; col < n - 3; col+=4) {
                __m256d A_vec1 = _mm256_load_pd(&A[pivot * n + col]);
                __m256d output_vec1 = _mm256_load_pd(&output[pivot * n + col]);

                __m256d row_vec1 = _mm256_load_pd(&A[row * n + col]);
                __m256d row_vec2 = _mm256_load_pd(&output[row * n + col]);
                __m256d row_vec3 = _mm256_load_pd(&A[(row + 1) * n + col]);
                __m256d row_vec4 = _mm256_load_pd(&output[(row + 1) * n + col]);
                __m256d row_vec5 = _mm256_load_pd(&A[(row + 2) * n + col]);
                __m256d row_vec6 = _mm256_load_pd(&output[(row + 2) * n + col]);
                __m256d row_vec7 = _mm256_load_pd(&A[(row + 3) * n + col]);
                __m256d row_vec8 = _mm256_load_pd(&output[(row + 3) * n + col]);

                _mm256_store_pd(&A[row * n + col], _mm256_fnmadd_pd(factor_vec, A_vec1, row_vec1));
                _mm256_store_pd(&output[row * n + col], _mm256_fnmadd_pd(factor_vec, output_vec1, row_vec2));
                _mm256_store_pd(&A[(row + 1) * n + col], _mm256_fnmadd_pd(factor_vec2, A_vec1, row_vec3));
                _mm256_store_pd(&output[(row + 1) * n + col], _mm256_fnmadd_pd(factor_vec2, output_vec1, row_vec4));
                _mm256_store_pd(&A[(row + 2) * n + col], _mm256_fnmadd_pd(factor_vec3, A_vec1, row_vec5));
                _mm256_store_pd(&output[(row + 2) * n + col], _mm256_fnmadd_pd(factor_vec3, output_vec1, row_vec6));
                _mm256_store_pd(&A[(row + 3) * n + col], _mm256_fnmadd_pd(factor_vec4, A_vec1, row_vec7));
                _mm256_store_pd(&output[(row + 3) * n + col], _mm256_fnmadd_pd(factor_vec4, output_vec1, row_vec8));
            }
        }
        for (; row < pivot; row++) {
            __m256d factor_vec = _mm256_set1_pd(A[row * n + pivot]);

            for (int col = 0; col < n - 3; col+=4) {
                __m256d A_vec1 = _mm256_load_pd(&A[pivot * n + col]);
                __m256d output_vec1 = _mm256_load_pd(&output[pivot * n + col]);
                __m256d row_vec1 = _mm256_load_pd(&A[row * n + col]);
                __m256d row_vec2 = _mm256_load_pd(&output[row * n + col]);

                _mm256_store_pd(&A[row * n + col], _mm256_fnmadd_pd(factor_vec, A_vec1, row_vec1));
                _mm256_store_pd(&output[row * n + col], _mm256_fnmadd_pd(factor_vec, output_vec1, row_vec2));
            }
        }
        for (row = pivot + 1; row < n - 3; row+=4) {
            __m256d factor_vec = _mm256_set1_pd(A[row * n + pivot]);
            __m256d factor_vec2 = _mm256_set1_pd(A[(row + 1) * n + pivot]);
            __m256d factor_vec3 = _mm256_set1_pd(A[(row + 2) * n + pivot]);
            __m256d factor_vec4 = _mm256_set1_pd(A[(row + 3) * n + pivot]);

            for (int col = 0; col < n - 3; col+=4) {
                __m256d A_vec1 = _mm256_loadu_pd(&A[pivot * n + col]);
                __m256d output_vec1 = _mm256_loadu_pd(&output[pivot * n + col]);

                __m256d row_vec1 = _mm256_loadu_pd(&A[row * n + col]);
                __m256d row_vec2 = _mm256_loadu_pd(&output[row * n + col]);
                __m256d row_vec3 = _mm256_loadu_pd(&A[(row + 1) * n + col]);
                __m256d row_vec4 = _mm256_loadu_pd(&output[(row + 1) * n + col]);
                __m256d row_vec5 = _mm256_loadu_pd(&A[(row + 2) * n + col]);
                __m256d row_vec6 = _mm256_loadu_pd(&output[(row + 2) * n + col]);
                __m256d row_vec7 = _mm256_loadu_pd(&A[(row + 3) * n + col]);
                __m256d row_vec8 = _mm256_loadu_pd(&output[(row + 3) * n + col]);

                _mm256_storeu_pd(&A[row * n + col], _mm256_fnmadd_pd(factor_vec, A_vec1, row_vec1));
                _mm256_storeu_pd(&output[row * n + col], _mm256_fnmadd_pd(factor_vec, output_vec1, row_vec2));
                _mm256_storeu_pd(&A[(row + 1) * n + col], _mm256_fnmadd_pd(factor_vec2, A_vec1, row_vec3));
                _mm256_storeu_pd(&output[(row + 1) * n + col], _mm256_fnmadd_pd(factor_vec2, output_vec1, row_vec4));
                _mm256_storeu_pd(&A[(row + 2) * n + col], _mm256_fnmadd_pd(factor_vec3, A_vec1, row_vec5));
                _mm256_storeu_pd(&output[(row + 2) * n + col], _mm256_fnmadd_pd(factor_vec3, output_vec1, row_vec6));
                _mm256_storeu_pd(&A[(row + 3) * n + col], _mm256_fnmadd_pd(factor_vec4, A_vec1, row_vec7));
                _mm256_storeu_pd(&output[(row + 3) * n + col], _mm256_fnmadd_pd(factor_vec4, output_vec1, row_vec8));
            }
        }
        for (; row < n; row++) {
            __m256d factor_vec = _mm256_set1_pd(A[row * n + pivot]);

            for (int col = 0; col < n - 3; col+=4) {
                __m256d A_vec1 = _mm256_loadu_pd(&A[pivot * n + col]);
                __m256d output_vec1 = _mm256_loadu_pd(&output[pivot * n + col]);
                __m256d row_vec1 = _mm256_loadu_pd(&A[row * n + col]);
                __m256d row_vec2 = _mm256_loadu_pd(&output[row * n + col]);

                _mm256_storeu_pd(&A[row * n + col], _mm256_fnmadd_pd(factor_vec, A_vec1, row_vec1));
                _mm256_storeu_pd(&output[row * n + col], _mm256_fnmadd_pd(factor_vec, output_vec1, row_vec2));
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
