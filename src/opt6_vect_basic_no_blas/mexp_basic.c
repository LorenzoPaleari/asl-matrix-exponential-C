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
    __m256d abs_mask = _mm256_set1_pd(-0.0);

    // Initialize Y
    double div_n = 1.0/n;
    FLOPS_ADD(2);
    double row_sum;
    for (int i = 0; i < n; i++) {
        row_sum = 0.0;
        int j;
        __m256d y_sum = _mm256_setzero_pd();
        __m256d y_sum_1 = _mm256_setzero_pd();
        __m256d y_sum_2 = _mm256_setzero_pd();
        __m256d y_sum_3 = _mm256_setzero_pd();
        for (j = 0; j < n - 15; j+=16) {
            __m256d y = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[i * n + j]));
            __m256d y_1 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[i * n + j + 4]));
            __m256d y_2 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[i * n + j + 8]));
            __m256d y_3 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[i * n + j + 12]));

            y_sum = _mm256_add_pd(y_sum, y);
            y_sum_1 = _mm256_add_pd(y_sum_1, y_1);
            y_sum_2 = _mm256_add_pd(y_sum_2, y_2);
            y_sum_3 = _mm256_add_pd(y_sum_3, y_3);
        }
        y_sum = _mm256_add_pd(y_sum, y_sum_1);
        y_sum_2 = _mm256_add_pd(y_sum_2, y_sum_3);
        y_sum = _mm256_add_pd(y_sum, y_sum_2);

        double y_sum_arr[4];
        _mm256_store_pd(y_sum_arr, y_sum);

        row_sum += y_sum_arr[0] + y_sum_arr[1] + y_sum_arr[2] + y_sum_arr[3];

        for (; j < n; j++) {
            row_sum += A[i * n + j] > 0 ? A[i * n + j] : -A[i * n + j];
        }
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
        int j;
        __m256d y_sum = _mm256_setzero_pd();
        __m256d y_sum_1 = _mm256_setzero_pd();
        __m256d y_sum_2 = _mm256_setzero_pd();
        __m256d y_sum_3 = _mm256_setzero_pd();

        for (j = 0; j < 2*n - 15; j+=16) {
            __m256d var0 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(Y + j));
            __m256d var1 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(Y + j + 4));
            __m256d var2 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(Y + j + 8));
            __m256d var3 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(Y + j + 12));

            y_sum = _mm256_add_pd(y_sum, var0);
            y_sum_1 = _mm256_add_pd(y_sum_1, var1);
            y_sum_2 = _mm256_add_pd(y_sum_2, var2);
            y_sum_3 = _mm256_add_pd(y_sum_3, var3);
        }

        y_sum = _mm256_add_pd(y_sum, y_sum_1);
        y_sum_2 = _mm256_add_pd(y_sum_2, y_sum_3);

        y_sum = _mm256_add_pd(y_sum, y_sum_2);

        y_sums_0 = y_sum[0] + y_sum[2];
        y_sums_1 = y_sum[1] + y_sum[3];

        for (; j < 2*n; j+=8) {
            y_sums_0 += Y[j] > 0 ? Y[j] : -Y[j];
            y_sums_1 += Y[j + 1] > 0 ? Y[j + 1] : -Y[j + 1];
            y_sums_2 += Y[j + 2] > 0 ? Y[j + 2] : -Y[j + 2];
            y_sums_3 += Y[j + 3] > 0 ? Y[j + 3] : -Y[j + 3];
            y_sums_4 += Y[j + 4] > 0 ? Y[j + 4] : -Y[j + 4];
            y_sums_5 += Y[j + 5] > 0 ? Y[j + 5] : -Y[j + 5];
            y_sums_6 += Y[j + 6] > 0 ? Y[j + 6] : -Y[j + 6];
            y_sums_7 += Y[j + 7] > 0 ? Y[j + 7] : -Y[j + 7];
        }
        y_sums_0 = y_sums_0 + y_sums_2 + y_sums_4 + y_sums_6;
        y_sums_1 = y_sums_1 + y_sums_3 + y_sums_5 + y_sums_7;

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
    __m256d abs_mask = _mm256_set1_pd(-0.0);
    double result = 0.0;
    double column_sum = 0.0;
    int i;
    for(i = 0; i < n - 7; i+=8){
        __m256d y_sum = _mm256_setzero_pd();
        __m256d y_sum_1 = _mm256_setzero_pd();
        for (int j = 0; j < n - 3; j+=4){
            __m256d A_vec = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[j*n + i]));
            __m256d A_vec2 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[(j + 1)*n + i]));
            __m256d A_vec3 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[(j + 2)*n + i]));
            __m256d A_vec4 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[(j + 3)*n + i]));
            __m256d A_vec5 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[j*n + i + 4]));
            __m256d A_vec6 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[(j + 1)*n + i + 4]));
            __m256d A_vec7 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[(j + 2)*n + i + 4]));
            __m256d A_vec8 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[(j + 3)*n + i + 4]));

            __m256d sum = _mm256_add_pd(A_vec, A_vec2);
            __m256d sum2 = _mm256_add_pd(A_vec3, A_vec4);
            __m256d sum3 = _mm256_add_pd(A_vec5, A_vec6);
            __m256d sum4 = _mm256_add_pd(A_vec7, A_vec8);

            sum = _mm256_add_pd(sum2, sum);
            sum3 = _mm256_add_pd(sum3, sum4);

            y_sum = _mm256_add_pd(y_sum, sum);
            y_sum_1 = _mm256_add_pd(y_sum_1, sum3);
        }
        __m256d max = _mm256_max_pd(y_sum, y_sum_1);

        double max_[4];
        _mm256_store_pd(max_, max);

        max_[0] = max_[0] > max_[1] ? max_[0] : max_[1];
        max_[2] = max_[2] > max_[3] ? max_[2] : max_[3];
        result = max_[0] > max_[2] ? max_[0] : max_[2];
    }
    for(; i < n; i++){
        for (int j = 0; j < n; j ++){
                column_sum += A[j*n + i] > 0 ? A[j*n + i] : -A[j*n + i];
        }
        result = column_sum > result ? column_sum : result;
    }
    return result;
}

// calculates the absolute value of a matrix whihch is the same matrix but all elements are positive
__attribute__((always_inline)) void m_abs(double *A, double *result, int n){
    __m256d abs_mask = _mm256_set1_pd(-0.0);
    for(int i = 0; i < n*n - 15; i+=16){
        _mm256_store_pd(&result[i], _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[i])));
        _mm256_store_pd(&result[i + 4], _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[i + 4])));
        _mm256_store_pd(&result[i + 8], _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[i + 8])));
        _mm256_store_pd(&result[i + 12], _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A[i + 12])));
    }
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
            __m256d scalar = _mm256_set1_pd(12.0);
            for(int i = 0; i < n - 3; i+=4) {
                for(int j = 0; j < n - 3; j+=4) {
                    __m256d y = _mm256_load_pd(&A2[i * n + j]);
                    __m256d y_1 = _mm256_load_pd(&A2[(i + 1) * n + j]);
                    __m256d y_2 = _mm256_load_pd(&A2[(i + 2) * n + j]);
                    __m256d y_3 = _mm256_load_pd(&A2[(i + 3) * n + j]);

                    __m256d y_5 = _mm256_mul_pd(y, scalar);
                    __m256d y_6 = _mm256_mul_pd(y_1, scalar);
                    __m256d y_7 = _mm256_mul_pd(y_2, scalar);
                    __m256d y_8 = _mm256_mul_pd(y_3, scalar);

                    _mm256_store_pd(&temp2[i * n + j], y_5);
                    _mm256_store_pd(&temp2[(i + 1) * n + j], y_6);
                    _mm256_store_pd(&temp2[(i + 2) * n + j], y_7);
                    _mm256_store_pd(&temp2[(i + 3) * n + j], y_8);
                    _mm256_store_pd(&temp1[i * n + j], y);
                    _mm256_store_pd(&temp1[(i + 1) * n + j], y_1);
                    _mm256_store_pd(&temp1[(i + 2) * n + j], y_2);
                    _mm256_store_pd(&temp1[(i + 3) * n + j], y_3);
                }
                temp1[i*n + i] += 60.0;
                temp1[(i + 1)*n + i + 1] += 60.0;
                temp1[(i + 2)*n + i + 2] += 60.0;
                temp1[(i + 3)*n + i + 3] += 60.0;
                temp2[i*n + i] += 120.0;
                temp2[(i + 1)*n + i + 1] += 120.0;
                temp2[(i + 2)*n + i + 2] += 120.0;
                temp2[(i + 3)*n + i + 3] += 120.0;
            }
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(2)
    

            for(int i = 0; i < n*n - 15; i+=16) {
                __m256d y = _mm256_load_pd(&temp4[i]);
                __m256d y_1 = _mm256_load_pd(&temp4[i + 4]);
                __m256d y_2 = _mm256_load_pd(&temp4[i + 8]);
                __m256d y_3 = _mm256_load_pd(&temp4[i + 12]);

                __m256d y_5 = _mm256_load_pd(&temp2[i]);
                __m256d y_6 = _mm256_load_pd(&temp2[i + 4]);
                __m256d y_7 = _mm256_load_pd(&temp2[i + 8]);
                __m256d y_8 = _mm256_load_pd(&temp2[i + 12]);

                __m256d y_9 = _mm256_add_pd(y, y_5);
                __m256d y_10 = _mm256_add_pd(y_1, y_6);
                __m256d y_11 = _mm256_add_pd(y_2, y_7);
                __m256d y_12 = _mm256_add_pd(y_3, y_8);

                __m256d y_13 = _mm256_sub_pd(y_5, y);
                __m256d y_14 = _mm256_sub_pd(y_6, y_1);
                __m256d y_15 = _mm256_sub_pd(y_7, y_2);
                __m256d y_16 = _mm256_sub_pd(y_8, y_3);

                _mm256_store_pd(&temp1[i], y_9);
                _mm256_store_pd(&temp1[i + 4], y_10);
                _mm256_store_pd(&temp1[i + 8], y_11);
                _mm256_store_pd(&temp1[i + 12], y_12);
                _mm256_store_pd(&temp2[i], y_13);
                _mm256_store_pd(&temp2[i + 4], y_14);
                _mm256_store_pd(&temp2[i + 8], y_15);
                _mm256_store_pd(&temp2[i + 12], y_16);
            }
            FLOPS_ADD_N2(2)

            //DISCOVER THE NUMBER TO USE BLOCKING
            solve_system(temp1, temp2, X, n);

            free(temp1);
            free(temp2);
            free(temp4);
            free(tempA);
            FLOPS_PRINT
            return 1;
        }
    }

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

            __m256d scalar = _mm256_set1_pd(420.0);
            __m256d scalar2 = _mm256_set1_pd(3360.0);
            __m256d scalar3 = _mm256_set1_pd(30.0);
            for(int i = 0; i < n; i+=4) {
                for(int j = 0; j < n; j+=4) {
                    __m256d y = _mm256_load_pd(&A2[i*n + j]);
                    __m256d y_1 = _mm256_load_pd(&A2[(i + 1)*n + j]);
                    __m256d y_2 = _mm256_load_pd(&A2[(i + 2)*n + j]);
                    __m256d y_3 = _mm256_load_pd(&A2[(i + 3)*n + j]);

                    __m256d y_4 = _mm256_load_pd(&A4[i*n + j]);
                    __m256d y_5 = _mm256_load_pd(&A4[(i + 1)*n + j]);
                    __m256d y_6 = _mm256_load_pd(&A4[(i + 2)*n + j]);
                    __m256d y_7 = _mm256_load_pd(&A4[(i + 3)*n + j]);

                    __m256d y_8 = _mm256_fmadd_pd(y, scalar, y_4);
                    __m256d y_9 = _mm256_fmadd_pd(y_1, scalar, y_5);
                    __m256d y_10 = _mm256_fmadd_pd(y_2, scalar, y_6);
                    __m256d y_11 = _mm256_fmadd_pd(y_3, scalar, y_7);

                    y_4 = _mm256_mul_pd(y_4, scalar3);
                    y_5 = _mm256_mul_pd(y_5, scalar3);
                    y_6 = _mm256_mul_pd(y_6, scalar3);
                    y_7 = _mm256_mul_pd(y_7, scalar3);

                    __m256d y_12 = _mm256_fmadd_pd(y, scalar2, y_4);
                    __m256d y_13 = _mm256_fmadd_pd(y_1, scalar2, y_5);
                    __m256d y_14 = _mm256_fmadd_pd(y_2, scalar2, y_6);
                    __m256d y_15 = _mm256_fmadd_pd(y_3, scalar2, y_7);

                    _mm256_store_pd(&temp1[i*n + j], y_8);
                    _mm256_store_pd(&temp1[(i + 1)*n + j], y_9);
                    _mm256_store_pd(&temp1[(i + 2)*n + j], y_10);
                    _mm256_store_pd(&temp1[(i + 3)*n + j], y_11);
                    _mm256_store_pd(&temp2[i*n + j], y_12);
                    _mm256_store_pd(&temp2[(i + 1)*n + j], y_13);
                    _mm256_store_pd(&temp2[(i + 2)*n + j], y_14);
                    _mm256_store_pd(&temp2[(i + 3)*n + j], y_15);
                }
                temp1[i*n + i] += 15120.0;
                temp1[(i + 1)*n + i + 1] += 15120.0;
                temp1[(i + 2)*n + i + 2] += 15120.0;
                temp1[(i + 3)*n + i + 3] += 15120.0;
                temp2[i*n + i] += 30240.0;
                temp2[(i + 1)*n + i + 1] += 30240.0;
                temp2[(i + 2)*n + i + 2] += 30240.0;
                temp2[(i + 3)*n + i + 3] += 30240.0;
            }

            for(int i = 0; i < n*n - 15; i+=16) {
                __m256d y = _mm256_load_pd(&temp4[i]);
                __m256d y_1 = _mm256_load_pd(&temp4[i + 4]);
                __m256d y_2 = _mm256_load_pd(&temp4[i + 8]);
                __m256d y_3 = _mm256_load_pd(&temp4[i + 12]);

                __m256d y_5 = _mm256_load_pd(&temp2[i]);
                __m256d y_6 = _mm256_load_pd(&temp2[i + 4]);
                __m256d y_7 = _mm256_load_pd(&temp2[i + 8]);
                __m256d y_8 = _mm256_load_pd(&temp2[i + 12]);

                __m256d y_9 = _mm256_add_pd(y, y_5);
                __m256d y_10 = _mm256_add_pd(y_1, y_6);
                __m256d y_11 = _mm256_add_pd(y_2, y_7);
                __m256d y_12 = _mm256_add_pd(y_3, y_8);

                __m256d y_13 = _mm256_sub_pd(y_5, y);
                __m256d y_14 = _mm256_sub_pd(y_6, y_1);
                __m256d y_15 = _mm256_sub_pd(y_7, y_2);
                __m256d y_16 = _mm256_sub_pd(y_8, y_3);

                _mm256_store_pd(&temp1[i], y_9);
                _mm256_store_pd(&temp1[i + 4], y_10);
                _mm256_store_pd(&temp1[i + 8], y_11);
                _mm256_store_pd(&temp1[i + 12], y_12);
                _mm256_store_pd(&temp2[i], y_13);
                _mm256_store_pd(&temp2[i + 4], y_14);
                _mm256_store_pd(&temp2[i + 8], y_15);
                _mm256_store_pd(&temp2[i + 12], y_16);
            }
            FLOPS_ADD_N2(2)

            //DISCOVER THE NUMBER TO USE BLOCKING
            solve_system(temp1, temp2, X, n);

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
            __m256d scalar = _mm256_set1_pd(277200.0);
            __m256d scalar2 = _mm256_set1_pd(1512.0);
            __m256d scalar3 = _mm256_set1_pd(56.0);
            __m256d scalar4 = _mm256_set1_pd(25200.0);
            __m256d scalar5 = _mm256_set1_pd(1995840.0);
            for(int i = 0; i < n; i+=4) {
                for(int j = 0; j < n; j+=4) {
                    __m256d y = _mm256_load_pd(&A2[i*n + j]);
                    __m256d y_1 = _mm256_load_pd(&A2[(i + 1)*n + j]);
                    __m256d y_2 = _mm256_load_pd(&A2[(i + 2)*n + j]);
                    __m256d y_3 = _mm256_load_pd(&A2[(i + 3)*n + j]);

                    __m256d y_4 = _mm256_load_pd(&A4[i*n + j]);
                    __m256d y_5 = _mm256_load_pd(&A4[(i + 1)*n + j]);
                    __m256d y_6 = _mm256_load_pd(&A4[(i + 2)*n + j]);
                    __m256d y_7 = _mm256_load_pd(&A4[(i + 3)*n + j]);

                    __m256d y_8 = _mm256_load_pd(&A6[i*n + j]);
                    __m256d y_9 = _mm256_load_pd(&A6[(i + 1)*n + j]);
                    __m256d y_10 = _mm256_load_pd(&A6[(i + 2)*n + j]);
                    __m256d y_11 = _mm256_load_pd(&A6[(i + 3)*n + j]);

                    __m256d y_12 = _mm256_fmadd_pd(y, scalar, y_8);
                    y_12 = _mm256_fmadd_pd(y_4, scalar2, y_12);
                    __m256d y_13 = _mm256_fmadd_pd(y_1, scalar, y_9);
                    y_13 = _mm256_fmadd_pd(y_5, scalar2, y_13);
                    __m256d y_14 = _mm256_fmadd_pd(y_2, scalar, y_10);
                    y_14 = _mm256_fmadd_pd(y_6, scalar2, y_14);
                    __m256d y_15 = _mm256_fmadd_pd(y_3, scalar, y_11);
                    y_15 = _mm256_fmadd_pd(y_7, scalar2, y_15);

                    y_8 = _mm256_mul_pd(y_8, scalar3);
                    y_9 = _mm256_mul_pd(y_9, scalar3);
                    y_10 = _mm256_mul_pd(y_10, scalar3);
                    y_11 = _mm256_mul_pd(y_11, scalar3);

                    __m256d y_16 = _mm256_fmadd_pd(y, scalar5, y_8);
                    y_16 = _mm256_fmadd_pd(y_4, scalar4, y_16);
                    __m256d y_17 = _mm256_fmadd_pd(y_1, scalar5, y_9);
                    y_17 = _mm256_fmadd_pd(y_5, scalar4, y_17);
                    __m256d y_18 = _mm256_fmadd_pd(y_2, scalar5, y_10);
                    y_18 = _mm256_fmadd_pd(y_6, scalar4, y_18);
                    __m256d y_19 = _mm256_fmadd_pd(y_3, scalar5, y_11);
                    y_19 = _mm256_fmadd_pd(y_7, scalar4, y_19);
                    
                    _mm256_store_pd(&temp1[i*n + j], y_12);
                    _mm256_store_pd(&temp1[(i + 1)*n + j], y_13);
                    _mm256_store_pd(&temp1[(i + 2)*n + j], y_14);
                    _mm256_store_pd(&temp1[(i + 3)*n + j], y_15);
                    _mm256_store_pd(&temp2[i*n + j], y_16);
                    _mm256_store_pd(&temp2[(i + 1)*n + j], y_17);
                    _mm256_store_pd(&temp2[(i + 2)*n + j], y_18);
                    _mm256_store_pd(&temp2[(i + 3)*n + j], y_19);
                }
                temp1[i*n + i] += 8648640.0;
                temp1[(i + 1) * n + i + 1] += 8648640.0;
                temp1[(i + 2) * n + i + 2] += 8648640.0;
                temp1[(i + 3) * n + i + 3] += 8648640.0;
                temp2[i*n + i] += 17297280.0;
                temp2[(i + 1) * n + i + 1] += 17297280.0;
                temp2[(i + 2) * n + i + 2] += 17297280.0;
                temp2[(i + 3) * n + i + 3] += 17297280.0;
            }

            for(int i = 0; i < n*n - 15; i+=16) {
                __m256d y = _mm256_load_pd(&temp4[i]);
                __m256d y_1 = _mm256_load_pd(&temp4[i + 4]);
                __m256d y_2 = _mm256_load_pd(&temp4[i + 8]);
                __m256d y_3 = _mm256_load_pd(&temp4[i + 12]);

                __m256d y_5 = _mm256_load_pd(&temp2[i]);
                __m256d y_6 = _mm256_load_pd(&temp2[i + 4]);
                __m256d y_7 = _mm256_load_pd(&temp2[i + 8]);
                __m256d y_8 = _mm256_load_pd(&temp2[i + 12]);

                __m256d y_9 = _mm256_add_pd(y, y_5);
                __m256d y_10 = _mm256_add_pd(y_1, y_6);
                __m256d y_11 = _mm256_add_pd(y_2, y_7);
                __m256d y_12 = _mm256_add_pd(y_3, y_8);

                __m256d y_13 = _mm256_sub_pd(y_5, y);
                __m256d y_14 = _mm256_sub_pd(y_6, y_1);
                __m256d y_15 = _mm256_sub_pd(y_7, y_2);
                __m256d y_16 = _mm256_sub_pd(y_8, y_3);

                _mm256_store_pd(&temp1[i], y_9);
                _mm256_store_pd(&temp1[i + 4], y_10);
                _mm256_store_pd(&temp1[i + 8], y_11);
                _mm256_store_pd(&temp1[i + 12], y_12);
                _mm256_store_pd(&temp2[i], y_13);
                _mm256_store_pd(&temp2[i + 4], y_14);
                _mm256_store_pd(&temp2[i + 8], y_15);
                _mm256_store_pd(&temp2[i + 12], y_16);
            }
            FLOPS_ADD_N2(2)

            //DISCOVER THE NUMBER TO USE BLOCKING
            solve_system(temp1, temp2, X, n);
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
            __m256d scalar = _mm256_set1_pd(302702400.0);
            __m256d scalar_1 = _mm256_set1_pd(2162160.0);
            __m256d scalar_2 = _mm256_set1_pd(3960.0);
            __m256d scalar_3 = _mm256_set1_pd(90.0);
            __m256d scalar_4 = _mm256_set1_pd(2075673600.0);
            __m256d scalar_5 = _mm256_set1_pd(30270240.0);
            __m256d scalar_6 = _mm256_set1_pd(110880.0);
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    __m256d y = _mm256_load_pd(&A2[i*n + j]);
                    __m256d y_1 = _mm256_load_pd(&A2[(i + 1)*n + j]);
                    __m256d y_2 = _mm256_load_pd(&A2[(i + 2)*n + j]);
                    __m256d y_3 = _mm256_load_pd(&A2[(i + 3)*n + j]);

                    __m256d y_4 = _mm256_load_pd(&A4[i*n + j]);
                    __m256d y_5 = _mm256_load_pd(&A4[(i + 1)*n + j]);
                    __m256d y_6 = _mm256_load_pd(&A4[(i + 2)*n + j]);
                    __m256d y_7 = _mm256_load_pd(&A4[(i + 3)*n + j]);

                    __m256d y_8 = _mm256_load_pd(&A6[i*n + j]);
                    __m256d y_9 = _mm256_load_pd(&A6[(i + 1)*n + j]);
                    __m256d y_10 = _mm256_load_pd(&A6[(i + 2)*n + j]);
                    __m256d y_11 = _mm256_load_pd(&A6[(i + 3)*n + j]);

                    __m256d y_12 = _mm256_load_pd(&A8[i*n + j]);
                    __m256d y_13 = _mm256_load_pd(&A8[(i + 1)*n + j]);
                    __m256d y_14 = _mm256_load_pd(&A8[(i + 2)*n + j]);
                    __m256d y_15 = _mm256_load_pd(&A8[(i + 3)*n + j]);

                    __m256d y_16 = _mm256_fmadd_pd(y, scalar, y_12);
                    y_16 = _mm256_fmadd_pd(y_4, scalar_1, y_16);
                    y_16 = _mm256_fmadd_pd(y_8, scalar_2, y_16);
                    __m256d y_17 = _mm256_fmadd_pd(y_1, scalar, y_13);
                    y_17 = _mm256_fmadd_pd(y_5, scalar_1, y_17);
                    y_17 = _mm256_fmadd_pd(y_9, scalar_2, y_17);
                    __m256d y_18 = _mm256_fmadd_pd(y_2, scalar, y_14);
                    y_18 = _mm256_fmadd_pd(y_6, scalar_1, y_18);
                    y_18 = _mm256_fmadd_pd(y_10, scalar_2, y_18);
                    __m256d y_19 = _mm256_fmadd_pd(y_3, scalar, y_15);
                    y_19 = _mm256_fmadd_pd(y_7, scalar_1, y_19);
                    y_19 = _mm256_fmadd_pd(y_11, scalar_2, y_19);

                    y_12 = _mm256_mul_pd(y_12, scalar_3);
                    y_13 = _mm256_mul_pd(y_13, scalar_3);
                    y_14 = _mm256_mul_pd(y_14, scalar_3);
                    y_15 = _mm256_mul_pd(y_15, scalar_3);

                    __m256d y_20 = _mm256_fmadd_pd(y, scalar_4, y_12);
                    y_20 = _mm256_fmadd_pd(y_4, scalar_5, y_20);
                    y_20 = _mm256_fmadd_pd(y_8, scalar_6, y_20);
                    __m256d y_21 = _mm256_fmadd_pd(y_1, scalar_4, y_13);
                    y_21 = _mm256_fmadd_pd(y_5, scalar_5, y_21);
                    y_21 = _mm256_fmadd_pd(y_9, scalar_6, y_21);
                    __m256d y_22 = _mm256_fmadd_pd(y_2, scalar_4, y_14);
                    y_22 = _mm256_fmadd_pd(y_6, scalar_5, y_22);
                    y_22 = _mm256_fmadd_pd(y_10, scalar_6, y_22);
                    __m256d y_23 = _mm256_fmadd_pd(y_3, scalar_4, y_15);
                    y_23 = _mm256_fmadd_pd(y_7, scalar_5, y_23);
                    y_23 = _mm256_fmadd_pd(y_11, scalar_6, y_23);
                    
                    _mm256_store_pd(&temp1[i*n + j], y_16);
                    _mm256_store_pd(&temp1[(i + 1)*n + j], y_17);
                    _mm256_store_pd(&temp1[(i + 2)*n + j], y_18);
                    _mm256_store_pd(&temp1[(i + 3)*n + j], y_19);
                    _mm256_store_pd(&temp2[i*n + j], y_20);
                    _mm256_store_pd(&temp2[(i + 1)*n + j], y_21);
                    _mm256_store_pd(&temp2[(i + 2)*n + j], y_22);
                    _mm256_store_pd(&temp2[(i + 3)*n + j], y_23);
                }
                temp1[i*n + i] += 8821612800.0;
                temp1[(i + 1)*n + i] += 8821612800.0;
                temp1[(i + 2)*n + i] += 8821612800.0;
                temp1[(i + 3)*n + i] += 8821612800.0;
                temp2[i*n + i] += 17643225600.0;
                temp2[(i + 1)*n + i] += 17643225600.0;
                temp2[(i + 2)*n + i] += 17643225600.0;
                temp2[(i + 3)*n + i] += 17643225600.0;
            }

            //6 matrixes
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(2)
            FLOPS_ADD_N2(4)
            FLOPS_ADD_N2(4)
            FLOPS_ADD_N2(3)

            for(int i = 0; i < n*n - 15; i+=16) {
                __m256d y = _mm256_load_pd(&temp4[i]);
                __m256d y_1 = _mm256_load_pd(&temp4[i + 4]);
                __m256d y_2 = _mm256_load_pd(&temp4[i + 8]);
                __m256d y_3 = _mm256_load_pd(&temp4[i + 12]);

                __m256d y_5 = _mm256_load_pd(&temp2[i]);
                __m256d y_6 = _mm256_load_pd(&temp2[i + 4]);
                __m256d y_7 = _mm256_load_pd(&temp2[i + 8]);
                __m256d y_8 = _mm256_load_pd(&temp2[i + 12]);

                __m256d y_9 = _mm256_add_pd(y, y_5);
                __m256d y_10 = _mm256_add_pd(y_1, y_6);
                __m256d y_11 = _mm256_add_pd(y_2, y_7);
                __m256d y_12 = _mm256_add_pd(y_3, y_8);

                __m256d y_13 = _mm256_sub_pd(y_5, y);
                __m256d y_14 = _mm256_sub_pd(y_6, y_1);
                __m256d y_15 = _mm256_sub_pd(y_7, y_2);
                __m256d y_16 = _mm256_sub_pd(y_8, y_3);

                _mm256_store_pd(&temp1[i], y_9);
                _mm256_store_pd(&temp1[i + 4], y_10);
                _mm256_store_pd(&temp1[i + 8], y_11);
                _mm256_store_pd(&temp1[i + 12], y_12);
                _mm256_store_pd(&temp2[i], y_13);
                _mm256_store_pd(&temp2[i + 4], y_14);
                _mm256_store_pd(&temp2[i + 8], y_15);
                _mm256_store_pd(&temp2[i + 12], y_16);
            }
            FLOPS_ADD_N2(2)

            //DISCOVER THE NUMBER TO USE BLOCKING
            solve_system(temp1, temp2, X, n);
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

    __m256d scal = _mm256_set1_pd(eta1);
    for (int i = 0; i < n*n - 15; i+=16){
        __m256d y = _mm256_load_pd(&tempA[i]);
        __m256d y_1 = _mm256_load_pd(&tempA[i + 4]);
        __m256d y_2 = _mm256_load_pd(&tempA[i + 8]);
        __m256d y_3 = _mm256_load_pd(&tempA[i + 12]);

        _mm256_store_pd(&tempA[i], _mm256_mul_pd(y, scal));
        _mm256_store_pd(&tempA[i + 4], _mm256_mul_pd(y_1, scal));
        _mm256_store_pd(&tempA[i + 8], _mm256_mul_pd(y_2, scal));
        _mm256_store_pd(&tempA[i + 12], _mm256_mul_pd(y_3, scal));
    }
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
    __m256d val = _mm256_set1_pd(value);
    __m256d val1 = _mm256_set1_pd(value2);
    __m256d val2 = _mm256_set1_pd(value3);
    __m256d val3 = _mm256_set1_pd(value4);
    for (int i = 0; i < n*n - 3; i+=4){
        __m256d y = _mm256_load_pd(&A[i]);
        __m256d y_1 = _mm256_load_pd(&A2[i]);
        __m256d y_2 = _mm256_load_pd(&A4[i]);
        __m256d y_3 = _mm256_load_pd(&A6[i]);

        y = _mm256_mul_pd(y, val);
        y_1 = _mm256_mul_pd(y_1, val1);
        y_2 = _mm256_mul_pd(y_2, val2);
        y_3 = _mm256_mul_pd(y_3, val3);

        _mm256_store_pd(&tempA[i], y);
        _mm256_store_pd(&A2[i], y_1);
        _mm256_store_pd(&A4[i], y_2);
        _mm256_store_pd(&A6[i], y_3);
    }
    FLOPS_ADD_N2(4)
    FLOPS_ADD(4)

    //P_M & Q_M
    __m256d scalar = _mm256_set1_pd(40840800.0);
    __m256d scalar2 = _mm256_set1_pd(16380.0);
    __m256d scalar3 = _mm256_set1_pd(182.0);
    __m256d scalar4 = _mm256_set1_pd(960960.0);
    __m256d scalar5 = _mm256_set1_pd(1323241920.0);
    for(int i = 0; i < n; i+=4) {
        for(int j = 0; j < n; j+=4) {
            __m256d y = _mm256_load_pd(&A2[i*n + j]);
            __m256d y_1 = _mm256_load_pd(&A2[(i + 1)*n + j]);
            __m256d y_2 = _mm256_load_pd(&A2[(i + 2)*n + j]);
            __m256d y_3 = _mm256_load_pd(&A2[(i + 3)*n + j]);

            __m256d y_4 = _mm256_load_pd(&A4[i*n + j]);
            __m256d y_5 = _mm256_load_pd(&A4[(i + 1)*n + j]);
            __m256d y_6 = _mm256_load_pd(&A4[(i + 2)*n + j]);
            __m256d y_7 = _mm256_load_pd(&A4[(i + 3)*n + j]);

            __m256d y_8 = _mm256_load_pd(&A6[i*n + j]);
            __m256d y_9 = _mm256_load_pd(&A6[(i + 1)*n + j]);
            __m256d y_10 = _mm256_load_pd(&A6[(i + 2)*n + j]);
            __m256d y_11 = _mm256_load_pd(&A6[(i + 3)*n + j]);

            __m256d y_12 = _mm256_fmadd_pd(y, scalar, y_8);
            y_12 = _mm256_fmadd_pd(y_4, scalar2, y_12);
            __m256d y_13 = _mm256_fmadd_pd(y_1, scalar, y_9);
            y_13 = _mm256_fmadd_pd(y_5, scalar2, y_13);
            __m256d y_14 = _mm256_fmadd_pd(y_2, scalar, y_10);
            y_14 = _mm256_fmadd_pd(y_6, scalar2, y_14);
            __m256d y_15 = _mm256_fmadd_pd(y_3, scalar, y_11);
            y_15 = _mm256_fmadd_pd(y_7, scalar2, y_15);

            y_8 = _mm256_mul_pd(y_8, scalar3);
            y_9 = _mm256_mul_pd(y_9, scalar3);
            y_10 = _mm256_mul_pd(y_10, scalar3);
            y_11 = _mm256_mul_pd(y_11, scalar3);

            __m256d y_16 = _mm256_fmadd_pd(y, scalar5, y_8);
            y_16 = _mm256_fmadd_pd(y_4, scalar4, y_16);
            __m256d y_17 = _mm256_fmadd_pd(y_1, scalar5, y_9);
            y_17 = _mm256_fmadd_pd(y_5, scalar4, y_17);
            __m256d y_18 = _mm256_fmadd_pd(y_2, scalar5, y_10);
            y_18 = _mm256_fmadd_pd(y_6, scalar4, y_18);
            __m256d y_19 = _mm256_fmadd_pd(y_3, scalar5, y_11);
            y_19 = _mm256_fmadd_pd(y_7, scalar4, y_19);
            
            _mm256_store_pd(&temp1[i*n + j], y_12);
            _mm256_store_pd(&temp1[(i + 1)*n + j], y_13);
            _mm256_store_pd(&temp1[(i + 2)*n + j], y_14);
            _mm256_store_pd(&temp1[(i + 3)*n + j], y_15);
            _mm256_store_pd(&temp2[i*n + j], y_16);
            _mm256_store_pd(&temp2[(i + 1)*n + j], y_17);
            _mm256_store_pd(&temp2[(i + 2)*n + j], y_18);
            _mm256_store_pd(&temp2[(i + 3)*n + j], y_19);
        }
    }
    FLOPS_ADD_N2(3)
    FLOPS_ADD_N2(3)
    FLOPS_ADD_N2(3)

    //Second part
    scalar = _mm256_set1_pd(1187353796428800.0);
    scalar2 = _mm256_set1_pd(10559470521600.0);
    scalar3 = _mm256_set1_pd(33522128640.0);
    scalar4 = _mm256_set1_pd(670442572800.0);
    scalar5 = _mm256_set1_pd(129060195264000.0);
    __m256d scalar6 = _mm256_set1_pd(7771770303897600.0);
    for(int i = 0; i < n; i+=4) {
        for(int j = 0; j < n; j+=4) {
            __m256d y = _mm256_load_pd(&A2[i*n + j]);
            __m256d y_1 = _mm256_load_pd(&A2[(i + 1)*n + j]);
            __m256d y_2 = _mm256_load_pd(&A2[(i + 2)*n + j]);
            __m256d y_3 = _mm256_load_pd(&A2[(i + 3)*n + j]);

            __m256d y_4 = _mm256_load_pd(&A4[i*n + j]);
            __m256d y_5 = _mm256_load_pd(&A4[(i + 1)*n + j]);
            __m256d y_6 = _mm256_load_pd(&A4[(i + 2)*n + j]);
            __m256d y_7 = _mm256_load_pd(&A4[(i + 3)*n + j]);

            __m256d y_8 = _mm256_load_pd(&A6[i*n + j]);
            __m256d y_9 = _mm256_load_pd(&A6[(i + 1)*n + j]);
            __m256d y_10 = _mm256_load_pd(&A6[(i + 2)*n + j]);
            __m256d y_11 = _mm256_load_pd(&A6[(i + 3)*n + j]);

            __m256d y_20 = _mm256_load_pd(&temp4[i*n + j]);
            __m256d y_21 = _mm256_load_pd(&temp4[(i + 1)*n + j]);
            __m256d y_22 = _mm256_load_pd(&temp4[(i + 2)*n + j]);
            __m256d y_23 = _mm256_load_pd(&temp4[(i + 3)*n + j]);
            __m256d y_24 = _mm256_load_pd(&temp1[i*n + j]);
            __m256d y_25 = _mm256_load_pd(&temp1[(i + 1)*n + j]);
            __m256d y_26 = _mm256_load_pd(&temp1[(i + 2)*n + j]);
            __m256d y_27 = _mm256_load_pd(&temp1[(i + 3)*n + j]);

            y_20 = _mm256_fmadd_pd(y_8, scalar3, y_20);
            y_21 = _mm256_fmadd_pd(y_9, scalar3, y_21);
            y_22 = _mm256_fmadd_pd(y_10, scalar3, y_22);
            y_23 = _mm256_fmadd_pd(y_11, scalar3, y_23);
            y_24 = _mm256_fmadd_pd(y_8, scalar4, y_24);
            y_25 = _mm256_fmadd_pd(y_9, scalar4, y_25);
            y_26 = _mm256_fmadd_pd(y_10, scalar4, y_26);
            y_27 = _mm256_fmadd_pd(y_11, scalar4, y_27);

            y_20 = _mm256_fmadd_pd(y, scalar, y_20);
            y_20 = _mm256_fmadd_pd(y_4, scalar2, y_20);
            y_21 = _mm256_fmadd_pd(y_1, scalar, y_21);
            y_21 = _mm256_fmadd_pd(y_5, scalar2, y_21);
            y_22 = _mm256_fmadd_pd(y_2, scalar, y_22);
            y_22 = _mm256_fmadd_pd(y_6, scalar2, y_22);
            y_23 = _mm256_fmadd_pd(y_3, scalar, y_23);
            y_23 = _mm256_fmadd_pd(y_7, scalar2, y_23);

            y_24 = _mm256_fmadd_pd(y, scalar6, y_24);
            y_24 = _mm256_fmadd_pd(y_4, scalar5, y_24);
            y_25 = _mm256_fmadd_pd(y_1, scalar6, y_25);
            y_25 = _mm256_fmadd_pd(y_5, scalar5, y_25);
            y_26 = _mm256_fmadd_pd(y_2, scalar6, y_26);
            y_26 = _mm256_fmadd_pd(y_6, scalar5, y_26);
            y_27 = _mm256_fmadd_pd(y_3, scalar6, y_27);
            y_27 = _mm256_fmadd_pd(y_7, scalar5, y_27);

            _mm256_store_pd(&temp4[i*n + j], y_20);
            _mm256_store_pd(&temp4[(i + 1)*n + j], y_21);
            _mm256_store_pd(&temp4[(i + 2)*n + j], y_22);
            _mm256_store_pd(&temp4[(i + 3)*n + j], y_23);
            _mm256_store_pd(&temp1[i*n + j], y_24);
            _mm256_store_pd(&temp1[(i + 1)*n + j], y_25);
            _mm256_store_pd(&temp1[(i + 2)*n + j], y_26);
            _mm256_store_pd(&temp1[(i + 3)*n + j], y_27);
        }
        temp4[i*n + i] += 32382376266240000.0;
        temp4[(i + 1)*n + i + 1] += 32382376266240000.0;
        temp4[(i + 2)*n + i + 2] += 32382376266240000.0;
        temp4[(i + 3)*n + i + 3] += 32382376266240000.0;
        temp1[i*n + i] += 64764752532480000.0;
        temp1[(i + 1)*n + i + 1] += 64764752532480000.0;
        temp1[(i + 2)*n + i + 2] += 64764752532480000.0;
        temp1[(i + 3)*n + i + 3] += 64764752532480000.0;
    }
    FLOPS_ADD_N2(4)
    FLOPS_ADD_N2(4)
    FLOPS_ADD_N2(4)
    FLOPS_ADD_N(2)

    //PROBABBLY 4 MATRIXES HERE

    for(int i = 0; i < n*n - 15; i+=16) {
        __m256d y = _mm256_load_pd(&temp1[i]);
        __m256d y_1 = _mm256_load_pd(&temp1[i + 4]);
        __m256d y_2 = _mm256_load_pd(&temp1[i + 8]);
        __m256d y_3 = _mm256_load_pd(&temp1[i + 12]);

        __m256d y_5 = _mm256_load_pd(&temp2[i]);
        __m256d y_6 = _mm256_load_pd(&temp2[i + 4]);
        __m256d y_7 = _mm256_load_pd(&temp2[i + 8]);
        __m256d y_8 = _mm256_load_pd(&temp2[i + 12]);

        __m256d y_9 = _mm256_add_pd(y, y_5);
        __m256d y_10 = _mm256_add_pd(y_1, y_6);
        __m256d y_11 = _mm256_add_pd(y_2, y_7);
        __m256d y_12 = _mm256_add_pd(y_3, y_8);

        __m256d y_13 = _mm256_sub_pd(y, y_5);
        __m256d y_14 = _mm256_sub_pd(y_1, y_6);
        __m256d y_15 = _mm256_sub_pd(y_2, y_7);
        __m256d y_16 = _mm256_sub_pd(y_3, y_8);

        _mm256_store_pd(&temp4[i], y_9);
        _mm256_store_pd(&temp4[i + 4], y_10);
        _mm256_store_pd(&temp4[i + 8], y_11);
        _mm256_store_pd(&temp4[i + 12], y_12);
        _mm256_store_pd(&temp1[i], y_13);
        _mm256_store_pd(&temp1[i + 4], y_14);
        _mm256_store_pd(&temp1[i + 8], y_15);
        _mm256_store_pd(&temp1[i + 12], y_16);
    }
    FLOPS_ADD_N2(2)

    //DISCOVER THE NUMBER TO USE BLOCKING
    solve_system(temp4, temp1, X, n);

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
