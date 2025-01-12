#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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


__attribute__((always_inline)) void matrix_vector_multiply(double *A, double *x, double *y, int n){
    double temp[4], temp1[4], temp2[4], temp3[4];
    for (int i = 0; i < n; i+=4){
        __m256d zero =  _mm256_setzero_pd();
        __m256d zero1 = _mm256_setzero_pd();
        __m256d zero2 = _mm256_setzero_pd();
        __m256d zero3 = _mm256_setzero_pd();
        for (int j = 0; j < n; j+=4){
            __m256d var0 = _mm256_load_pd(A + i*n + j);
            __m256d var1 = _mm256_load_pd(A + (i+1)*n + j);
            __m256d var2 = _mm256_load_pd(A + (i+2)*n + j);
            __m256d var3 = _mm256_load_pd(A + (i+3)*n + j);

            __m256d var4 = _mm256_load_pd(x + j);
            
            zero = _mm256_fmadd_pd(var0, var4, zero);
            zero1 = _mm256_fmadd_pd(var1, var4, zero1);
            zero2 = _mm256_fmadd_pd(var2, var4, zero2);
            zero3 = _mm256_fmadd_pd(var3, var4, zero3);
        }

        _mm256_store_pd(temp, zero);
        _mm256_store_pd(temp1, zero1);
        _mm256_store_pd(temp2, zero2);
        _mm256_store_pd(temp3, zero3);
        y[i] = temp[0] + temp[1] + temp[2] + temp[3];
        y[i+1] = temp1[0] + temp1[1] + temp1[2] + temp1[3];
        y[i+2] = temp2[0] + temp2[1] + temp2[2] + temp2[3];
        y[i+3] = temp3[0] + temp3[1] + temp3[2] + temp3[3];
    }
    FLOPS_ADD_N2(2)
}

__attribute__((always_inline)) int is_triangular(double *A, int n) {
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++)
            if (A[i * n + j] != 0.0)
                return 0;
    }
    return 1;
}

__attribute__((always_inline)) void luDecomposition(double *Q, int *Piv, int n) {
    double temp;
    int max_row_index;

    FLOPS_ADD_N(2);
    for (int pivot = 0; pivot < n; pivot++){
        
        //Take the maximum value in the column i and save the index
        max_row_index = pivot;
        for (int j = pivot + 1; j < n; j++){
            if (fabs(Q[j*n + pivot]) > fabs(Q[max_row_index*n + pivot])){
                max_row_index = j;
            }
        }
        FLOPS_ADD_N(2);

        //Swap the columns, do it to increse ILP
        if (max_row_index != pivot){
            int j;
            for (j = 0; j < n; j++){
                temp = Q[pivot*n + j];
                Q[pivot*n + j] = Q[max_row_index*n + j];
                Q[max_row_index*n + j] = temp;
            }
        }

        Piv[pivot] = max_row_index;
        //Divide the pivot column by the pivot element
        FLOPS_ADD(n - pivot - 1)

        //Update the rest of the matrix by subtracting the pivot column multiplied by the pivot row
        for (int j = pivot + 1; j < n; j++){
            Q[j*n + pivot] /= Q[pivot*n + pivot];
            for (int k = pivot + 1; k < n; k++){
                Q[j*n + k] -= Q[j*n + pivot] * Q[pivot*n + k];
                FLOPS_ADD(2)
            }
        }
    }
}

__attribute__((always_inline)) void last_part(double *P, double *Q, int n){
    // Forward substitution

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            for (int k = 0; k < j; k++){
                P[j*n + i] -= P[k*n + i] * Q[j*n + k];
            }
        }
    }

// Backward substitution
    for (int i = n - 1; i >= 0; i--) {
        for (int k = n - 1; k > i; k--) {
            for (int j = 0; j < n; j++) {
                P[i*n + j] -= Q[i*n + k] * P[k*n + j];
            }
        }
        for (int j = 0; j < n; j++) {
            P[i*n + j] /= Q[i*n + i];
        }
    }

    FLOPS_ADD_N2(1)
    FLOPS_ADD_N3(2)
}

__attribute__((always_inline)) void solve_system(double *P, double *Q, int n) {
    int *Pivot = malloc(n * sizeof(int));
    double temp;

    luDecomposition(Q, Pivot, n);

    //LAPACKE_dgesv(LAPACK_ROW_MAJOR,n, n, Q, n, Piv, P2, n);
    //LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, Q_T, n, Piv);

    // Swap P rows accordingly to PIvot
    FLOPS_ADD_N(1)
    for (int i = 0; i < n; i++){
        if (Pivot[i] != i){
            int j;
            for (j = 0; j < n; j++){
                temp = P[i*n + j];
                P[i*n + j] = P[Pivot[i]*n + j];
                P[Pivot[i]*n + j] = temp;
            }
        }
    }

    last_part(P, Q, n);
    free(Pivot);
}

__attribute__((always_inline)) double norm_sum_abs(double *Yg_A, int n) {
    __m256d abs_mask = _mm256_set1_pd(-0.0);
    double est = 0.0;
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
        __m256d var0 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(Yg_A + j));
        __m256d var1 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(Yg_A + j + 4));
        __m256d var2 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(Yg_A + j + 8));
        __m256d var3 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(Yg_A + j + 12));

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
        y_sums_0 += Yg_A[j] > 0 ? Yg_A[j] : -Yg_A[j];
        y_sums_1 += Yg_A[j + 1] > 0 ? Yg_A[j + 1] : -Yg_A[j + 1];
        y_sums_2 += Yg_A[j + 2] > 0 ? Yg_A[j + 2] : -Yg_A[j + 2];
        y_sums_3 += Yg_A[j + 3] > 0 ? Yg_A[j + 3] : -Yg_A[j + 3];
        y_sums_4 += Yg_A[j + 4] > 0 ? Yg_A[j + 4] : -Yg_A[j + 4];
        y_sums_5 += Yg_A[j + 5] > 0 ? Yg_A[j + 5] : -Yg_A[j + 5];
        y_sums_6 += Yg_A[j + 6] > 0 ? Yg_A[j + 6] : -Yg_A[j + 6];
        y_sums_7 += Yg_A[j + 7] > 0 ? Yg_A[j + 7] : -Yg_A[j + 7];
    }
    y_sums_0 = y_sums_0 + y_sums_2 + y_sums_4 + y_sums_6;
    y_sums_1 = y_sums_1 + y_sums_3 + y_sums_5 + y_sums_7;

    FLOPS_ADD_N(4);


    FLOPS_ADD(1);
    if (y_sums_0 < y_sums_1)
    {
        est = y_sums_1;
    }else{
        est = y_sums_0;
    }
    return est;
}


/*****************************
* mexp_basic.c
******************************/

int mexp_no_blas(double *A, double *A2, double *A4, double *A6, double *A8, double *A10, double **X, int n){
    FLOPS_RESET
    
    //Generic variable
    __m256d abs_mask = _mm256_set1_pd(-0.0);
    double row_sum;
    double column_sum;
    double d4, d6, d8, d10;
    double eta1, eta3, eta4, eta5;
    int s;
    double est;
    double s_temp;
    double normA, alpha;
    double ell;
    srand(24);
    
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

    size = n*2*sizeof(double);
    if (size % alignment != 0) {
        size += alignment - (size % alignment);
    }
    double *Yg = (double *)aligned_alloc(alignment, size);
    double *support = (double *)aligned_alloc(alignment, size);

    size = n*sizeof(double);
    if (size % alignment != 0) {
        size += alignment - (size % alignment);
    }
    double *temp3 = (double *)aligned_alloc(alignment, size);
    double *temp5 = (double *)aligned_alloc(alignment, size);
    double *Yg_A = (double *)aligned_alloc(alignment, size);

    int *ind_histg = (int *)malloc((n + 1) * sizeof(int));

    FLOPS_ADD(3)
    //Some precomputed constants
    double u = pow(2, -53);
    int ind_hist_leng = 0;
    double div_n_for_normest = 1.0/n;
    double d6_factor = 1.0/6.0;

    //COMPUTE NORM(A) AND the transpose of the absolute of A to use later for norm
    column_sum = 0;
    double temp_T[8][4];
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

            _mm256_store_pd(temp_T[0], A_vec);
            _mm256_store_pd(temp_T[1], A_vec2);
            _mm256_store_pd(temp_T[2], A_vec3);
            _mm256_store_pd(temp_T[3], A_vec4);
            _mm256_store_pd(temp_T[4], A_vec5);
            _mm256_store_pd(temp_T[5], A_vec6);
            _mm256_store_pd(temp_T[6], A_vec7);
            _mm256_store_pd(temp_T[7], A_vec8);

            for (int k = 0; k < 8; k++){
                tempA[(i + k)*n + j] = temp_T[k][0];
                tempA[(i + k)*n + j + 1] = temp_T[k][1];
                tempA[(i + k)*n + j + 2] = temp_T[k][2];
                tempA[(i + k)*n + j + 3] = temp_T[k][3];
            }
        }
        __m256d max = _mm256_max_pd(y_sum, y_sum_1);

        double max_[4];
        _mm256_store_pd(max_, max);

        max_[0] = max_[0] > max_[1] ? max_[0] : max_[1];
        max_[2] = max_[2] > max_[3] ? max_[2] : max_[3];
        normA = max_[0] > max_[2] ? max_[0] : max_[2];
    }
    for(; i < n; i++){
        for (int j = 0; j < n; j ++){
                tempA[i*n + j] = A[j*n + i] > 0 ? A[j*n + i] : -A[j*n + i];
                column_sum += tempA[i*n + j];
        }
        normA = column_sum > normA ? column_sum : normA;
    }
    FLOPS_ADD_N(8)
    FLOPS_ADD_N2(2)


    /*****************************
    * NORMEST STARTS M=0
    ******************************/

    for (int i = 0; i < n; i++) {
        row_sum = 0.0;
        int j;
        __m256d y_sum = _mm256_setzero_pd();
        __m256d y_sum_1 = _mm256_setzero_pd();
        __m256d y_sum_2 = _mm256_setzero_pd();
        __m256d y_sum_3 = _mm256_setzero_pd();
        for (j = 0; j < n - 15; j+=16) {
            __m256d y = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A2[i * n + j]));
            __m256d y_1 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A2[i * n + j + 4]));
            __m256d y_2 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A2[i * n + j + 8]));
            __m256d y_3 = _mm256_andnot_pd(abs_mask, _mm256_load_pd(&A2[i * n + j + 12]));

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
            row_sum += A2[i * n + j] > 0 ? A2[i * n + j] : -A2[i * n + j];
        }
        long long n = (long long)(row_sum * div_n_for_normest);
        Yg[i * 2] = (double)n;
    }
    FLOPS_ADD_N2(2)
    FLOPS_ADD_N(1)

    for (int i = 1; i < 2; i++) {
        int col_idx = rand() % n;
        for (int j = 0; j < n; j++) {
            Yg[2 * j + i] = A2[j * n + col_idx];
        }
        ind_histg[ind_hist_leng++] = col_idx;
    }

    FLOPS_ADD_N2(4);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 2, n, 1.0, A2, n, Yg, 2, 0.0, support, 2);

    est = norm_sum_abs(Yg, n);
    /*****************************
    * NORMEST ENDS M=2
    ******************************/

    d4 = pow(est, 0.25);

    /*****************************
    * NORMEST STARTS M=3
    ******************************/
    FLOPS_ADD_N2(4);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 2, n, 1.0, A2, n, Yg, 2, 0.0, support, 2);

    est = norm_sum_abs(support, n);

    free(Yg);
    free(support);
    free(ind_histg);
    /*****************************
    * NORMEST ENDS M=3
    ******************************/


    d6 = pow(est, d6_factor);
    eta1 = d4 > d6 ? d4 : d6;
    FLOPS_ADD(3)

   /*****************************
    * NORMEST STARTS M=0
    ******************************/
    for (int i = 0; i < n; i++)
        temp3[i] = 1.0;

    matrix_vector_multiply(tempA, temp3, Yg_A, n);
    matrix_vector_multiply(tempA, Yg_A, temp3, n);
    matrix_vector_multiply(tempA, temp3, Yg_A, n);
    matrix_vector_multiply(tempA, Yg_A, temp3, n);
    matrix_vector_multiply(tempA, temp3, Yg_A, n);
    matrix_vector_multiply(tempA, Yg_A, temp3, n);
    matrix_vector_multiply(tempA, temp3, Yg_A, n);

    est = 0;
    for (int i = 0; i < n; i++){
        est = Yg_A[i] > est ? Yg_A[i] : est;
    }
    FLOPS_ADD_N(1)

    FLOPS_ADD(1)
    if (eta1 <= 0.01495585217958292){
        alpha = est;
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
            solve_system(temp1, temp2, n);
            if (*X != NULL) {
                free(*X);
            }
            *X = temp1;

            free(temp2);
            free(temp3);
            free(temp5);
            free(temp4);
            free(tempA);
            free(Yg_A);
            FLOPS_PRINT;
            return 1;
        }
    }
    
    FLOPS_ADD(1)

    matrix_vector_multiply(tempA, Yg_A, temp3, n);
    matrix_vector_multiply(tempA, temp3, Yg_A, n);
    matrix_vector_multiply(tempA, Yg_A, temp5, n);
    matrix_vector_multiply(tempA, temp5, Yg_A, n);

    est = 0;
    for (int i = 0; i < n; i++){
        est = Yg_A[i] > est ? Yg_A[i] : est;
    }
    FLOPS_ADD_N(1)

    if (eta1 <= 0.2539398330063230){
        alpha = est;
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
            solve_system(temp1, temp2, n);

            if (*X != NULL) {
                free(*X);
            }
            *X = temp1;

            free(temp2);
            free(temp3);
            free(temp4);
            free(temp5);
            free(tempA);
            free(Yg_A);
            FLOPS_PRINT;
            return 2;
        }
    }

    FLOPS_ADD(1)

    est = 0;
    for (int i = 0; i < n; i++){
        est = temp3[i] > est ? temp3[i] : est;
    }
    FLOPS_ADD_N(1)

    d8 = pow(est, 0.125);
    eta3 = d6 > d8 ? d6 : d8;
    FLOPS_ADD(2)

    matrix_vector_multiply(tempA, Yg_A, temp3, n);
    matrix_vector_multiply(tempA, temp3, Yg_A, n);
    matrix_vector_multiply(tempA, Yg_A, temp3, n);
    matrix_vector_multiply(tempA, temp3, Yg_A, n);

    est = 0;
    for (int i = 0; i < n; i++){
        est = Yg_A[i] > est ? Yg_A[i] : est;
    }
    FLOPS_ADD_N(1)

    if (eta3 <= 0.9504178996162932){
        alpha = est;
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
            solve_system(temp1, temp2, n);
            if (*X != NULL) {
                free(*X);
            }
            *X = temp1;
            free(temp2);
            free(temp3);
            free(temp4);
            free(tempA);
            free(temp5);
            free(Yg_A);
            FLOPS_PRINT;
            return 3;
        }
    }

    matrix_vector_multiply(tempA, Yg_A, temp3, n);
    matrix_vector_multiply(tempA, temp3, Yg_A, n);
    matrix_vector_multiply(tempA, Yg_A, temp3, n);
    matrix_vector_multiply(tempA, temp3, Yg_A, n);

    est = 0;
    for (int i = 0; i < n; i++){
        est = Yg_A[i] > est ? Yg_A[i] : est;
    }
    FLOPS_ADD_N(1)

    FLOPS_ADD(1)
    if (eta3 <= 2.097847961257068){
        alpha = est;
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
            solve_system(temp1, temp2, n);
            if (*X != NULL) {
                free(*X);
            }
            *X = temp1;
            free(temp2);
            free(temp3);
            free(temp4);
            free(tempA);
            free(temp5);
            free(Yg_A);
            FLOPS_PRINT;
            return 3;
        }
    }
   
    // Up to line 27
    // Continues below
    FLOPS_ADD_N(1);
    est = 0;
    for (int i = 0; i < n; i++){
        est = temp5[i] > est ? temp5[i] : est;
    }
    d10 = pow(est, 0.1);
    free(temp5);
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

    eta3 = pow(2, -s*19);
    FLOPS_ADD(2)
    for (int i = 0; i < n; i++)
        Yg_A[i] = Yg_A[i] * eta3;
    FLOPS_ADD_N(1)
        
    eta1 = pow(2, -s);
    normA = normA * eta1;
    for (int i = 0; i < n*n; i++)
        tempA[i] = tempA[i] * eta1;
    FLOPS_ADD_N2(1)
    FLOPS_ADD(2)

    matrix_vector_multiply(tempA, Yg_A, temp3, n);
    matrix_vector_multiply(tempA, temp3, Yg_A, n);
    matrix_vector_multiply(tempA, Yg_A, temp3, n);
    matrix_vector_multiply(tempA, temp3, Yg_A, n);
    matrix_vector_multiply(tempA, Yg_A, temp3, n);
    matrix_vector_multiply(tempA, temp3, Yg_A, n);
    matrix_vector_multiply(tempA, Yg_A, temp3, n);
    matrix_vector_multiply(tempA, temp3, Yg_A, n);

    est = 0;
    for (int i = 0; i < n; i++){
        est = Yg_A[i] > est ? Yg_A[i] : est;
    }
    FLOPS_ADD_N(1)

    alpha = ceil(log2(est/(u * 113250775606021113483283660800000000.0 * normA)) * 0.03846153846);
    ell = alpha > 0 ? alpha : 0;
    s += ell;
    FLOPS_ADD(8)

    free(Yg_A);
    free(temp3);
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
    solve_system(temp4, temp1, n);
    if (*X != NULL) {
        free(*X);
    }
    *X = temp4;

    if (s==0){ 
        free(temp1);
        free(temp2);
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
        if (s % 2 == 1)
            free(temp4);
        free(tempA);
        FLOPS_PRINT
        return 5;
    } else {
        free(temp1);
        free(temp2);
        free(tempA);
        FLOPS_PRINT
        return 6;
    }
}