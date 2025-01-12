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

__attribute__((always_inline)) void luDecomposition(double *Q, int *Piv, int n, int block_size, int start_from) {
    double max_value, temp, temp2, temp3, temp4;
    int max_row_index;
    FLOPS_ADD(block_size * 2);
    for (int pivot = start_from; pivot < start_from + block_size; pivot++){
        
        //Take the maximum value in the column i and save the index
        max_value = Q[pivot*n + pivot] > 0 ? Q[pivot*n + pivot] : -Q[pivot*n + pivot];
        max_row_index = pivot;
        for (int j = pivot + 1; j < n; j++){
            temp = Q[j*n + pivot] > 0 ? Q[j*n + pivot] : -Q[j*n + pivot];
            if (temp > max_value){
                max_value = temp;
                max_row_index = j;
            }
        }
        FLOPS_ADD_N(2);

        //Swap the columns, do it to increse ILP
        if (max_row_index != pivot){
            int j;
            for (j = 0; j < n - 3; j+=4){
                temp = Q[pivot*n + j];
                Q[pivot*n + j] = Q[max_row_index*n + j];
                Q[max_row_index*n + j] = temp;
                temp2 = Q[pivot*n + j + 1];
                Q[pivot*n + j + 1] = Q[max_row_index*n + j + 1];
                Q[max_row_index*n + j + 1] = temp2;
                temp3 = Q[pivot*n + j + 2];
                Q[pivot*n + j + 2] = Q[max_row_index*n + j + 2];
                Q[max_row_index*n + j + 2] = temp3;
                temp4 = Q[pivot*n + j + 3];
                Q[pivot*n + j + 3] = Q[max_row_index*n + j + 3];
                Q[max_row_index*n + j + 3] = temp4;
            }
            for (; j < n; j++){
                temp = Q[pivot*n + j];
                Q[pivot*n + j] = Q[max_row_index*n + j];
                Q[max_row_index*n + j] = temp;
            }
        }

        Piv[pivot] = max_row_index;
        double shit = 1.0 / Q[pivot*n + pivot];
        //Divide the pivot column by the pivot element
        FLOPS_ADD(n - pivot - 1)
        for (int j = pivot + 1; j < n; j++){
            Q[j*n + pivot] *= shit;
        }

        //Update the rest of the matrix by subtracting the pivot column multiplied by the pivot row
        int j;
        for (j = pivot + 1; j < n - 3; j+=4){
            __m256d var0 = _mm256_broadcast_sd(Q + j*n + pivot);
            __m256d var1 = _mm256_broadcast_sd(Q + (j + 1)*n + pivot);
            __m256d var2 = _mm256_broadcast_sd(Q + (j + 2)*n + pivot);
            __m256d var3 = _mm256_broadcast_sd(Q + (j + 3)*n + pivot);
            int k;
            for (k = pivot + 1; k < start_from + block_size - 3; k+=4){
                __m256d var4 = _mm256_loadu_pd(Q + pivot*n + k);

                __m256d var8 = _mm256_loadu_pd(Q + j*n + k);
                __m256d var9 = _mm256_loadu_pd(Q + (j + 1)*n + k);
                __m256d var10 = _mm256_loadu_pd(Q + (j + 2)* n + k);
                __m256d var11 = _mm256_loadu_pd(Q + (j + 3)* n + k);
                
                __m256d var12 = _mm256_fnmadd_pd(var0, var4, var8);
                __m256d var13 = _mm256_fnmadd_pd(var1, var4, var9);
                __m256d var14 = _mm256_fnmadd_pd(var2, var4, var10);
                __m256d var15 = _mm256_fnmadd_pd(var3, var4, var11);

                _mm256_storeu_pd(&Q[j*n + k], var12);
                _mm256_storeu_pd(&Q[(j + 1)*n + k], var13);
                _mm256_storeu_pd(&Q[(j + 2)*n + k], var14);
                _mm256_storeu_pd(&Q[(j + 3)*n + k], var15);
                FLOPS_ADD(8)
            }
            for (; k < start_from + block_size; k++){
                Q[j*n + k] -= Q[j*n + pivot] * Q[pivot*n + k];
                Q[(j + 1)*n + k] -= Q[(j + 1)*n + pivot] * Q[pivot*n + k];
                Q[(j + 2)*n + k] -= Q[(j + 2)*n + pivot] * Q[pivot*n + k];
                Q[(j + 3)*n + k] -= Q[(j + 3)*n + pivot] * Q[pivot*n + k];
                FLOPS_ADD(8)
            }
        }
        for (; j < n; j++){
            __m256d var0 = _mm256_broadcast_sd(Q + j*n + pivot);
            int k;
            for (k = pivot + 1; k < start_from + block_size - 3; k+=4){
                __m256d var1 = _mm256_loadu_pd(Q + pivot*n + k);
                __m256d var2 = _mm256_loadu_pd(Q + j*n + k);
                _mm256_storeu_pd(&Q[j*n + k], _mm256_fnmadd_pd(var0, var1, var2));
                FLOPS_ADD(2)
            }
            for (; k < start_from + block_size; k++){
                Q[j*n + k] -= Q[j*n + pivot] * Q[pivot*n + k];
                FLOPS_ADD(2)
            }
        }
    }
}

__attribute__((always_inline))  void last_part(double *P, double *Q, int n){
    // Forward substitution

    for (int j = 0; j < n; j++){
        int k;
        for (k = 0; k < j - 3; k+=4){
            __m256d var0 = _mm256_broadcast_sd(Q + j*n + k);
            __m256d var1 = _mm256_broadcast_sd(Q + j*n + k + 1);
            __m256d var2 = _mm256_broadcast_sd(Q + j*n + k + 2);
            __m256d var3 = _mm256_broadcast_sd(Q + j*n + k + 3);
            for (int i = 0; i < n; i+=4){
                __m256d var4 = _mm256_load_pd(P + k*n + i);
                __m256d var5 = _mm256_load_pd(P + (k + 1)*n + i);
                __m256d var6 = _mm256_load_pd(P + (k + 2)*n + i);
                __m256d var7 = _mm256_load_pd(P + (k + 3)*n + i);

                __m256d var8 = _mm256_load_pd(P + j*n + i);

                __m256d var9 = _mm256_mul_pd(var0, var4);
                var9 = _mm256_fmadd_pd(var1, var5, var9);
                var9 = _mm256_fmadd_pd(var2, var6, var9);
                var9 = _mm256_fmadd_pd(var3, var7, var9);

                _mm256_store_pd(&P[j*n + i], _mm256_sub_pd(var8, var9));
            }
        }
        for (; k < j; k++){
            __m256d var0 = _mm256_broadcast_sd(Q + j*n + k);
            for (int i = 0; i < n; i+=4){
                __m256d var1 = _mm256_load_pd(P + k*n + i);
                __m256d var2 = _mm256_load_pd(P + j*n + i);
                _mm256_store_pd(&P[j*n + i], _mm256_fnmadd_pd(var0, var1, var2));
            }
        }
    }

// Backward substitution
    for (int i = n - 1; i >= 0; i--) {
        int k;
        for (k = n - 1; k > i + 3; k-=4) {
            __m256d var0 = _mm256_broadcast_sd(Q + i*n + k);
            __m256d var1 = _mm256_broadcast_sd(Q + i*n + k - 1);
            __m256d var2 = _mm256_broadcast_sd(Q + i*n + k - 2);
            __m256d var3 = _mm256_broadcast_sd(Q + i*n + k - 3);
            for (int j = 0; j < n; j+=4) {
                __m256d var4 = _mm256_load_pd(P + k*n + j);
                __m256d var5 = _mm256_load_pd(P + (k - 1)*n + j);
                __m256d var6 = _mm256_load_pd(P + (k - 2)*n + j);
                __m256d var7 = _mm256_load_pd(P + (k - 3)*n + j);

                __m256d var8 = _mm256_load_pd(P + i*n + j);

                __m256d var9 = _mm256_mul_pd(var0, var4);
                var9 = _mm256_fmadd_pd(var1, var5, var9);
                var9 = _mm256_fmadd_pd(var2, var6, var9);
                var9 = _mm256_fmadd_pd(var3, var7, var9);

                _mm256_store_pd(&P[i*n + j], _mm256_sub_pd(var8, var9));
            }
        }
        for (; k > i; k--){
            __m256d var0 = _mm256_broadcast_sd(Q + i*n + k);
            for (int j = 0; j < n; j+=4){
                __m256d var1 = _mm256_load_pd(P + k*n + j);
                __m256d var2 = _mm256_load_pd(P + i*n + j);
                _mm256_store_pd(&P[i*n + j], _mm256_fnmadd_pd(var0, var1, var2));
            }
        }
        __m256d var0 = _mm256_broadcast_sd(Q + i*n + i);
        for (int j = 0; j < n; j+=4) {
            __m256d var1 = _mm256_load_pd(P + i*n + j);
            _mm256_store_pd(&P[i*n + j], _mm256_div_pd(var1, var0));
        }
    }

    FLOPS_ADD_N2(1)
    FLOPS_ADD_N3(2)
}

__attribute__((always_inline)) void solve_system(double *P, double *Q, int n) {
    int *Pivot = malloc(n * sizeof(int));
    double temp, temp2, temp3, temp4;
    int block_size;
    int size;

    if (n < 256)
        block_size = n;
    else
        block_size = 128;

    for (int i = 0; i < n; i+=block_size){
        luDecomposition(Q, Pivot, n, block_size, i);

        if (i < n - block_size){
            for (int k = i + 1; k < i + block_size; k++) {
                int j;
                for (j = i; j < k - 3; j+=4) {
                    __m256d var0 = _mm256_broadcast_sd(Q + k*n + j);
                    __m256d var1 = _mm256_broadcast_sd(Q + k*n + j + 1);
                    __m256d var2 = _mm256_broadcast_sd(Q + k*n + j + 2);
                    __m256d var3 = _mm256_broadcast_sd(Q + k*n + j + 3);
                    for (int l = i + block_size; l < n; l+=4) {
                        __m256d var4 = _mm256_load_pd(Q + k*n + l);

                        __m256d var5 = _mm256_load_pd(Q + j*n + l);
                        __m256d var6 = _mm256_load_pd(Q + (j + 1)*n + l);
                        __m256d var7 = _mm256_load_pd(Q + (j + 2)*n + l);
                        __m256d var8 = _mm256_load_pd(Q + (j + 3)*n + l);
                        
                        __m256d var9 = _mm256_mul_pd(var0, var5);
                        var9 = _mm256_fmadd_pd(var1, var6, var9);
                        var9 = _mm256_fmadd_pd(var2, var7, var9);
                        var9 = _mm256_fmadd_pd(var3, var8, var9);

                        _mm256_store_pd(&Q[k*n + l], _mm256_sub_pd(var4, var9));
                        //Q[k*n + l] -= Q[k*n + j] * Q[j*n + l];
                        FLOPS_ADD(8)
                    }
                }
                for (; j < k; j++) {
                    __m256d var0 = _mm256_broadcast_sd(Q + k*n + j);
                    for (int l = i + block_size; l < n; l+=4) {
                        __m256d var1 = _mm256_load_pd(Q + k*n + l);
                        __m256d var2 = _mm256_load_pd(Q + j*n + l);
                        
                        _mm256_store_pd(&Q[k*n + l], _mm256_fnmadd_pd(var0, var2, var1));
                        //Q[k*n + l] -= Q[k*n + j] * Q[j*n + l];
                        FLOPS_ADD(2)
                    }
                }
            }
        

            size = n - i - block_size;
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            size, size, block_size, -1.0, Q + (i + block_size) * n + i , n,
            Q + i + block_size + i*n, n, 1.0, Q + (i + block_size)*n + i + block_size, n);
            FLOPS_ADD(size*size*block_size*2)
        }
    }

    //LAPACKE_dgesv(LAPACK_ROW_MAJOR,n, n, Q, n, Piv, P2, n);
    //LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, Q_T, n, Piv);

    // Swap P rows accordingly to PIvot
    FLOPS_ADD_N(1)
    for (int i = 0; i < n; i++){
        if (Pivot[i] != i){
            int j;
            for (j = 0; j < n - 3; j+=4){
                temp = P[i*n + j];
                P[i*n + j] = P[Pivot[i]*n + j];
                P[Pivot[i]*n + j] = temp;
                temp2 = P[i*n + j + 1];
                P[i*n + j + 1] = P[Pivot[i]*n + j + 1];
                P[Pivot[i]*n + j + 1] = temp2;
                temp3 = P[i*n + j + 2];
                P[i*n + j + 2] = P[Pivot[i]*n + j + 2];
                P[Pivot[i]*n + j + 2] = temp3;
                temp4 = P[i*n + j + 3];
                P[i*n + j + 3] = P[Pivot[i]*n + j + 3];
                P[Pivot[i]*n + j + 3] = temp4;
            }
            for (; j < n; j++){
                temp = P[i*n + j];
                P[i*n + j] = P[Pivot[i]*n + j];
                P[Pivot[i]*n + j] = temp;
            }
        }
    }

    last_part(P, Q, n);
    free(Pivot);
}

int mexp_no_blas(double *A, double *A2, double *A4, double *A6, double *A8, double *A10, double **X, int n){
    FLOPS_RESET
    solve_system(A2, A, n);
    FLOPS_PRINT
}