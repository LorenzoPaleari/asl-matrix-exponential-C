#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

void luDecomposition(double *Q, int *Piv, int n) {
    double max_value, temp, temp2, temp3, temp4;
    int max_row_index;

    FLOPS_ADD_N(2);
    for (int pivot = 0; pivot < n; pivot++){
        
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
            for (j = 0; j < n; j++){
                temp = Q[pivot*n + j];
                Q[pivot*n + j] = Q[max_row_index*n + j];
                Q[max_row_index*n + j] = temp;
            }
        }

        Piv[pivot] = max_row_index;
        double shit = 1.0 / Q[pivot*n + pivot];
        //Divide the pivot column by the pivot element
        FLOPS_ADD(n - pivot - 1)

        //Update the rest of the matrix by subtracting the pivot column multiplied by the pivot row
        for (int j = pivot + 1; j < n; j++){
            Q[j*n + pivot] *= shit;
            for (int k = pivot + 1; k < n; k++){
                Q[j*n + k] -= Q[j*n + pivot] * Q[pivot*n + k];
                FLOPS_ADD(2)
            }
        }
    }
}

void last_part(double *P, double *Q, int n){
    // Forward substitution

    for (int j = 0; j < n; j++){
        for (int k = 0; k < j; k++){
            for (int i = 0; i < n; i++){
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

void solve_system(double *P, double *Q, int n) {
    int *Pivot = malloc(n * sizeof(int));
    double temp, temp2, temp3, temp4;
    int block_size;
    int size;

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

int mexp_no_blas(double *A, double *A2, double *A4, double *A6, double *A8, double *A10, double **X, int n){
    FLOPS_RESET
    solve_system(A2, A, n);
    FLOPS_PRINT
}