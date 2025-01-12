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

void m_transpose(double *A, double *result, int n){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            result[i*n + j] = A[j*n + i];
        }
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

void luDecomposition(double *Q, int *Piv, int n) {
    double temp;
    int max_row_index;

    FLOPS_ADD_N_(2);
    for (int pivot = 0; pivot < n; pivot++){
        
        //Take the maximum value in the column i and save the index
        max_row_index = pivot;
        for (int j = pivot + 1; j < n; j++){
            if (fabs(Q[j*n + pivot]) > fabs(Q[max_row_index*n + pivot])){
                max_row_index = j;
            }
        }
        FLOPS_ADD_N_(2);

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
        FLOPS_ADD_(n - pivot - 1)

        //Update the rest of the matrix by subtracting the pivot column multiplied by the pivot row
        for (int j = pivot + 1; j < n; j++){
            Q[j*n + pivot] /= Q[pivot*n + pivot];
            for (int k = pivot + 1; k < n; k++){
                Q[j*n + k] -= Q[j*n + pivot] * Q[pivot*n + k];
                FLOPS_ADD_(2)
            }
        }
    }
}

void last_part(double *P, double *Q, int n){
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

    FLOPS_ADD_N2_(1)
    FLOPS_ADD_N3_(2)
}

void solve_system(double *P, double *Q, int n) {
    FLOPS_RESET_
    int *Pivot = malloc(n * sizeof(int));
    double temp, temp2, temp3, temp4;
    int block_size;
    int size;

    luDecomposition(Q, Pivot, n);

    //LAPACKE_dgesv(LAPACK_ROW_MAJOR,n, n, Q, n, Piv, P2, n);
    //LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, Q_T, n, Piv);

    // Swap P rows accordingly to PIvot
    FLOPS_ADD_N_(1)
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
    FLOPS_PRINT_
}