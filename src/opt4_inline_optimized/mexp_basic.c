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

/*****************************
* utils.c
******************************/

__attribute__((always_inline)) int is_triangular(double *A, int n) {
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++)
            if (A[i * n + j] != 0.0)
                return 0;
    }
    return 1;
}

__attribute__((always_inline)) bool invert_matrix(const double *input, double *output, int n) {
    size_t alignment = 32;
    size_t size = n * n * sizeof(double);
    if (size % alignment != 0) {
        size += alignment - (size % alignment);
    }
    double *A = (double *)aligned_alloc(alignment, size);
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
    size_t alignment = 32;
    size_t size = n * n * sizeof(double);
    if (size % alignment != 0) {
        size += alignment - (size % alignment);
    }
    double *Q_inv = (double *)aligned_alloc(alignment, size);
    if (*X == NULL) {
        *X = (double *)aligned_alloc(alignment, size);
    }
    FLOPS_ADD(6)
    if (!invert_matrix(Q, Q_inv, n)) {
        printf("Error: matrix is singular and cannot be inverted.\n");
        return;
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, Q_inv, n, P, n, 0.0, *X, n);
    FLOPS_ADD_N3(2)
    free(Q_inv);
}

__attribute__((always_inline)) bool update_vec(double* Y, int n, int* ind_hist, int ind_hist_len, double* A){
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
        return false;
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
            } while (is_in_history && limit < 5);
            FLOPS_ADD(2);
            ind_hist[ind_hist_len++] = pick;
            for (int j = 0; j < n; j++) {
                Y[2 * j + i] = A[j * n + pick];
            }
        }
    }
    return true;
}

__attribute__((always_inline)) double norm_sum_abs(double *Yg_A, int n) {
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
    for (j = 0; j < 2*n - 7; j+=8) {
        y_sums_0 += Yg_A[j] > 0 ? Yg_A[j] : -Yg_A[j];
        y_sums_1 += Yg_A[j + 1] > 0 ? Yg_A[j + 1] : -Yg_A[j + 1];
        y_sums_2 += Yg_A[j + 2] > 0 ? Yg_A[j + 2] : -Yg_A[j + 2];
        y_sums_3 += Yg_A[j + 3] > 0 ? Yg_A[j + 3] : -Yg_A[j + 3];
        y_sums_4 += Yg_A[j + 4] > 0 ? Yg_A[j + 4] : -Yg_A[j + 4];
        y_sums_5 += Yg_A[j + 5] > 0 ? Yg_A[j + 5] : -Yg_A[j + 5];
        y_sums_6 += Yg_A[j + 6] > 0 ? Yg_A[j + 6] : -Yg_A[j + 6];
        y_sums_7 += Yg_A[j + 7] > 0 ? Yg_A[j + 7] : -Yg_A[j + 7];
    }
    for (; j < 2*n; j+=2) {
        y_sums_0 += Yg_A[j] > 0 ? Yg_A[j] : -Yg_A[j];
        y_sums_1 += Yg_A[j + 1] > 0 ? Yg_A[j + 1] : -Yg_A[j + 1];
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

int mexp(double *A, double **X, int n){
    FLOPS_RESET
    
    //Generic variable
    double row_sum, row_sum2, row_sum3, row_sum4, row_sum5, row_sum6, row_sum7, row_sum8;
    double column_sum, column_sum2, column_sum3, column_sum4, column_sum5, column_sum6, column_sum7, column_sum8;
    double *A2, *A4, *A6, *A8;
    double d4, d6, d8, d10;
    double eta1, eta3, eta4, eta5;
    int s;
    double est;
    double est_old = 0;
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
    A2 = NULL;
    A4 = NULL;
    A6 = NULL;
    A8 = NULL;

    FLOPS_ADD(3)
    //Some precomputed constants
    double u = pow(2, -53);
    int ind_hist_leng = 0;
    double div_n_for_normest = 1.0/n;
    double d6_factor = 1.0/6.0;

    //COMPUTE NORM(A) AND the transpose of the absolute of A to use later for norm
    column_sum = 0.0;
    column_sum2 = 0.0;
    column_sum3 = 0.0;
    column_sum4 = 0.0;
    column_sum5 = 0.0;
    column_sum6 = 0.0;
    column_sum7 = 0.0;
    column_sum8 = 0.0;
    int i;
    for(i = 0; i < n - 127; i+=128)
        for (int j = 0; j < n - 127; j +=128)
            for (int k = i; k < i+128; k++)
                for (int l = j; l < j+128; l++)
                    tempA[l*n + k] = A[k*n + l] > 0 ? A[k*n + l] : -A[k*n + l];
    for(; i < n; i++)
        for (int j = 0; j < n; j ++)
                tempA[j*n + i] = A[i*n + j] > 0 ? A[i*n + j] : -A[i*n + j];
    for (int j = 0; j < n; j++) {
        int i;
        for (i = 0; i < n - 7; i+=8) {
            column_sum += A[i * n + j] > 0 ? A[i * n + j] : -A[i * n + j];
            column_sum2 += A[(i + 1) * n + j] > 0 ? A[(i + 1) * n + j] : -A[(i + 1) * n + j];
            column_sum3 += A[(i + 2) * n + j] > 0 ? A[(i + 2) * n + j] : -A[(i + 2) * n + j];
            column_sum4 += A[(i + 3) * n + j] > 0 ? A[(i + 3) * n + j] : -A[(i + 3) * n + j];
            column_sum5 += A[(i + 4) * n + j] > 0 ? A[(i + 4) * n + j] : -A[(i + 4) * n + j];
            column_sum6 += A[(i + 5) * n + j] > 0 ? A[(i + 5) * n + j] : -A[(i + 5) * n + j];
            column_sum7 += A[(i + 6) * n + j] > 0 ? A[(i + 6) * n + j] : -A[(i + 6) * n + j];
            column_sum8 += A[(i + 7) * n + j] > 0 ? A[(i + 7) * n + j] : -A[(i + 7) * n + j];
        }
        for (; i < n; i++) {
            column_sum += A[i * n + j] > 0 ? A[i * n + j] : -A[i * n + j];
        }
        column_sum = column_sum + column_sum2 + column_sum3 + column_sum4 + column_sum5 + column_sum6 + column_sum7 + column_sum8;
        normA = column_sum > normA ? column_sum : normA;
    }
    FLOPS_ADD(1)
    FLOPS_ADD_N(1)
    FLOPS_ADD_N2(3)

    // Calculate the size of the memory to allocate
    size = n * n * sizeof(double);
    // Ensure the size is a multiple of the alignment
    if (size % alignment != 0) {
        size += alignment - (size % alignment);
    }
    // Allocate the memory
    A2 = (double *)aligned_alloc(alignment, size);
    if (A2 == NULL) {
        printf("Error: memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A, n, A, n, 0.0, A2, n);
    FLOPS_ADD_N3(2)
    FLOPS_ADD(5)

    /*****************************
    * NORMEST STARTS M=0
    ******************************/

    for (int i = 0; i < n; i++) {
        row_sum = 0.0;
        row_sum2 = 0.0;
        row_sum3 = 0.0;
        row_sum4 = 0.0;
        row_sum5 = 0.0;
        row_sum6 = 0.0;
        row_sum7 = 0.0;
        row_sum8 = 0.0;
        int j;
        for (j = 0; j < n - 7; j+=8) {
            row_sum += A2[i * n + j] > 0 ? A2[i * n + j] : -A2[i * n + j];
            row_sum2 += A2[i * n + j + 1] > 0 ? A2[i * n + j + 1] : -A2[i * n + j + 1];
            row_sum3 += A2[i * n + j + 2] > 0 ? A2[i * n + j + 2] : -A2[i * n + j + 2];
            row_sum4 += A2[i * n + j + 3] > 0 ? A2[i * n + j + 3] : -A2[i * n + j + 3];
            row_sum5 += A2[i * n + j + 4] > 0 ? A2[i * n + j + 4] : -A2[i * n + j + 4];
            row_sum6 += A2[i * n + j + 5] > 0 ? A2[i * n + j + 5] : -A2[i * n + j + 5];
            row_sum7 += A2[i * n + j + 6] > 0 ? A2[i * n + j + 6] : -A2[i * n + j + 6];
            row_sum8 += A2[i * n + j + 7] > 0 ? A2[i * n + j + 7] : -A2[i * n + j + 7];
        }
        for (; j < n; j++) {
            row_sum += A2[i * n + j] > 0 ? A2[i * n + j] : -A2[i * n + j];
        }
        row_sum += row_sum2 + row_sum3 + row_sum4 + row_sum5 + row_sum6 + row_sum7 + row_sum8;
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
        Yg_A[i] = 1.0;

    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp3, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp3, 1, 0.0, Yg_A, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp3, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp3, 1, 0.0, Yg_A, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp3, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp3, 1, 0.0, Yg_A, 1);

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
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    temp1[i*n + j] = A2[i*n + j];
                    temp2[i*n + j] = A2[i*n + j] * 12.0;
                    if (i == j) {
                        temp1[i*n + j] += 60.0;
                        temp2[i*n + j] += 120.0;
                    }
                }
            }
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(3)
        
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A, n, temp1, n, 0.0, temp4, n);
            FLOPS_ADD_N3(2)

            for(int i = 0; i < n*n; i++) {
                temp1[i] = temp2[i] + temp4[i];
                temp2[i] = temp2[i] - temp4[i];
            }
            FLOPS_ADD_N2(2)

            //DISCOVER THE NUMBER TO USE BLOCKING
            solve_system(temp1, temp2, X, n);

            free(temp1);
            free(temp2);
            free(temp3);
            free(temp5);
            free(temp4);
            free(Yg_A);
            free(tempA);
            free(A2);
            FLOPS_PRINT;
            return 1;
        }
    }
    
    size = n * n * sizeof(double);
    if (size % alignment != 0) {
        size += alignment - (size % alignment);
    }
    // Allocate the memory
    A4 = (double *)aligned_alloc(alignment, size);
    if (A4 == NULL) {
        printf("Error: memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A2, n, A2, n, 0.0, A4, n);
    FLOPS_ADD_N3(2)
    FLOPS_ADD(6)

    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp3, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp3, 1, 0.0, Yg_A, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp5, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp5, 1, 0.0, Yg_A, 1);

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
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    temp1[i*n + j] = A2[i*n + j] * 420;
                    temp2[i*n + j] = A2[i*n + j] * 3360.0;
                    if (i == j) {
                        temp1[i*n + j] += 15120.0;
                        temp2[i*n + j] += 30240.0;
                    }
                }
            }
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(3)

            for(int i = 0; i < n*n; i++) {
                temp1[i] += A4[i];
                temp2[i] += A4[i] * 30.0;
            }
            FLOPS_ADD_N2(3)

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A, n, temp1, n, 0.0, temp4, n);
            FLOPS_ADD_N3(2)

            for(int i = 0; i < n*n; i++) {
                temp1[i] = temp2[i] + temp4[i];
                temp2[i] = temp2[i] - temp4[i];
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
            free(A2);
            free(Yg_A);
            free(A4);
            FLOPS_PRINT;
            return 2;
        }
    }

    size = n * n * sizeof(double);
    if (size % alignment != 0) {
        size += alignment - (size % alignment);
    }
    // Allocate the memory
    A6 = (double *)aligned_alloc(alignment, size);
    if (A6 == NULL) {
        printf("Error: memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A4, n, A2, n, 0.0, A6, n);
    FLOPS_ADD_N3(2)
    FLOPS_ADD(6)

    est = 0;
    for (int i = 0; i < n; i++){
        est = temp3[i] > est ? temp3[i] : est;
    }
    FLOPS_ADD_N(1)

    d8 = pow(est, 0.125);
    eta3 = d6 > d8 ? d6 : d8;
    FLOPS_ADD(2)

    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp3, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp3, 1, 0.0, Yg_A, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp3, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp3, 1, 0.0, Yg_A, 1);

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
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    temp1[i*n + j] = A2[i*n + j] * 277200;
                    temp2[i*n + j] = A2[i*n + j] * 1995840.0;
                    if (i == j) {
                        temp1[i*n + j] += 8648640.0;
                        temp2[i*n + j] += 17297280.0;
                    }
                }
            }
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(3)

            for(int i = 0; i < n*n; i++) {
                temp1[i] += A4[i] * 1512.0;
                temp2[i] += A4[i] * 25200.0;
            }
            FLOPS_ADD_N2(4)

            for(int i = 0; i < n*n; i++) {
                temp1[i] += A6[i];
                temp2[i] += A6[i] * 56.0;
            }
            FLOPS_ADD_N2(3)

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A, n, temp1, n, 0.0, temp4, n);
            FLOPS_ADD_N3(2)

            for(int i = 0; i < n*n; i++) {
                temp1[i] = temp2[i] + temp4[i];
                temp2[i] = temp2[i] - temp4[i];
            }
            FLOPS_ADD_N2(2)

            //DISCOVER THE NUMBER TO USE BLOCKING
            solve_system(temp1, temp2, X, n);
            free(temp1);
            free(temp2);
            free(temp3);
            free(temp4);
            free(tempA);
            free(A2);
            free(A4);
            free(Yg_A);
            free(temp5);
            free(A6);
            FLOPS_PRINT;
            return 3;
        }
    }

    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp3, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp3, 1, 0.0, Yg_A, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp3, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp3, 1, 0.0, Yg_A, 1);

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
            size = n * n * sizeof(double);
            if (size % alignment != 0) {
                size += alignment - (size % alignment);
            }
            // Allocate the memory
            A8 = (double *)aligned_alloc(alignment, size);
            if (A8 == NULL) {
                printf("Error: memory allocation failed.\n");
                exit(EXIT_FAILURE);
            }
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        n, n, n, 1.0, A4, n, A4, n, 0.0, A8, n);
            FLOPS_ADD_N3(2)
            FLOPS_ADD(5)

            //6 matrixes
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    temp1[i*n + j] = A2[i*n + j] * 302702400;
                    temp2[i*n + j] = A2[i*n + j] * 2075673600.0;
                    if (i == j) {
                        temp1[i*n + j] += 8821612800.0;
                        temp2[i*n + j] += 17643225600.0;
                    }
                }
            }
            FLOPS_ADD_N2(2)
            FLOPS_ADD_N(3)

            for(int i = 0; i < n*n; i++) {
                temp1[i] += A4[i] * 2162160.0;
                temp2[i] += A4[i] * 30270240.0;
            }
            FLOPS_ADD_N2(4)

            for(int i = 0; i < n*n; i++) {
                temp1[i] += A6[i]*3960.0;
                temp2[i] += A6[i] * 110880.0;
            }
            FLOPS_ADD_N2(4)

            for(int i = 0; i < n*n; i++) {
                temp1[i] += A8[i];
                temp2[i] += A8[i] * 90.0;
            }

            FLOPS_ADD_N2(3)

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A, n, temp1, n, 0.0, temp4, n);
            FLOPS_ADD_N3(2)

            for(int i = 0; i < n*n; i++) {
                temp1[i] = temp2[i] + temp4[i];
                temp2[i] = temp2[i] - temp4[i];
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
            free(Yg_A);
            free(A2);
            free(A4);
            free(A6);
            free(A8);
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
    for (int i = 0; i < n; i++)
        Yg_A[i] = Yg_A[i] * eta3;
    FLOPS_ADD_N(1)
    FLOPS_ADD(2)


    eta1 = pow(2, -s);
    normA = normA * eta1;
    for (int i = 0; i < n*n; i++)
        tempA[i] = tempA[i] * eta1;
    FLOPS_ADD_N2(1)
    FLOPS_ADD(2)

    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp3, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp3, 1, 0.0, Yg_A, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp3, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp3, 1, 0.0, Yg_A, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp3, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp3, 1, 0.0, Yg_A, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, Yg_A, 1, 0.0, temp3, 1);
    FLOPS_ADD_N2(2);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, 1, n, 1.0, tempA, n, temp3, 1, 0.0, Yg_A, 1);

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
    for (int i = 0; i < n*n; i++){
        tempA[i] = A[i] * value;
        A2[i] = A2[i] * value2;
        A4[i] = A4[i] * value3;
        A6[i] = A6[i] * value4;
    }
    FLOPS_ADD_N2(4)
    FLOPS_ADD(4)

    //P_M & Q_M
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            temp1[i*n + j] = A2[i*n + j] * 40840800.0 + A4[i*n + j] * 16380.0;
            temp2[i*n + j] = A2[i*n + j] * 1323241920.0 + A4[i*n + j] * 960960.0;
        }
    }
    FLOPS_ADD_N2(3)
    FLOPS_ADD_N2(3)

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            temp1[i*n + j] += A6[i*n + j];
            temp2[i*n + j] += A6[i*n + j] * 182.0;
        }
    }
    FLOPS_ADD_N2(3)
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A6, n, temp1, n, 0.0, temp4, n);  //Temp1 --> Temp4
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A6, n, temp2, n, 0.0, temp1, n);  //Temp2 --> Temp1
    FLOPS_ADD_N3(4)

    //Second part
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            temp4[i*n + j] += A2[i*n + j] * 1187353796428800.0;
            temp1[i*n + j] += A2[i*n + j] * 7771770303897600.0;
        }
    }
    FLOPS_ADD_N2(4)

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            temp4[i*n + j] += A4[i*n + j] * 10559470521600.0;
            temp1[i*n + j] += A4[i*n + j] * 129060195264000.0;
        }
    }
    FLOPS_ADD_N2(4)

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            temp4[i*n + j] += A6[i*n + j] * 33522128640.0;
            temp1[i*n + j] += A6[i*n + j] * 670442572800.0;
        }
    }
    FLOPS_ADD_N2(4)

    for (int i = 0; i < n; i++){
        temp4[i*n + i] += 32382376266240000.0;
        temp1[i*n + i] += 64764752532480000.0;
    }
    FLOPS_ADD_N(2)

    //PROBABBLY 4 MATRIXES HERE

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,      //Temp2 --> Temp1
                n, n, n, 1.0, tempA, n, temp4, n, 0.0, temp2, n);  //Temp1 --> Temp2
    FLOPS_ADD_N3(2)

    for(int i = 0; i < n*n; i++) {
        temp4[i] = temp1[i] + temp2[i];
        temp1[i] = temp1[i] - temp2[i];
    }
    FLOPS_ADD_N2(2)

    //DISCOVER THE NUMBER TO USE BLOCKING
    solve_system(temp4, temp1, X, n);

    if (s==0){ 
        free(temp1);
        free(temp2);
        free(temp4);
        free(tempA);
        free(A2);
        free(A4);
        free(A6);
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

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, *X, n, *X, n, 0.0, temp1, n);
            FLOPS_ADD_N3(2)

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

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, temp1, n, temp1, n, 0.0, *X, n);
            FLOPS_ADD_N3(2)

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

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, *X, n, *X, n, 0.0, temp1, n);
            FLOPS_ADD_N3(2)

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
        free(A2);
        free(A4);
        free(A6);
        free(tempA);
        FLOPS_PRINT
        return 5;
    } else {
        int i;
        for (i = 0; i < s-1; i+=2) {
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, *X, n, *X, n, 0.0, temp1, n);
            FLOPS_ADD_N3(2)
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, temp1, n, temp1, n, 0.0, *X, n);
            FLOPS_ADD_N3(2)
        }
        for (; i < s; i++){
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, *X, n, *X, n, 0.0, temp1, n);
            FLOPS_ADD_N3(2)
            *X = temp1;
        }
        if (s % 2 == 0)
            free(temp1);
        free(temp2);
        free(temp4);
        free(A2);
        free(A4);
        free(A6);
        free(tempA);
        FLOPS_PRINT
        return 6;
    }
}
