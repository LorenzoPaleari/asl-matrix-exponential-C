#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

#ifdef FLOPS
double flops_opt__ = 0;
double flops_opt_N__ = 0;
double flops_opt_N2__ = 0;
double flops_opt_N3__ = 0;
#define FLOPS_RESET__ flops_opt__ = 0; flops_opt_N__ = 0; flops_opt_N2__ = 0; flops_opt_N3__ = 0;
#define FLOPS_ADD__(x) flops_opt__ += x;
#define FLOPS_ADD_N__(x) flops_opt_N__ += x;
#define FLOPS_ADD_N2__(x) flops_opt_N2__ += x;
#define FLOPS_ADD_N3__(x) flops_opt_N3__ += x;
#define FLOPS_PRINT__ printf("NORMS \nFlops: %f\n", flops_opt__); printf("Flops N: %f\n", flops_opt_N__); printf("Flops N^2: %f\n", flops_opt_N2__); printf("Flops N^3: %f\n", flops_opt_N3__);
#else
#define FLOPS_RESET__
#define FLOPS_ADD__(x)
#define FLOPS_ADD_N__(x)
#define FLOPS_ADD_N2__(x)
#define FLOPS_ADD_N3__(x)
#define FLOPS_PRINT__
#endif

double _norm1est(double *A, int n, int m) {
    FLOPS_RESET__;
    int max_iter = 5;
    srand(24);
    double est_old = 0;
    int *ind_hist = (int *)malloc((n + 1) * sizeof(int));
    int ind_hist_len = 0;
    double *Y = (double *)malloc(n * 2 * sizeof(double));

    // Initialize Y
    double div_n = 1.0/n;
    FLOPS_ADD__(1);
    for (int i = 0; i < n; i++) {
        double row_sum = 0;
        for (int j = 0; j < n; j++) {
            row_sum += A[i * n + j];
        }
        long long n = (long long)(row_sum * div_n);
        Y[i * 2] = (double)n;
    }
    FLOPS_ADD_N2__(1)
    FLOPS_ADD_N__(1)

    // Random columns
    for (int i = 1; i < 2; i++) {
        int col_idx = rand() % n;
        for (int j = 0; j < n; j++) {
            Y[2 * j + i] = A[j * n + col_idx];
        }
        ind_hist[ind_hist_len++] = col_idx;
    }
    FLOPS_ADD_N__(2)


    for (int k = 0; k < max_iter; k++) {
        
        for (int _ = 0; _ < m - 1; _++) {
                FLOPS_ADD_N2__(4);
                FLOPS_ADD__(5);
                matrix_multiply_blas(n, 2, n, A, Y, &Y);
        }
        
        
        double y_sums_0 = 0;
        double y_sums_1 = 0;
        FLOPS_ADD_N__(2);
        for (int j = 0; j < n; j++) {
           y_sums_0 += Y[2 * j] > 0 ? Y[2 * j] : -Y[2 * j];
        }
        for (int j = 0; j < n; j++) {
           y_sums_1 += Y[2 * j + 1] > 0 ? Y[2 * j + 1] : -Y[2 * j + 1];
        }
        FLOPS_ADD_N__(2);
        
        double est; 
        FLOPS_ADD__(1);
        if (y_sums_0 < y_sums_1)
        {
            est = y_sums_1;
        }else{
            est = y_sums_0;
        }
        

        FLOPS_ADD__(1);
        if (est <= est_old) {
            est = est_old;
            break;
        }
        est_old = est;

    bool all_in_hist = true;
    for (int i = 0; i < 2; i++) {
        bool in_hist = false;
        for (int j = 0; j < ind_hist_len; j++) {
            FLOPS_ADD__(1);
            if (ind_hist[j] == i) {
                in_hist = true;
                break;
            }
        }
        FLOPS_ADD__(1);
        if (!in_hist) {
            all_in_hist = false;
            break;
        }
    }
    FLOPS_ADD__(1);
    if (all_in_hist) {
        break;
    } else {
        int pick;
        for (int i = 0; i < 2; i++) {
            int is_in_history;
            int limit = 0;
            do {
                is_in_history = 0;
                FLOPS_ADD_N__(1);
                pick = rand() % n;
                FLOPS_ADD__(1);
                limit++;
                for (int j = 0; j < ind_hist_len; j++) {
                    FLOPS_ADD__(1);
                    if (ind_hist[j] == pick) {
                        is_in_history = 1;
                        break;
                    }
                }
            } while (is_in_history && limit < max_iter);
            FLOPS_ADD__(2);
            ind_hist[ind_hist_len++] = pick;
            for (int j = 0; j < n; j++) {
                Y[2 * j + i] = A[j * n + pick];
            }
        }
    }
}
free(ind_hist);
free(Y);
FLOPS_PRINT__
return est_old;
}

int main2() {
    unsigned int seed = 12;
    srand(seed);

    int n = 20;
    double A[20][20] = {
        {0.54, 0.91, 0.63, 0.27, 0.97, 0.62, 0.71, 0.48, 0.58, 0.95, 0.49, 0.10, 0.78, 0.34, 0.87, 0.15, 0.26, 0.13, 0.47, 0.85},
        {0.22, 0.16, 0.42, 0.67, 0.51, 0.27, 0.86, 0.18, 0.97, 0.35, 0.74, 0.91, 0.19, 0.25, 0.38, 0.82, 0.60, 0.69, 0.49, 0.45},
        {0.43, 0.94, 0.52, 0.19, 0.38, 0.71, 0.82, 0.46, 0.57, 0.13, 0.89, 0.12, 0.68, 0.24, 0.47, 0.25, 0.26, 0.43, 0.77, 0.15},
        {0.92, 0.86, 0.72, 0.37, 0.01, 0.62, 0.81, 0.48, 0.18, 0.55, 0.49, 0.90, 0.28, 0.34, 0.17, 0.95, 0.86, 0.23, 0.47, 0.55},
        {0.02, 0.66, 0.42, 0.27, 0.51, 0.77, 0.36, 0.18, 0.37, 0.35, 0.24, 0.61, 0.19, 0.85, 0.68, 0.32, 0.20, 0.49, 0.79, 0.45},
        {0.53, -0.61, 0.63, 0.37, 0.97, 0.92, 0.71, 0.48, 0.28, 0.15, 0.59, 0.70, 0.78, 0.84, 0.87, 0.65, 0.46, 0.33, 0.87, 0.25},
        {0.72, 0.26, 0.42, 0.97, 0.81, 0.57, 0.86, 0.78, 0.47, 0.35, 0.74, 0.21, 0.69, 0.25, 0.38, 0.62, 0.60, 0.79, 0.29, 0.35},
        {0.23, 0.16, 0.62, 0.67, 0.51, 0.47, 0.56, 0.58, 0.47, 0.55, 0.39, 0.11, 0.19, 0.85, 0.48, 0.12, 0.26, 0.43, 0.77, 0.95},
        {0.43, 0.84, 0.92, 0.19, 0.38, 0.11, 0.82, 0.76, 0.17, 0.73, 0.89, 0.32, 0.28, 0.54, 0.47, 0.25, 0.76, 0.63, 0.17, 0.65},
        {0.62, 0.56, 0.32, 0.67, 0.01, 0.62, 0.41, 0.58, 0.58, 0.55, 0.49, 0.40, 0.98, 0.64, 0.37, 0.95, 0.46, 0.53, 0.77, 0.25},
        {0.92, 0.36, 0.72, 0.87, 0.91, 0.27, 0.16, 0.78, 0.97, 0.95, 0.74, 0.51, 0.29, 0.65, 0.78, 0.82, 0.30, 0.69, 0.19, 0.45},
        {0.73, 0.26, 0.82, 0.77, 0.51, 0.97, 0.66, 0.88, 0.57, 0.35, 0.24, 0.41, 0.89, 0.35, 0.38, 0.52, 0.20, 0.49, 0.79, 0.75},
        {0.53, 0.21, 0.33, 0.87, 0.97, 0.62, 0.61, 0.58, 0.48, 0.35, 0.49, 0.90, 0.48, 0.24, 0.27, 0.55, 0.76, 0.43, 0.47, 0.25},
        {0.72, 0.76, 0.62, 0.47, 0.81, 0.27, 0.56, 0.38, 0.97, 0.45, 0.54, 0.21, 0.49, 0.25, 0.38, 0.62, 0.40, 0.69, 0.79, 0.45},
        {0.83, 0.86, 0.12, 0.67, 0.51, 0.47, 0.56, 0.18, 0.67, 0.75, 0.34, 0.91, 0.19, 0.65, 0.38, 0.82, 0.60, 0.29, 0.69, 0.49},
        {0.43, 0.64, 0.22, 0.29, 0.38, 0.71, 0.82, 0.26, 0.77, 0.83, 0.89, 0.12, 0.68, 0.24, 0.17, 0.25, 0.36, 0.23, 0.47, 0.55},
        {0.92, 0.76, 0.72, 0.57, 0.41, 0.62, 0.21, 0.88, 0.78, 0.75, 0.79, 0.70, 0.58, 0.44, 0.17, 0.15, 0.76, 0.33, 0.27, 0.95},
        {0.02, 0.46, 0.12, 0.27, 0.71, 0.37, 0.36, 0.18, 0.37, 0.35, 0.74, 0.31, 0.19, 0.85, 0.58, 0.82, 0.80, 0.49, 0.79, 0.15},
        {0.53, 0.61, 0.93, 0.77, 0.57, 0.32, 0.71, 0.48, 0.28, 0.15, 0.59, 0.60, 0.78, 0.54, 0.87, 0.65, 0.26, 0.73, 0.47, 0.25},
        {0.22, 0.96, 0.42, 0.57, 0.21, 0.77, 0.86, 0.18, 0.47, 0.35, 0.24, 0.61, 0.19, 0.65, 0.38, 0.32, 0.60, 0.69, 0.29, 0.45}
        };

    // Print A matrix
    //printf("Matrix A:\n");
    //print_matrix(n, A);

    double est = _norm1est( (double*)A, n, 300);

    printf("est: %lf\n", est);

    //est = frobenius_norm_cblas(A, n);
    //printf("est Frob: %lf\n", est);

    return 0;
}
