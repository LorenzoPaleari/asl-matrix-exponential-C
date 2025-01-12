#pragma once
#include <random>

template<typename T>
void rands(T * m, size_t n, double lo = .0, double up = 10.0, bool is_triangular = false)
{
    std::random_device rd;
    std::mt19937 gen{24};
    std::uniform_real_distribution<T> dist(lo / n * 8.0, up / n * 8.0);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            m[i*n + j] = (is_triangular && i > j) ? .0 : dist(gen);
}

template<typename T>
void build(T **a, int n)
{
    *a = static_cast<T *>(aligned_alloc(32, n * n * sizeof(T)));
}

template<typename T>
void destroy(T* m)
{
    free(m);
}

template<typename T>
T nrm_sqr_diff(T *x, T *y, int n) {
    T nrm_sqr = 0.0;
    for(int i = 0; i < n*n; i++) {
        nrm_sqr += (x[i] - y[i]) * (x[i] - y[i]);
    }
    
    if (isnan(nrm_sqr)) {
      nrm_sqr = INFINITY;
    }
    
    return nrm_sqr;
}

// Calculates the 1-norm of a matrix whihch is the maximum absolute column sum
template<typename T>
T norm_1(T *A, int n){
    T result = 0.0;
    T temp[n];
    
    for(int j = 0; j < n; j++){
        temp[j] = 0.0;
        for(int i = 0; i < n; i++)
            temp[j] += fabs(A[i*n + j]);
    }

    for(int i = 0; i < n; i++)
        result = fmax(temp[i], result);
    
    return result;
}

template<typename T>
T rel_err(T *x, T *y, int n) {
    T *tmp;
    build((T **) &tmp, n);
    T rel_err = 0.0;
    
    for(int i = 0; i < n*n; i++) {
        tmp[i] = x[i] - y[i];
    }
    
    rel_err = norm_1((T *) tmp, n) / norm_1((T *) x, n);
    
    destroy((T *) tmp);
    
    return rel_err;
}
