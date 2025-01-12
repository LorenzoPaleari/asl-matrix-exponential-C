#pragma once

#include <string>

/* Prototype of the algorithm's function */
typedef int(*comp_func)(double *A, double **X, int n);
typedef int(*comp_func_no_blas)(double *A, double *A2, double *A4, double *A6, double *A8, double *A10, double **X, int n);

void add_function(comp_func f, std::string name);
void add_function_no_blas(comp_func_no_blas f, std::string name);
void register_functions();
void register_functions_no_blas();

// Imported from C
// extern "C" ...
extern "C" int mexp_model(double *A, double **X, int n);
extern "C" int mexp(double *A, double **X, int n);
extern "C" int mexp_no_blas(double *A, double *A2, double *A4, double *A6, double *A8, double *A10, double **X, int n);