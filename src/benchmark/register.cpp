#include <immintrin.h>
#include "common.h"

/*
* Called by the driver to register the functions
*/
void register_functions() {
  //add_function(&mexp, "mexp");
  //add_function(&mexp_model, "mexp_model");
}

void register_functions_no_blas() {
  add_function_no_blas(&mexp_no_blas, "mexp_no_blas");
}