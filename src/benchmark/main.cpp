/**
*      _________   _____________________  ____  ______
*     / ____/   | / ___/_  __/ ____/ __ \/ __ \/ ____/
*    / /_  / /| | \__ \ / / / /   / / / / / / / __/
*   / __/ / ___ |___/ // / / /___/ /_/ / /_/ / /___
*  /_/   /_/  |_/____//_/  \____/\____/_____/_____/
*
*  http://www.acl.inf.ethz.ch/teaching/fastcode
*  How to Write Fast Numerical Code 263-2300 - ETH Zurich
*  Copyright (C) 2019 
*                   Tyler Smith        (smitht@inf.ethz.ch) 
*                   Alen Stojanov      (astojanov@inf.ethz.ch)
*                   Gagandeep Singh    (gsingh@inf.ethz.ch)
*                   Markus Pueschel    (pueschel@inf.ethz.ch)
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program. If not, see http://www.gnu.org/licenses/.
*/

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib> 

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "tsc_x86.h"
#include "utils.h"
#include "common.h"
#include "mkl.h"

using namespace std;

#define NR 32
#define CYCLES_REQUIRED 1e8
#define REP 1
#define EPS (1e-3)
#define EPS_TEST (1e-12)

/* Global vars, used to keep track of testing functions */
vector<comp_func> userFuncs;
vector<comp_func_no_blas> userFuncs_no_blas;
vector<string> funcNames;
int numFuncs = 0;

/*
* Calls the python script to compute the matrix X_base
*/
int get_python_computation(double *A, double *X_base, int n, string cmd, char *tmp_file) {
  cout << "> Writing matrix to file ... " << flush;
  ofstream out(tmp_file, ios::binary);
  out.write(reinterpret_cast<const char*>(A), n * n * sizeof(double));
  out.close();
  cout << "✅" << endl;
  
  cout << "> Executing python script ... " << flush;
  int ret = system(cmd.c_str());
  if (ret != 0) {
    cout << "❌ Error calling Python script!\n" << endl;
    return 1;
  } else {
    cout << "✅" << endl;
  }

  cout << "> Reading matrix from file ... " << flush;
  ifstream in(tmp_file, ios::binary);
  in.read(reinterpret_cast<char*>(X_base), n * n * sizeof(double));
  in.close();
  cout << "✅" << endl;
  return 0;
}

/*
* Registers a function to be tested by the driver program.
* Registers a string description of the function as well.
*/
void add_function(comp_func f, string name) {
    userFuncs.push_back(f);
    funcNames.push_back(name);
    numFuncs++;
}
void add_function_no_blas(comp_func_no_blas f, string name) {
    userFuncs_no_blas.push_back(f);
    funcNames.push_back(name);
    numFuncs++;
}

/*
* Computes the number of cycles required per function
*/
double perf_test(comp_func f, double* A, double** X, int n) {
    double cycles = 0.;
    size_t num_runs = 1;
    // double multiplier = 1;
    myInt64 start, end;
    
    if (n < 2000)
      num_runs = 1;
    if (n < 1000)
      num_runs = 10;
    if (n < 500)
      num_runs = 100;
    if (n < 200)
      num_runs = 250;
    if (n < 100)
      num_runs = 500;
    if (n < 50)
      num_runs = 1000;
    
    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    // do {
    //     num_runs = num_runs * multiplier;
    //     start = start_tsc();
    //     for (size_t i = 0; i < num_runs; i++)
    //         f(A, X, n);
    //     end = stop_tsc(start);

    //     cycles = (double)end;
    //     multiplier = (CYCLES_REQUIRED) / (cycles);
        
    // } while (multiplier > 2);

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i)
            f(A, X, n);
        end = stop_tsc(start);
        cout << "." << flush;

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;
    }
    total_cycles /= REP;
    cycles = total_cycles;
    return  cycles;
}

double perf_test_no_blas(comp_func_no_blas f, double* A, double *A2, double *A4, double *A6, double *A8, double *A10, double **X, int n) {
    double cycles = 0.;
    size_t num_runs = 1;
    // double multiplier = 1;
    myInt64 start, end;
    
    if (n < 2000)
      num_runs = 10;
    if (n < 1000)
      num_runs = 50;
    if (n < 500)
      num_runs = 100;
    if (n < 200)
      num_runs = 250;
    if (n < 100)
      num_runs = 500;
    if (n < 50)
      num_runs = 1000;
    
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i)
            f(A, A2, A4, A6, A8, A10, X, n);
        end = stop_tsc(start);
        cout << "." << flush;

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;
    }
    total_cycles /= REP;
    cycles = total_cycles;
    return  cycles;
}

int benchmark_standard(int argc, char* argv[]){
  cout << endl;
  if (argc != 6) {
    cout << "Usage: " << argv[0] << " <N> <pyhton_path> <pyhton_script> <tmp_file> <test_dir>" << endl;
    return 1;
  }
  double perf;
  int i, return_value, fails;
  fails = 0;
  int n = atoi(argv[1]);
  double *A, *X, *X_base, *X_test;
  X = NULL;
  X_test = NULL;
  // allocates memory for matrix X_base and A
  build((double**) &X_base, n);
  build((double**) &A, n);
  ofstream outFile;

  cout << "> Cheking output directory ... " << flush;
  string cmd = "ls " + string(argv[5])  + " >/dev/null 2>&1";
  int ret = system(cmd.c_str());
  if (ret != 0) {
    cout << "creating ... " << flush;
    cmd = "mkdir " + string(argv[5])  + " >/dev/null 2>&1";
    system(cmd.c_str());
  }
  cout << "✅" << endl;

  cout << "> Registering functions ... " << flush;
  register_functions();
  if (numFuncs == 0){
    cout << endl;
    cout << "❌ No functions registered - nothing for driver to do." << endl;
    return 1;
  }
  cout << "✅" << endl;
  cout << "  "<< numFuncs << " functions registered" << endl;
  
  // Opens the output file
  #if !FLOPS && !ROOFLINE
  outFile.open(string(argv[5]) + "/" + to_string(n) + ".csv");
  outFile << "Function_name, Return_value, Error, Cycles" << endl;
  #endif

  // Bounds used for random generation of matrix A
  int test_num = 6;
  double low_bound_list[test_num] = {-10, -30, -1, -0.1, -0.01, -0.00001};
  double up_bound_list[test_num] = {10, 30, 1, .1, .01, .00001};
  bool is_triangular_list[test_num] = {false, true, false, false, false, false};
for (int t = 0; t < test_num; t++){
    cout << "\n  Performing test " << t + 1 << "/" << test_num << "\n" << endl;
    // Creates a test matrix given the bounds
    rands((double*) A, n, low_bound_list[t], up_bound_list[t], is_triangular_list[t]);

    #ifndef ROOFLINE
    //Calls the python script to compute the matrix X_base
    cmd = string(argv[2]) + " " + string(argv[3]) + " " + string(argv[4]);
    if(get_python_computation(A, X_base, n, cmd, argv[4]))
      return 1;
    
    // Computes the matrix X_test using the model implementation
    mexp_model(A, &X_test, n);
    #endif

    for (i = 0; i < numFuncs; i++) {
      comp_func f = userFuncs[i];

      // Checks the given function for validity
      cout << "\n> Checks '" << funcNames[i] << "' function for validity ... " << flush;
      return_value = f(A, &X, n);

      #ifndef ROOFLINE
      double error_rel = rel_err((double*) X_base, (double*) X, n);
      double error_test = rel_err((double*) X_test, (double*) X, n);
      if (error_rel > EPS || error_test > EPS_TEST) {
        cout << "❌" << endl;
        fails++;
      } else {
        cout << "✅" << endl;
      }
      cout << "  Return value: " << return_value << " - Error: " << error_rel << " - Error test base: " << error_test << endl;
      #endif

      #if !FLOPS && !ROOFLINE
      outFile << funcNames[i] << ", " << return_value << ", " << error_rel << ", ";
      // Computes the number of cycles required
      cout << "> Performance calculation ." << flush;
      perf = perf_test(f, A, &X, n);
      cout << " ✅\n  Cycles: " << perf << "\n" << endl;
      outFile << perf << endl;
      #endif
    }
  }
  
  #if !FLOPS && !ROOFLINE
  outFile.close();
  #endif

  // Deallocates memory
  destroy((double*) A);
  destroy((double*) X);
  destroy((double*) X_base);
  destroy((double*) X_test);
  return fails;
}

int benchmark_no_blas(int argc, char* argv[]){
  cout << endl;
  if (argc != 6) {
    cout << "Usage: " << argv[0] << " <N> <pyhton_path> <pyhton_script> <tmp_file> <test_dir>" << endl;
    return 1;
  }
  double perf;
  int i, return_value, fails;
  int n = atoi(argv[1]);
  double *A, *A2, *A4, *A6, *A8, *A10, *X;
  X = NULL;
  // allocates memory for matrix X_base and A
  build((double**) &A, n);
  build((double**) &A2, n);
  build((double**) &A4, n);
  build((double**) &A6, n);
  build((double**) &A8, n);
  build((double**) &A10, n);
  build((double**) &X, n);
  ofstream outFile;

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A, n, A, n, 0.0, A2, n);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A2, n, A2, n, 0.0, A4, n);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A4, n, A2, n, 0.0, A6, n);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A4, n, A4, n, 0.0, A8, n);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A6, n, A4, n, 0.0, A10, n);

  cout << "> Cheking output directory ... " << flush;
  string cmd = "ls " + string(argv[5])  + " >/dev/null 2>&1";
  int ret = system(cmd.c_str());
  if (ret != 0) {
    cout << "creating ... " << flush;
    cmd = "mkdir " + string(argv[5])  + " >/dev/null 2>&1";
    system(cmd.c_str());
  }
  cout << "✅" << endl;

  cout << "> Registering functions ... " << flush;
  register_functions_no_blas();
  if (numFuncs == 0){
    cout << endl;
    cout << "❌ No functions registered - nothing for driver to do." << endl;
    return 1;
  }
  cout << "✅" << endl;
  cout << "  "<< numFuncs << " functions registered" << endl;
  
  // Opens the output file
  #if !FLOPS && !ROOFLINE
  outFile.open(string(argv[5]) + "/" + to_string(n) + ".csv");
  outFile << "Function_name, Return_value, Cycles" << endl;
  #endif

  // Bounds used for random generation of matrix A
  int test_num = 6;
  double low_bound_list[test_num] = {-10, -30, -1, -0.1, -0.01, -0.00001};
  double up_bound_list[test_num] = {10, 30, 1, .1, .01, .00001};
  bool is_triangular_list[test_num] = {false, true, false, false, false, false};
for (int t = 0; t < test_num; t++){
    cout << "\n  Performing test " << t + 1 << "/" << test_num << "\n" << endl;
    // Creates a test matrix given the bounds
    rands((double*) A, n, low_bound_list[t], up_bound_list[t], is_triangular_list[t]);
    rands((double*) X, n, low_bound_list[t], up_bound_list[t], is_triangular_list[t]);

    for (i = 0; i < numFuncs; i++) {
      comp_func_no_blas f = userFuncs_no_blas[i];

      // Checks the given function for validity
      cout << "\n> Checks '" << funcNames[i] << "' function for validity ... " << flush;
      return_value = f(A, A2, A4, A6, A8, A10, &X, n);

      #if !FLOPS && !ROOFLINE
      outFile << funcNames[i] << ", " << return_value << ", ";
      // Computes the number of cycles required
      cout << "> Performance calculation ." << flush;
      perf = perf_test_no_blas(f, A, A2, A4, A6, A8, A10, &X, n);
      cout << " ✅\n  Cycles: " << perf << "\n" << endl;
      outFile << perf << endl;
      #endif
    }
  }
  
  #if !FLOPS && !ROOFLINE
  outFile.close();
  #endif

  // Deallocates memory
  destroy((double*) A);
  destroy((double*) A2);
  destroy((double*) A4);
  destroy((double*) A6);
  destroy((double*) A8);
  destroy((double*) A10);
  destroy((double*) X);
  return fails;
}

int main(int argc, char* argv[]) {
  #ifdef NO_BLAS
    return benchmark_no_blas(argc, argv);
  #else
    return benchmark_standard(argc, argv);
  #endif
}