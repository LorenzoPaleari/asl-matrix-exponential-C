# High-Performance Matrix Exponential in C

This repository contains a highly optimized C implementation of the matrix exponential algorithm introduced by Awad H. Al-Mohy and Nicholas J. Higham. Designed for **high performance**, it employs advanced optimizations such as blocking, vectorization, and LU decomposition to deliver results competitive with state-of-the-art tools like SciPy.

### Authors
- [Alejandro Cuadrón Lafuente](https://github.com/AlexCuadron)
- [Antonino Orofino](https://github.com/antoorofino)
- [Lorenzo Paleari](https://github.com/LorenzoPaleari)
- [Julian Sainz Martinez](https://github.com/jellothere)

## Project Overview

Developed as part of the **Advanced System Lab** course at ETH Zürich (2023), this project tackles computational challenges in matrix exponential calculation using advanced optimizations. For detailed results, see our [`report`](./docs/Report.pdf).

### Key Features

- **Algorithm**: Al-Mohy & Higham's scaling and squaring method.
- **Optimizations**:
  - Strength reduction and precomputation.
  - Instruction-level parallelism (ILP) and loop unrolling.
  - AVX2-based vectorization for key operations.
  - Blocking techniques for efficient memory utilization.
  - Transition from Gaussian elimination to LU decomposition for complex system solving.

### Motivation

Efficient matrix exponential computation is crucial in various fields such as:

- **Differential equations** (linear and partial).
- **Quantum mechanics**.
- **Control theory**.
- **Network analysis**.

Existing implementations often suffer from overscaling issues or are not optimized for high-performance systems. This project provides a fully optimized C implementation to address these gaps.

## Installation

To compile and run the project, ensure the following dependencies are available:

- **C Compiler**: GCC 10.2.1 or later.
- **BLAS Library**: Intel Math Kernel Library (MKL) or equivalent.
- **Intel Tools**: For roofline analysis and performance evaluation.

Ensure that your environment is properly configured for AVX2 and MKL optimizations to achieve peak performance.

## Results

This implementation achieves:

- **Runtime reduction**: From ~5s to **850ms** for $1024 \times 1024$ matrices.
- **Peak performance**: 10 flops/cycle, nearing BLAS's `dgemm` peak of 14 flops/cycle.
- **Competitiveness**: Matches SciPy's `linalg.expm` up to $512 \times 512$ matrices.
