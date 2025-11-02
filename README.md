# dft-fd — Finite-Difference Kohn–Sham DFT (C++)

A compact but serious Kohn–Sham DFT prototype built the *engineering* way:

- **Discretization:** uniform finite differences on a 3D grid (2nd-order Laplacian, periodic BCs).
- **Linear algebra:** **Eigen** sparse matrices.
- **Eigen-solver:** **ARPACK++** iterative solver (`DFT_USE_ARPACKPP`), with a small **Eigen** dense fallback for tiny demos.
- **SCF:** density mixing (linear), LDA-exchange (Dirac) as the XC demo.
- **Performance & portability:** **Kokkos** (optional) for portable density builds, **TBB** for high-level parallel loops, C++20 threads/`std::async` for background I/O.
- **Correctness ergonomics:** **Boost.Units** at API boundaries; **Boost.Math** for stable special functions (e.g. `cbrt`).

This repo is intentionally minimal yet professional: clean layout, clear extension points (Poisson/Hartree, better XC, Pulay/Broyden mixing, pseudopotentials, higher-order stencils), and CI-friendly CMake.

---

## Build

**Dependencies**
- C++20 compiler
- [Eigen 3.4+](https://eigen.tuxfamily.org)
- [oneTBB](https://github.com/oneapi-src/oneTBB)
- [Boost](https://www.boost.org/) (headers; uses Units and Math)
- Optional: [Kokkos](https://github.com/kokkos/kokkos) (set `-DUSE_KOKKOS=ON`)
- Optional: [ARPACK++](https://github.com/opencollab/arpackpp) + ARPACK/BLAS/LAPACK (set `-DUSE_ARPACK=ON`)

```bash
git clone https://github.com/yourname/dft-fd.git
cd dft-fd
cmake -B build -S . -DUSE_ARPACK=ON -DUSE_KOKKOS=ON \
      -DARPACKPP_INCLUDE_DIR=/path/to/arpackpp/include
cmake --build build -j
