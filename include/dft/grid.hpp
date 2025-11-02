#pragma once
#include <array>
#include <cstddef>
#include <vector>
#include <cassert>

namespace dft {

struct Grid3D {
  std::array<size_t,3> n;
  std::array<double,3> L;
  std::array<double,3> h;

  Grid3D(std::array<size_t,3> n_, std::array<double,3> L_)
      : n(n_), L(L_) {
    for(int d=0; d<3; ++d)
      h[d] = L[d] / static_cast<double>(n[d]);
  }

  size_t size() const { return n[0]*n[1]*n[2]; }

  size_t index(size_t i, size_t j, size_t k) const {
    return k + n[2]*(j + n[1]*i);
  }

  size_t pbc(int i, int dim) const {
    int N = static_cast<int>(n[dim]);
    int ii = i % N;
    if(ii < 0) ii += N;
    return static_cast<size_t>(ii);
  }
};

} // namespace dft
