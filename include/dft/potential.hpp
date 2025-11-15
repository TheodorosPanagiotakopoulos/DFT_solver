#pragma once
#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include "grid.hpp"

namespace dft { namespace pot {

struct HarmonicPlusCoulomb {
  double omega;
  std::vector<std::array<double,3>> centers;
  std::vector<double> charges;

  HarmonicPlusCoulomb(double omega_,
                      const std::vector<std::array<double,3>>& centers_,
                      const std::vector<double>& charges_)
      : omega(omega_), centers(centers_), charges(charges_) {}

  double operator()(double x, double y, double z) const {
    double V = 0.5 * omega * omega * (x*x + y*y + z*z);
    for (size_t a = 0; a < centers.size(); ++a) {
      double dx = x - centers[a][0];
      double dy = y - centers[a][1];
      double dz = z - centers[a][2];
      double r = std::sqrt(dx*dx + dy*dy + dz*dz + 1e-6);
      V += -charges[a] / r;
    }
    return V;
  }
};

inline double lda_exchange(double n) {
  if (n <= 1e-12) return 0.0;
  const double Cx = std::pow(3.0 / M_PI, 1.0 / 3.0);
  return -Cx * std::pow(n, 1.0 / 3.0);
}

inline double lda_correlation(double n) {
  if (n <= 1e-12) return 0.0;
  double rs = std::pow(3.0 / (4.0 * M_PI * n), 1.0 / 3.0);
  double A = 0.0311, B = -0.048, C = 0.0020, D = -0.0116;
  double gamma = -0.1423, beta1 = 1.0529, beta2 = 0.3334;
  double eps_c;
  if (rs < 1.0)
    eps_c = A * std::log(rs) + B + C * rs * std::log(rs) + D * rs;
  else
    eps_c = gamma / (1.0 + beta1 * std::sqrt(rs) + beta2 * rs);
  return eps_c;
}

inline std::vector<double> solve_poisson_jacobi(
    const Grid3D& g,
    const std::vector<double>& density,
    int max_iters = 5000,
    double tol = 1e-6)
{
  size_t N = g.size();
  std::vector<double> V(N, 0.0), Vnew(N, 0.0);
  double hx2 = g.h[0]*g.h[0];
  double hy2 = g.h[1]*g.h[1];
  double hz2 = g.h[2]*g.h[2];
  double denom = 2.0 * (1.0/hx2 + 1.0/hy2 + 1.0/hz2);

  for (int it = 0; it < max_iters; ++it) {
    double max_diff = 0.0;
    for (size_t i = 0; i < g.n[0]; ++i)
    for (size_t j = 0; j < g.n[1]; ++j)
    for (size_t k = 0; k < g.n[2]; ++k) {
      size_t ip = g.pbc(i+1,0), im = g.pbc(i-1,0);
      size_t jp = g.pbc(j+1,1), jm = g.pbc(j-1,1);
      size_t kp = g.pbc(k+1,2), km = g.pbc(k-1,2);
      size_t idx = g.index(i,j,k);

      double lap =
        (V[g.index(ip,j,k)] + V[g.index(im,j,k)]) / hx2 +
        (V[g.index(i,jp,k)] + V[g.index(i,jm,k)]) / hy2 +
        (V[g.index(i,j,kp)] + V[g.index(i,j,km)]) / hz2;

      Vnew[idx] = (lap + 4.0 * M_PI * density[idx]) / denom;

      double diff = std::fabs(Vnew[idx] - V[idx]);
      if (diff > max_diff) max_diff = diff;
    }

    std::swap(V, Vnew);
    if (max_diff < tol) {
      std::cout << "Hartree solver converged in " << it + 1 << " iterations.\n";
      break;
    }
  }

  return V;
}

template<typename ExtV>
std::vector<double> build_local_potential(
    const Grid3D& g,
    const std::vector<double>& density,
    const ExtV& vext)
{
  std::vector<double> V(g.size(), 0.0);
  auto Vhartree = solve_poisson_jacobi(g, density);

  double x0 = 0.5 * g.L[0];
  double y0 = 0.5 * g.L[1];
  double z0 = 0.5 * g.L[2];

  for (size_t i = 0; i < g.n[0]; ++i)
  for (size_t j = 0; j < g.n[1]; ++j)
  for (size_t k = 0; k < g.n[2]; ++k) {
    size_t idx = g.index(i, j, k);
    double x = i * g.h[0] - x0;
    double y = j * g.h[1] - y0;
    double z = k * g.h[2] - z0;
    double n = density[idx];

    V[idx] = vext(x, y, z)
           + Vhartree[idx]
           + lda_exchange(n)
           + lda_correlation(n);
  }

  return V;
}

}}
