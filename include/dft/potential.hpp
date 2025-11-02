#pragma once
#include <cmath>
#include <vector>
#include <boost/math/special_functions/cbrt.hpp>
#include "grid.hpp"

namespace dft { namespace pot {

struct HarmonicOscillator {
  double omega;
  explicit HarmonicOscillator(double omega_) : omega(omega_) {}
  double operator()(double x,double y,double z) const {
    return 0.5 * omega * omega * (x*x + y*y + z*z);
  }
};

inline double lda_exchange(double n) {
  if(n <= 0.0) return 0.0;
  const double Cx = std::pow(3.0/M_PI,1.0/3.0);
  return -Cx * boost::math::cbrt(n);
}

template<typename ExtV>
std::vector<double> build_local_potential(const Grid3D& g,
                                          const std::vector<double>& density,
                                          const ExtV& vext)
{
  std::vector<double> V(g.size());
  double x0 = 0.5 * g.L[0];
  double y0 = 0.5 * g.L[1];
  double z0 = 0.5 * g.L[2];
  for(size_t i=0;i<g.n[0];++i)
  for(size_t j=0;j<g.n[1];++j)
  for(size_t k=0;k<g.n[2];++k){
    size_t idx = g.index(i,j,k);
    double x = i*g.h[0] - x0;
    double y = j*g.h[1] - y0;
    double z = k*g.h[2] - z0;
    V[idx] = vext(x,y,z) + lda_exchange(density[idx]);
  }
  return V;
}

}} // namespace dft::pot
