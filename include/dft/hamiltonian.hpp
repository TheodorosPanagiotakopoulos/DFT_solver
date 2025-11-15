#pragma once
#include <Eigen/Sparse>
#include <vector>
#include "grid.hpp"

namespace dft {

using SpMat = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

struct HamiltonianFD {
  SpMat H;
  void assemble(const Grid3D& g, const std::vector<double>& V){
    size_t N = g.size();
    std::vector<Triplet> triplets;
    triplets.reserve(7*N);
    for(size_t i=0;i<g.n[0];++i){
      size_t ip=g.pbc(i+1,0), im=g.pbc(i-1,0);
      for(size_t j=0;j<g.n[1];++j){
        size_t jp=g.pbc(j+1,1), jm=g.pbc(j-1,1);
        for(size_t k=0;k<g.n[2];++k){
          size_t kp=g.pbc(k+1,2), km=g.pbc(k-1,2);
          size_t idx=g.index(i,j,k);
          double hx=g.h[0], hy=g.h[1], hz=g.h[2];
          double lap=-(2.0/(hx*hx)+2.0/(hy*hy)+2.0/(hz*hz));
          triplets.emplace_back(idx, idx, -0.5*lap + V[idx]);
          triplets.emplace_back(idx, g.index(ip,j,k), -0.5/(hx*hx));
          triplets.emplace_back(idx, g.index(im,j,k), -0.5/(hx*hx));
          triplets.emplace_back(idx, g.index(i,jp,k), -0.5/(hy*hy));
          triplets.emplace_back(idx, g.index(i,jm,k), -0.5/(hy*hy));
          triplets.emplace_back(idx, g.index(i,j,kp), -0.5/(hz*hz));
          triplets.emplace_back(idx, g.index(i,j,km), -0.5/(hz*hz));
        }
      }
    }
    H.resize(N,N);
    H.setFromTriplets(triplets.begin(),triplets.end());
    H.makeCompressed();
  }
};

}
