#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include "hamiltonian.hpp"
#include "eigensolver.hpp"
#include "potential.hpp"
#include "parallel.hpp"

namespace dft {

struct SCFOptions {
  int nelectrons = 2;
  int max_iters  = 50;
  double mixing  = 0.3;
  double energy_tol = 1e-6;
  double density_tol= 1e-5;
  int n_bands() const { return (nelectrons+1)/2; }
};

struct SCFResult {
  double Etot=0;
  int iters=0;
  std::vector<double> density;
  Eigen::VectorXd eigenvalues;
};

template<typename ExtV>
SCFResult run_scf(const Grid3D& g, const ExtV& vext, const SCFOptions& opt){
  size_t N = g.size();
  std::vector<double> n(N,1e-6), n_old(N,0.0);

  auto norm_diff=[&](const std::vector<double>& a,const std::vector<double>& b){
    double s=0.0;
    for(size_t i=0;i<N;++i){ double d=a[i]-b[i]; s+=d*d; }
    return std::sqrt(s);
  };

  double Etot_old=1e100;
  SCFResult out; out.density=n;

  for(int it=0; it<opt.max_iters; ++it){
    auto V = pot::build_local_potential(g,n,vext);
    HamiltonianFD H; H.assemble(g,V);
    int nb = opt.n_bands();
    auto ks = solve_lowest_dense(H.H, nb);

    std::vector<double> n_new(N,0.0);
    par::for_index(N,[&](size_t i){
      double acc=0.0;
      for(int b=0;b<nb;++b)
        acc += ks.eigenvectors(i,b)*ks.eigenvectors(i,b);
      n_new[i]=2.0*acc;
    });

    n_old=n;
    for(size_t i=0;i<N;++i)
      n[i]=(1.0-opt.mixing)*n[i]+opt.mixing*n_new[i];

    double e_kin=0.0; for(int b=0;b<nb;++b) e_kin+=ks.eigenvalues[b]*2.0;
    double e_pot=0.0; for(size_t i=0;i<N;++i) e_pot+=n[i]*V[i];
    double e_xc =0.0; for(size_t i=0;i<N;++i) e_xc +=n[i]*pot::lda_exchange(n[i]);
    double Etot=e_kin+e_pot+e_xc;

    double dn = norm_diff(n,n_old);
    std::cout<<"Iter "<<it+1<<"  E="<<Etot<<"  dE="<<(Etot-Etot_old)<<"  ||dn||="<<dn<<"\n";

    if(std::abs(Etot-Etot_old)<opt.energy_tol && dn<opt.density_tol){
      out.Etot=Etot; out.iters=it+1; out.density=n; out.eigenvalues=ks.eigenvalues;
      return out;
    }
    Etot_old=Etot;
  }
  out.Etot=Etot_old; out.density=n;
  return out;
}

}
