#pragma once
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <stdexcept>

namespace dft {

struct Eigenpairs {
  Eigen::VectorXd eigenvalues;
  Eigen::MatrixXd eigenvectors;
};

inline Eigenpairs solve_lowest_dense(const Eigen::SparseMatrix<double>& H, int k){
  Eigen::MatrixXd Hd = Eigen::MatrixXd(H);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hd);
  if(es.info()!=Eigen::Success) throw std::runtime_error("eigensolve failed");
  Eigenpairs out;
  out.eigenvalues = es.eigenvalues().head(k);
  out.eigenvectors = es.eigenvectors().leftCols(k);
  return out;
}

}

