/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "optRot/node_S3.h"
#include "optRot/numeric_helpers.h"
#include "optRot/vmf.h"
#include "optRot/vmf_mm.h"
#include "optRot/bound.h"
#include "optRot/upper_bound_indep_S3.h"

namespace OptRot {

class UpperBoundConvexS3 : public Bound<NodeS3> {
 public:
  UpperBoundConvexS3(const vMFMM<3>& vmf_mm_A, const vMFMM<3>& vmf_mm_B);
  virtual double Evaluate(const NodeS3& node);
 private:
  const vMFMM<3>& vmf_mm_A_;
  const vMFMM<3>& vmf_mm_B_;
};

Eigen::Matrix<double,4,4> BuildM(const Eigen::Vector3d& u, const
    Eigen::Vector3d& v);

template<uint32_t D>
bool FindLambda(const Eigen::Matrix<double, D,D>& A, const
    Eigen::Matrix<double, D,D>& B, double* lambda) {
  Eigen::GeneralizedEigenSolver<Eigen::Matrix<double,D,D>> ges(A, B, true);
  // eigenvalues are alphas/betas
  if ((ges.betas().array().abs() > 1e-6).all()) {
//    std::cout << "FindLambda: non singular EVs" << std::endl;
    uint32_t id_max = 0;
    Eigen::Matrix<double, D, 1> ev = ges.eigenvalues().real();
    double ev_max = ev.maxCoeff(&id_max);
//    Eigen::Matrix<double,D,1> alpha = ges.eigenvectors().col(id_max);
    Eigen::ColPivHouseholderQR<Eigen::Matrix<double,D,D>> qr(A-ev_max*B);
    if (qr.rank() < D) {
//      std::cout << "FindLambda: cannot find eigen vector rank " << qr.rank() << " < " << D << std::endl;
      return false;
    }
//    std::cout << "FindLambda: can find eigen vector." << std::endl;
    Eigen::Matrix<double,D,1> alpha = qr.solve(Eigen::Matrix<double,D,1>::Zero());
//    std::cout << "FindLambda: alphas = " << alpha.transpose() << std::endl;
    if ((alpha.array() >= 0.).all() || (alpha.array() <= 0.).all()) {
//      std::cout << "FindLambda: lambda = " << ev_max << std::endl;
      *lambda = ev_max;
      return true; 
    }
  }
  return false;
}

}
