/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <Eigen/Dense>
#include "optRot/numeric_helpers.h"

namespace OptRot {

/// vMF distribution templated on the dimension.
template <uint32_t D>
class vMF {
 public:
  vMF(const Eigen::Matrix<double, D, 1>& mu, double tau, double pi);
  vMF(const vMF<D>& vmf) = default;
  ~vMF() = default;
  double GetPi() const {return pi_;}
  double GetTau() const {return tau_;}
  const Eigen::Matrix<double, D, 1>& GetMu() const {return mu_;}
  double GetLogZ() const;
 private:
  Eigen::Matrix<double, D, 1> mu_;
  double tau_;
  double pi_;
};

template<uint32_t D>
double ComputeLogvMFtovMFcost(const vMF<D>& vmf_A, const vMF<D>& vmF_B, 
  const Eigen::Matrix<double, D, 1>& mu_B_prime);


}
#include "optRot/vmf_impl.h"