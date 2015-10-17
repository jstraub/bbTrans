/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include "optRot/numeric_helpers.h"

namespace OptRot {

double LogSumExp(const Eigen::VectorXd& x) {
  const double x_max = x.maxCoeff();
  return log((x.array() - x_max).exp().sum()) + x_max;
};

double SumExp(const Eigen::VectorXd& x) {
  const double x_max = x.maxCoeff();
  return (x.array() - x_max).exp().sum() * exp(x_max);
};

/// Compute log((exp(z) - exp(-z)) / z)
double ComputeLog2SinhOverZ(double z) {
  if (fabs(z) < 1.e-6) 
    return log(2.);
  else if (z < 50.) 
    return - log(z) + log(exp(2.*z) -1) - z;
  else
    return - log(z) + z;
};

}