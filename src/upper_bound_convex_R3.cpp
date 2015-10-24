/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "optRot/upper_bound_convex_R3.h"

namespace OptRot {

UpperBoundConvexR3::UpperBoundConvexR3(const
    std::vector<Normal<3>>& gmmA, const std::vector<Normal<3>>& gmmB, 
    const Eigen::Quaterniond& q) {
  ComputeGmmT(gmmA, gmmB, gmmT_, q);
}

double UpperBoundConvexR3::Evaluate(const NodeR3& node) {
  Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
  Eigen::Vector3d b = Eigen::Vector3d::Zero();
  double c = 0.;
  for (auto& gT : gmmT_) {
    Eigen::Vector3d tU = FindMinTranslationInNode(gT.GetOmega(),
        gT.GetXi(), node);
    Eigen::Vector3d tL = FindMaxTranslationInNode(gT.GetOmega(),
        gT.GetXi(), node);
    double L = -0.5*(tL-gT.GetMu()).transpose() *
      gT.GetSigmaLDLT().solve(tL-gT.GetMu());
    double U = -0.5*(tU-gT.GetMu()).transpose() *
      gT.GetSigmaLDLT().solve(tU-gT.GetMu());
    double g = (1.-exp(L-U))*exp(U)/(U-L);
    double h = (U*exp(L-U)-L)*exp(U)/(U-L);
    double D = gT.GetPi() / (sqrt((2.*M_PI)*(2.*M_PI)*(2.*M_PI)) *
        exp(gT.GetLogDetSigma()));
    A -= 0.5*D*g*gT.GetOmega();
    b += D*g*gT.GetXi();
    c += D*(h-0.5*g*(gT.GetMu().transpose()*gT.GetXi())(0));
  }
  Eigen::Vector3d t = FindMinTranslationInNode(-A, 0.5*b, node);
  double ub = (t.transpose()*A*t)(0) + (b.transpose()*t)(0) + c;
  return ub;
}

Eigen::Vector3d FindMaxTranslationInNode(const Eigen::Matrix3d& A, 
    const Eigen::Vector3d& b, const NodeR3& node) {
  // Check corners of box.
  Eigen::VectorXd vals(8);
  for (uint32_t i=0; i<8; ++i) {
    Eigen::Vector3d t;
    node.GetBox().GetCorner(i, t);
    vals(i) = (t.transpose()*A*t - 2.*t.transpose()*b)(0);
  }
  uint32_t id_max = 0;
  vals.maxCoeff(&id_max);
  Eigen::Vector3d t;
  node.GetBox().GetCorner(id_max, t);
  return t;
}

}
