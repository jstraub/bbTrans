/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "optRot/lower_bound_R3.h"

namespace OptRot {

LowerBoundR3::LowerBoundR3(const
    std::vector<Normal<3>>& gmmA, const std::vector<Normal<3>>& gmmB, 
    const Eigen::Quaterniond& q) {
  ComputeGmmT(gmmA, gmmB, gmmT_, q);
}

double LowerBoundR3::Evaluate(const NodeR3& node) {
  Eigen::Matrix<double,3,9> xs;
  Eigen::Matrix<double,9,1> lbs;
  Evaluate(node, xs, lbs);
  return lbs.maxCoeff();
}

double LowerBoundR3::EvaluateAndSet(NodeR3& node) {
  Eigen::Matrix<double,3,9> xs;
  Eigen::Matrix<double,9,1> lbs;
  Evaluate(node, xs, lbs);
  uint32_t id_max = 0;
  double lb = lbs.maxCoeff(&id_max);
  node.SetLB(lb);
  node.SetLbArgument(xs.col(id_max));
  return lb;
}

void LowerBoundR3::Evaluate(const NodeR3& node,
    Eigen::Matrix<double,3,9>& xs, Eigen::Matrix<double,9,1>& lbs) {
  xs.col(0) = node.GetBox().GetCenter();
  Eigen::Vector3d c;
  for (uint32_t i=0; i<8; ++i) {
    node.GetBox().GetCorner(i,c);
    xs.col(i+1) = c;
  }
  lbs = Eigen::VectorXd::Zero(9);
  for (uint32_t i=0; i<8; ++i) {
    for (auto& gT : gmmT_)
      lbs(i) += gT.GetPi() * gT.pdf(xs.col(i));
  }
}

void ComputeGmmT( const std::vector<Normal<3>>& gmmA, const
    std::vector<Normal<3>>& gmmB, std::vector<Normal<3>>& gmmT, const
    Eigen::Quaterniond& q) {
  gmmT.reserve(gmmA.size() * gmmB.size());
  Eigen::Matrix3d R = q.toRotationMatrix();
  for (auto& gA : gmmA) 
    for (auto& gB : gmmB) {
      gmmT.push_back(
          Normal<3>(gB.GetMu() - R*gA.GetMu(),
            R*gA.GetSigma()*R.transpose() + gB.GetSigma(), 
            gB.GetPi()*gA.GetPi()));
    }
}

}
