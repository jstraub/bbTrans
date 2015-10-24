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
  // TODO might get a slightly tighter bound by taking the max value
  // over the corners of the box.
  double lb = 0.;
  Eigen::Vector3d x = node.GetBox().GetCenter();
  for (auto& gT : gmmT_)
    lb += gT.GetPi() * gT.pdf(x);
  return lb;
}

void ComputeGmmT( const std::vector<Normal<3>>& gmmA, const
    std::vector<Normal<3>>& gmmB, std::vector<Normal<3>>& gmmT, const
    Eigen::Quaterniond& q) {
  gmmT.reserve(gmmA.size() * gmmB.size());
  Eigen::Matrix3d R = q.toRotationMatrix();
  for (auto& gA : gmmA) 
    for (auto& gB : gmmB) {
      gmmT.push_back(
//          Normal<3>(R*gB.GetMu() - gA.GetMu(),
//            gA.GetSigma() + R*gB.GetSigma()*R.transpose(), 
//            gB.GetPi()*gA.GetPi()));
          Normal<3>(R.transpose()*gB.GetMu() - gA.GetMu(),
            gA.GetSigma() + R.transpose()*gB.GetSigma()*R, 
            gB.GetPi()*gA.GetPi()));
    }
}

}
