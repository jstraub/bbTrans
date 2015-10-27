/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "optRot/node_R3.h"
#include "optRot/numeric_helpers.h"
#include "optRot/normal.h"
#include "optRot/bound.h"

namespace OptRot {

class LowerBoundR3 : public Bound<NodeR3> {
 public:
  LowerBoundR3(const std::vector<Normal<3>>& gmmA, const
      std::vector<Normal<3>>& gmmB, const Eigen::Quaterniond& q);
  virtual ~LowerBoundR3() = default;
  virtual double Evaluate(const NodeR3& node);
 private:
  std::vector<Normal<3>> gmmT_;
};

void ComputeGmmT( const std::vector<Normal<3>>& gmmA, const
    std::vector<Normal<3>>& gmmB, std::vector<Normal<3>>& gmmT, const
    Eigen::Quaterniond& q);

}
