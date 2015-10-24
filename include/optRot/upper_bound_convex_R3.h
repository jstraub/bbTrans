/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "optRot/node.h"
#include "optRot/numeric_helpers.h"
#include "optRot/normal.h"
#include "optRot/bound.h"
#include "optRot/lower_bound_R3.h"
#include "optRot/upper_bound_indep_R3.h"

namespace OptRot {

class UpperBoundConvexR3 : public Bound<NodeR3> {
 public:
  UpperBoundConvexR3(const std::vector<Normal<3>>& gmm_A, const
      std::vector<Normal<3>>& gmm_B, const Eigen::Quaterniond& q);
  virtual ~UpperBoundConvexR3() = default;
  virtual double Evaluate(const NodeR3& node);
 private:
  std::vector<Normal<3>> gmmT_;
};

Eigen::Vector3d FindMaxTranslationInNode(const Eigen::Matrix3d& A, 
    const Eigen::Vector3d& b, const NodeR3& node);

}

