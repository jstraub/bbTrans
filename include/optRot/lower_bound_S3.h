/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "optRot/node_S3.h"
#include "optRot/numeric_helpers.h"
#include "optRot/vmf.h"
#include "optRot/vmf_mm.h"
#include "optRot/bound.h"

namespace OptRot {

class LowerBoundS3 : public Bound<NodeS3> {
 public:
  LowerBoundS3(const vMFMM<3>& vmf_mm_A, const vMFMM<3>& vmf_mm_B);
  virtual ~LowerBoundS3() = default;
  virtual double Evaluate(const NodeS3& node);
  virtual double EvaluateAndSet(NodeS3& node);
 private:
  void Evaluate(const NodeS3& node, std::vector<Eigen::Quaterniond>& qs,
      Eigen::Matrix<double,5,1>& lbs);
  const vMFMM<3>& vmf_mm_A_;
  const vMFMM<3>& vmf_mm_B_;
};

}
