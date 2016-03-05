/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "bbTrans/node_TpS3.h"
#include "bbTrans/numeric_helpers.h"
#include "bbTrans/bound.h"
#include "bbTrans/lower_bound_TpS3.h"
#include "bbTrans/upper_bound_indep_S3.h"

namespace bb {

class UpperBoundIndepTpS3 : public Bound<NodeTpS3> {
 public:
  UpperBoundIndepTpS3(UpperBoundIndepS3& boundS3);
  virtual ~UpperBoundIndepTpS3() = default;
  virtual double Evaluate(const NodeTpS3& node);
  virtual double EvaluateAndSet(NodeTpS3& node);
 private:
  UpperBoundIndepS3& boundS3_;
};

}

