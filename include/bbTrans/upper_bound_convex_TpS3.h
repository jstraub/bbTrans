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
#include "bbTrans/upper_bound_convex_S3.h"

namespace bb {

class UpperBoundConvexTpS3 : public Bound<NodeTpS3> {
 public:
  UpperBoundConvexTpS3(UpperBoundConvexS3& boundS3);
  virtual ~UpperBoundConvexTpS3() = default;
  virtual double Evaluate(const NodeTpS3& node);
  virtual double EvaluateAndSet(NodeTpS3& node);
 private:
  UpperBoundConvexS3& boundS3_;
};

}

