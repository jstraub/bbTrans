/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "bbTrans/node_S3.h"
#include "bbTrans/node_TpS3.h"
#include "bbTrans/numeric_helpers.h"
#include "bbTrans/bound.h"
#include "bbTrans/lower_bound_S3.h"

namespace bb {

class LowerBoundTpS3 : public Bound<NodeTpS3> {
 public:
  LowerBoundTpS3(LowerBoundS3& boundS3);
  virtual ~LowerBoundTpS3() = default;
  // find the highest lower bound over the 5 S3 nodes managed by
  // NodeTpS3
  virtual double Evaluate(const NodeTpS3& node);
  virtual double EvaluateAndSet(NodeTpS3& node);
 private:
  LowerBoundS3& boundS3_;
};

}
