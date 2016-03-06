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

template<class NodeLin>
class LowerBoundLin : public Bound<NodeLin> {
 public:
  LowerBoundLin(LowerBoundS3& boundS3);
  virtual ~LowerBoundLin() = default;
  virtual double Evaluate(const NodeLin& node);
  virtual double EvaluateAndSet(NodeLin& node);
 private:
  LowerBoundS3& boundS3_;
};
typedef LowerBoundTpS3 LowerBoundLin<NodeTpS3>;
typedef LowerBoundAA   LowerBoundLin<NodeAA>;
}
#include "bbTrans/lower_bound_Lin_impl.hpp"
