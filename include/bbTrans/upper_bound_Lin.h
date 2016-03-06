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

template<class UpperBound, class NodeLin>
class UpperBoundConvexLin : public Bound<NodeLin> {
 public:
  UpperBoundConvexLin(UpperBound& boundS3);
  virtual ~UpperBoundConvexLin() = default;
  virtual double Evaluate(const NodeLin& node);
  virtual double EvaluateAndSet(NodeLin& node);
 private:
  UpperBound& boundS3_;
};

}

#include "bbTrans/upper_bound_Lin_impl.hpp"
