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

template<class NodeLin>
class UpperBoundIndepLin : public Bound<NodeLin> {
 public:
  UpperBoundIndepLin(UpperBoundIndepS3& boundS3);
  virtual ~UpperBoundIndepLin() = default;
  virtual double Evaluate(const NodeLin& node);
  virtual double EvaluateAndSet(NodeLin& node);
 private:
  UpperBoundIndepS3& boundS3_;
};

}

