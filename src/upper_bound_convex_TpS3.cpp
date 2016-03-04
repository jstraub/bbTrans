/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "bbTrans/upper_bound_convex_TpS3.h"

namespace bb {

UpperBoundConvexTpS3::UpperBoundConvexTpS3(BS3<NodeS3>& boundS3) 
  : boundS3_(boundS3)
{ }

double UpperBoundConvexTpS3::Evaluate(const NodeTpS3& node) {
  double ub = 1e99;
  for (uint32_t i=0; i<5; ++i) {
    double ub_i = boundS3_.evaluate(NodeTpS3.GetNodeS3(i));
    if (ub_i < ub)
      ub = ub_i;
  }
  return ub;
}

double UpperBoundConvexTpS3::EvaluateAndSet(NodeTpS3& node) {
  double ub = Evaluate(node);
  node.SetUB(ub);
  return ub;
}

}
