/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "bbTrans/lower_bound_TpS3.h"

namespace bb {

LowerBoundTpS3::LowerBoundTpS3(BS3<NodeS3>& boundS3) 
  : boundS3_(boundS3)
{ }

double LowerBoundTpS3::Evaluate(const NodeTpS3& node) {
  double lb = -1e99;
  for (uint32_t i=0; i<5; ++i) {
    double lb_i = boundS3_.evaluate(NodeTpS3.GetNodeS3(i));
    if (lb_i > lb)
      lb = lb_i;
  }
  return lb;
}

double LowerBoundTpS3::EvaluateAndSet(NodeTpS3& node) {
  double lb = -1e99;
  uint32_t id_max = 0;
  for (uint32_t i=0; i<5; ++i) {
    double lb_i = boundS3_.Evaluate(NodeTpS3.GetNodeS3(i));
    if (lb_i > lb) {
      lb = lb_i;
      id_max = i;
    }
  }
//  Eigen::Matrix<double,3,9> xs;
//  Eigen::Matrix<double,9,1> lbs;
//  Evaluate(node, xs, lbs);
//  uint32_t id_max = 0;
////  double lb = lbs.maxCoeff(&id_max);
//  double lb = lbs(id_max);
//  node.SetLbArgument(xs.col(id_max));
  node.SetLB(lb);
  // Set the LB argument in the S3 node
  boundS3_.EvaluateAndSet(NodeTpS3.GetNodeS3(id_max));
  // Copy the LB argument over to the TpS3 node
  node.SetLbArgument(NodeTpS3.GetNodeS3(id_max).GetLbArgument());
  return lb;
}

}
