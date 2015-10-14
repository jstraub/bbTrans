/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include "optRot/upper_bound_convexity_log.h"

namespace OptRot {

UpperBoundConvexityLog::UpperBoundConvexityLog(const vMFMM<3>&
    vmf_mm_A, const vMFMM<3>& vmf_mm_B) 
  : vmf_mm_A_(vmf_mm_A), vmf_mm_B_(vmf_mm_B)
{}

double UpperBoundConvexityLog::Evaluate(const Node& node) {
  return 0.;
}

}
