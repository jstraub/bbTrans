/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <vector>
#include <list>
#include "optRot/node.h"
#include "optRot/bound.h"

namespace OptRot {

class BranchAndBound {
 public:
  BranchAndBound(Bound& lower_bound, Bound& upper_bound);
  ~BranchAndBound() = default;
  Node Compute(std::list<Node>& nodes);
 private:
  Bound& lower_bound_;
  Bound& upper_bound_;
};

}
