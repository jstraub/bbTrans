/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <vector>
#include "optRot/node.h"
#include "optRot/bound.h"

namespace OptRot {

class UpperBoundLog : public Bound {

};

class UpperBoundConvexityLog : public Bound {

};

class BranchAndBound {
 public:
  BranchAndBound();
  ~BranchAndBound() = default;
 private:
  
};

}
