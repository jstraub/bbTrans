/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <memory>
#include <vector>
#include <list>
#include <fstream>
#include "optRot/node.h"
#include "optRot/bound.h"

namespace OptRot {

template <class Node>
class BranchAndBound {
 public:
  BranchAndBound(Bound<Node>& lower_bound, Bound<Node>& upper_bound);
  ~BranchAndBound() = default;
  Node Compute(std::list<Node>& nodes, double eps, uint32_t max_it);
 private:
  Bound<Node>& lower_bound_;
  Bound<Node>& upper_bound_;
  uint32_t BoundAndPrune(std::list<Node>& nodes, double& lb, double&
      ub, double eps);
  void WriteStats(std::ofstream& out, std::list<Node>& nodes, double
      lb, double ub);
  void WriteNodes(std::ofstream& out, std::list<Node>& nodes, double
      lb, double ub);
};
}
#include "optRot/branch_and_bound_impl.h"
