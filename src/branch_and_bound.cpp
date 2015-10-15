/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include "optRot/branch_and_bound.h"

namespace OptRot {

BranchAndBound::BranchAndBound(Bound& lower_bound, Bound& upper_bound)
  : lower_bound_(lower_bound), upper_bound_(upper_bound)
{}

Node BranchAndBound::Compute(std::list<Node>& nodes) {
  double lb = -1.e6;
  double ub = 1.e6;
  double eps = 1.0e-3 * M_PI / 180.;
  uint32_t it = 0;
  uint32_t max_it = 10000;
  Node node_star = nodes.front();
  for (auto& node : nodes) {
    node.SetLB(lower_bound_.Evaluate(node));
    node.SetUB(upper_bound_.Evaluate(node));
  }

  uint32_t n_nodes = std::distance(nodes.begin(), nodes.end());
  while (it < max_it && node_star.GetBoundGap() > eps && n_nodes > 0) {
//  while (it < max_it && fabs(ub - lb) > eps && n_nodes > 0) {
    nodes.remove_if(IsPrunableNode(lb));
    auto node_i = std::max_element(nodes.begin(), nodes.end(), LessThanNodeUB());
    node_star = *std::max_element(nodes.begin(), nodes.end(), LessThanNodeLB());
    double lbn = node_i->GetLB();
    double ubn = node_i->GetUB();
//    std::cout << "@" << it << " # " << n_nodes << ": cur " 
//      << lbn << " < " << ubn << "\t global "
//      << lb << " < " << ub << " " << " |.| " << fabs(ub - lb) << "\t selected "
//      << node_star.GetLB() << " < " << node_star.GetUB() << " "
//      << node_star.GetLevel() << " " << node_star.GetIds()[0] 
//      << std::endl;
    if (lbn > lb) {
      // Remember this node as the current best.
      lb = lbn;
      ub = ubn;
    } else {
      // Branch and check resulting nodes.
      std::vector<Node> new_nodes = node_i->Branch();
      for (auto n : new_nodes) {
        n.SetUB(upper_bound_.Evaluate(n));
        if (ubn > lb) {
          // Remember this node since we cannot prune it.
          n.SetLB(lower_bound_.Evaluate(n));
          nodes.push_back(n);
          ++n_nodes;
        }
      }
      // Pop the examined node since we are branching here
      nodes.erase(node_i); 
      --n_nodes;
    }
    ++it;
  }
  std::cout << "@" << it << " # " << n_nodes << ": global " 
    << lb << " < " << ub << " " << " |.| " << fabs(ub - lb) << "\t selected "
    << node_star.GetLB() << " < " << node_star.GetUB() << " "
    << node_star.GetLevel() << " " << node_star.GetIds()[0] << "\t" 
    << node_star.GetTetrahedron().GetCenter().transpose() 
    << std::endl;
  return node_star;
}

}
