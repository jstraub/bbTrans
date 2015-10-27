/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "optRot/node.h"

namespace OptRot {

BaseNode::BaseNode(uint32_t lvl, std::vector<uint32_t> ids) :
  lvl_(lvl), ids_(ids), lb_(-1e12), ub_(1e12)
{}

BaseNode::BaseNode(uint32_t lvl, std::vector<uint32_t> ids, double lb,
    double ub) : lvl_(lvl), ids_(ids), lb_(lb), ub_(ub) {
}

BaseNode::BaseNode(const BaseNode& node) : lvl_(node.GetLevel()),
  ids_(node.GetIds()), lb_(node.GetLB()), ub_(node.GetUB()) {
}

uint64_t BaseNode::GetIdAtLevel(uint32_t lvl) const {
  uint64_t id = 0;
  uint32_t factor = 1;
  for (uint32_t i=0; i< std::min(ids_.size(), size_t(lvl)); ++i) {
    id += ids_[i]*factor;
    factor *= GetBranchingFactor(i);
  } 
  return id;
}

}
