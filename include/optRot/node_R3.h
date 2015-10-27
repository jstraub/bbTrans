/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <stdint.h>
#include <iostream>
#include <vector>
#include <list>
#include "optRot/node.h"
#include "optRot/box.h"

namespace OptRot {

class NodeR3 : public BaseNode {
 public:
  NodeR3(const Box& box, uint32_t lvl,
      std::vector<uint32_t> ids);
  NodeR3(const NodeR3& node);
  virtual ~NodeR3() = default;
  virtual std::vector<NodeR3> Branch() const;
  const Box& GetBox() const { return box_;}
 protected:
  Box box_;
};

std::list<NodeR3> GenerateNotesThatTessellateR3(const Eigen::Vector3d&
    min, const Eigen::Vector3d& max, double max_side_len); 
}
