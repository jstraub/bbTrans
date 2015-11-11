/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <stdint.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include "optRot/node.h"
#include "optRot/box.h"

namespace OptRot {

class NodeR3 : public BaseNode {
 public:
  NodeR3(const Box& box, std::vector<uint32_t> ids);
  NodeR3(const NodeR3& node);
  virtual ~NodeR3() = default;
  virtual std::vector<NodeR3> Branch() const;
  const Box& GetBox() const { return box_;}
  const Eigen::Vector3d& GetLbArgument() const {return t_lb_;}
  void SetLbArgument(const Eigen::Vector3d& t) {t_lb_ = t;}
  virtual uint32_t GetBranchingFactor(uint32_t i) const { return 8;}
  virtual std::string ToString() const {
    std::stringstream out; out << GetBox().GetCenter().transpose();
    return out.str();
  };
  virtual std::string Serialize() const {
    std::stringstream out; 
    Eigen::Vector3d c;
    for (uint32_t i=0; i<8; ++i) {
      GetBox().GetCorner(i, c);
      out << c(0) << " " << c(1) << " " << c(2) << std::endl;
    }
    return out.str();
  };
  std::string GetSpace() const { return "R3"; }
  double GetVolume() const { return box_.GetVolume();}
 protected:
  Box box_;
  Eigen::Vector3d t_lb_;
};

std::list<NodeR3> GenerateNotesThatTessellateR3(const Eigen::Vector3d&
    min, const Eigen::Vector3d& max, double max_side_len); 
}
