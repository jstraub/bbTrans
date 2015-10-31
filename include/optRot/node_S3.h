/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <stdint.h>
#include <vector>
#include <list>
#include <sstream>
#include <string>
#include "optRot/node.h"
#include "optRot/tetrahedron.h"
#include "optRot/s3_tessellation.h"

namespace OptRot {

class NodeS3 : public BaseNode {
 public:
  NodeS3(const Tetrahedron4D& tetrahedron, std::vector<uint32_t> ids);
  NodeS3(const NodeS3& node);
  virtual ~NodeS3() = default;
  virtual std::vector<NodeS3> Branch() const;
  const Tetrahedron4D& GetTetrahedron() const { return tetrahedron_;}
  void SetLbArgument(const Eigen::Quaterniond& q) {q_lb_ = q;}
  Eigen::Quaterniond GetLbArgument() const
  {return q_lb_;}
  virtual uint32_t GetBranchingFactor(uint32_t i) const { return i==0? 600 : 8;}
  virtual std::string ToString() const {
    std::stringstream out; out << GetTetrahedron().GetCenter().transpose();
    return out.str();
  };
  std::string GetSpace() const { return "S3"; }
  double GetVolume() const { return tetrahedron_.GetVolume();}
 protected:
  Tetrahedron4D tetrahedron_;
  Eigen::Quaterniond q_lb_;
};

std::list<NodeS3> GenerateNotesThatTessellateS3();
}
