/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <stdint.h>
#include <vector>
#include <list>
#include <sstream>
#include <string>

#include "manifold/S.h"
#include "bbTrans/node_S3.h"
#include "bbTrans/box.h"
#include "bbTrans/tetrahedron.h"
#include "bbTrans/s3_tessellation.h"

namespace bb {

class NodeTpS3 : public BaseNode {
 public:
  NodeTpS3(const Box& box, std::vector<uint32_t> ids);
  NodeS3(const NodeS3& node);
  virtual ~NodeS3() = default;
  virtual std::vector<NodeS3> Branch() const;
//  const Tetrahedron4D& GetTetrahedron() const { return tetrahedron_;}
  void SetLbArgument(const Eigen::Quaterniond& q) {q_lb_ = q;}
  Eigen::Quaterniond GetLbArgument() const {return q_lb_;}
  virtual uint32_t GetBranchingFactor(uint32_t i) const { return 8;}
  virtual std::string ToString() const;
  virtual std::string Serialize() const;
  std::string GetSpace() const { return "TpS3"; }
  double GetVolume() const { return tetrahedron_.GetVolume();}
 protected:
  NodeR3 nodeTpS3_;
  std::vector<NodeS3> nodeS3s_;
  Eigen::Quaterniond q_lb_;
};

std::list<NodeS3> GenerateNotesThatTessellateS3();
}
