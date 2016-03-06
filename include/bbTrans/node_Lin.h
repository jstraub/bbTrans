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
#include "bbTrans/node_R3.h"
#include "bbTrans/node_S3.h"
#include "bbTrans/box.h"
#include "bbTrans/tetrahedron.h"

namespace bb {

class NodeLin : public BaseNode {
 public:
  NodeLin(const Box& box, std::vector<uint32_t> ids);
  NodeLin(const NodeLin& node);
  virtual ~NodeLin() = default;

  /// Give access to the interior Tetrahedron of the partitining of the cube
  const Tetrahedron4D& GetTetrahedron() const {return nodeS3s_[4].GetTetrahedron();}
  void SetLbArgument(const Eigen::Quaterniond& q) {q_lb_ = q;}
  Eigen::Quaterniond GetLbArgument() const {return q_lb_;}
  virtual uint32_t GetBranchingFactor(uint32_t i) const { return 8;}
  virtual std::string ToString() const;
  virtual std::string Serialize() const;
  std::string GetSpace() const { return "Lin"; }
  double GetVolume() const;
  const NodeS3& GetNodeS3(uint32_t i) const { return nodeS3s_[i]; }
  NodeS3& GetNodeS3(uint32_t i) { return nodeS3s_[i]; }
 protected:
  NodeR3 nodeLin_;
  std::vector<NodeS3> nodeS3s_;
  Eigen::Quaterniond q_lb_;
  virtual void Linearize(const Box& box);
  virtual Tetrahedron4D TetraFromBox(const Box& box, uint32_t i0, uint32_t i1,
    uint32_t i2, uint32_t i3) = 0;
};
}
