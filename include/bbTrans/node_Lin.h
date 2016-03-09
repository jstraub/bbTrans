/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <stdint.h>
#include <vector>
#include <list>
#include <sstream>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Geometry>

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

  void SetLbArgument(const Eigen::Quaterniond& q) {q_lb_ = q;}
  Eigen::Quaterniond GetLbArgument() const {return q_lb_;}
  Eigen::Quaterniond GetCenter() const;
  virtual uint32_t GetBranchingFactor(uint32_t i) const { return 8;}
  virtual std::string ToString() const;
  virtual std::string Serialize() const;
  std::string GetSpace() const { return "Lin"; }
  double GetVolume() const;
  NodeS3 GetNodeS3() const;
  const std::vector<Eigen::Quaterniond>& GetQuaternions() const { return qs_; }
  std::vector<Eigen::Quaterniond>& GetQuaternions() { return qs_; }
 protected:
  NodeR3 nodeLin_;
  std::vector<Eigen::Quaterniond> qs_;
  Eigen::Quaterniond q_lb_;
  virtual void Linearize(const Box& box);
  virtual Eigen::Quaterniond Project(const Eigen::Vector3d& c) const = 0;
};
}
