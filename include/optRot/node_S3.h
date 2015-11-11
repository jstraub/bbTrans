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
    std::stringstream out; 
    out << GetTetrahedron().GetCenter().transpose() << std::endl
      << GetTetrahedron().GetVertex(0).transpose() << std::endl
      << GetTetrahedron().GetVertex(1).transpose() << std::endl
      << GetTetrahedron().GetVertex(2).transpose() << std::endl
      << GetTetrahedron().GetVertex(3).transpose() << std::endl;
    out << "pairwise angles: ";
    for (uint32_t i=0; i < 4; ++i) 
      for (uint32_t j=0; j < 4; ++j) 
        if(i!=j)
          out << i << "," <<j<< ": "
            << GetTetrahedron().GetVertexQuaternion(i).angularDistance(
                GetTetrahedron().GetVertexQuaternion(j)) *180./M_PI<< " ";
    return out.str();
  };
  virtual std::string Serialize() const {
    std::stringstream out; 
    Eigen::Vector4d v;
    for (uint32_t i=0; i<4; ++i) {
      v = GetTetrahedron().GetVertex(i);
      out << v(0) << " " << v(1) << " " << v(2) << " " << v(3) << std::endl;
    }
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
