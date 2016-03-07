/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "bbTrans/node_TpS3.h"

namespace bb {

NodeTpS3::NodeTpS3(const Box& box, std::vector<uint32_t> ids) 
  : NodeLin(box, ids)
{ 
  Linearize(box);
}

NodeTpS3::NodeTpS3(const NodeTpS3& node) 
  : NodeLin(node)
{ }

Eigen::Quaterniond NodeTpS3::Project(const Eigen::Vector3d& c) const {
  // In accordance with the other tessellation approaches
  static Eigen::Vector4d north(1.,0.,0.,0.);
  static S<double,4> s3(north);
  Eigen::Vector4d qvec = s3.Exp(s3.ToAmbient(c)).vector();
  return Eigen::Quaterniond(qvec(0),qvec(1),qvec(2),qvec(3));
}

std::vector<NodeTpS3> NodeTpS3::Branch() const {
  std::vector<NodeTpS3> nodes;
  nodes.reserve(8);
  std::vector<NodeR3> boxs = nodeLin_.Branch();
  for (uint32_t i=0; i < boxs.size(); ++i) {
    std::vector<uint32_t> ids(this->ids_);
    ids.push_back(i);
    nodes.push_back(NodeTpS3(boxs[i].GetBox(), ids));
  }
  return nodes;
}


}
