/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "bbTrans/node_TpS3.h"

namespace bb {

NodeTpS3::NodeTpS3(const Box& box, std::vector<uint32_t> ids) 
  : BaseNode(ids), nodeTpS3_(box, ids)
{ 
  // subdivide box in TpS into 4 tetrahedra
  // https://www.ics.uci.edu/~eppstein/projects/tetra/
  S<double,4> s3;
  nodeS3s_.reserve(5);
  std::vector<Eigen::Vector3d> cs(4);
  Eigen::Vector4d north;
  north << 0.,0.,0.,1.;
  Eigen::Matrix4d qs;
  // NodeS3 1: 0 4 5 7
  box.GetCorner(0, cs[0]); box.GetCorner(4, cs[1]);
  box.GetCorner(5, cs[2]); box.GetCorner(7, cs[3]);
  for (uint32_t i=0; i<4; ++i)
    qs.col(i) = s3.Exp(north, cs[i]);
  Tetrahedron4D t(qs);
  std::vector<uint32_t> idsInternal(1,0);
  nodeS3s_.push_back(NodeS3(t, idsInternal));
  // NodeS3 2: 1 4 5 6
  box.GetCorner(1, cs[0]); box.GetCorner(4, cs[1]);
  box.GetCorner(5, cs[2]); box.GetCorner(6, cs[3]);
  for (uint32_t i=0; i<4; ++i)
    qs.col(i) = s3.Exp(north, cs[i]);
  t = Tetrahedron4D(qs);
  idsInternal = std::vector<uint32_t>(1,1);
  nodeS3s_.push_back(NodeS3(t, idsInternal));
  // NodeS3 3: 2 4 6 7
  box.GetCorner(2, cs[0]); box.GetCorner(4, cs[1]);
  box.GetCorner(6, cs[2]); box.GetCorner(7, cs[3]);
  for (uint32_t i=0; i<4; ++i)
    qs.col(i) = s3.Exp(north, cs[i]);
  t = Tetrahedron4D(qs);
  idsInternal = std::vector<uint32_t>(1,2);
  nodeS3s_.push_back(NodeS3(t, idsInternal));
  // NodeS3 4: 3 5 6 7
  box.GetCorner(3, cs[0]); box.GetCorner(5, cs[1]);
  box.GetCorner(6, cs[2]); box.GetCorner(7, cs[3]);
  for (uint32_t i=0; i<4; ++i)
    qs.col(i) = s3.Exp(north, cs[i]);
  t = Tetrahedron4D(qs);
  idsInternal = std::vector<uint32_t>(1,3);
  nodeS3s_.push_back(NodeS3(t, idsInternal));
  // NodeS3 5: 0 1 2 3
  box.GetCorner(0, cs[0]); box.GetCorner(1, cs[1]);
  box.GetCorner(2, cs[2]); box.GetCorner(3, cs[3]);
  for (uint32_t i=0; i<4; ++i)
    qs.col(i) = s3.Exp(north, cs[i]);
  t = Tetrahedron4D(qs);
  idsInternal = std::vector<uint32_t>(1,4);
  nodeS3s_.push_back(NodeS3(t, idsInternal));
}

NodeTpS3::NodeTpS3(const NodeTpS3& node) 
  : BaseNode(node.GetIds(), node.GetLB(), node.GetUB()),
  nodeTpS3_(node.nodeTpS3_), nodeS3s_(node.nodeS3s_), q_lb_(node.q_lb_)
{ }

std::vector<NodeTpS3> NodeTpS3::Branch() const {
  std::vector<NodeTpS3> nodes;
  nodes.reserve(8);
  std::vector<NodeR3> boxs = nodeTpS3_.Branch();
  for (uint32_t i=0; i < boxs.size(); ++i) {
    std::vector<uint32_t> ids(this->ids_);
    ids.push_back(i);
    nodes.push_back(NodeTpS3(boxs[i].GetBox(), ids));
  }
  return nodes;
}

double NodeTpS3::GetVolume() const { 
  double V = 0.;
  for (uint32_t i=0; i<5; ++i)
    V += nodeS3s_[i].GetVolume();
  return V;
}

std::string NodeTpS3::ToString() const {
  return nodeTpS3_.ToString();
};

std::string NodeTpS3::Serialize() const {
  return nodeTpS3_.Serialize();
};


}
