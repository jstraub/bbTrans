/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "bbTrans/node_AA.h"

namespace bb {

NodeAA::NodeAA(const Box& box, std::vector<uint32_t> ids) 
  : BaseNode(ids), nodeAA_(box, ids)
{ 
  // subdivide box in AA space into 4 tetrahedra
  // https://www.ics.uci.edu/~eppstein/projects/tetra/
  nodeS3s_.reserve(5);
  // NodeS3 1: 0 4 5 7
  Tetrahedron4D t = TetraFromBox(box, 0, 4, 5, 7);
  std::vector<uint32_t> idsInternal(1,0);
  nodeS3s_.push_back(NodeS3(t, idsInternal));
  // NodeS3 2: 1 4 5 6
  t = TetraFromBox(box, 1, 4, 5, 6);
  idsInternal = std::vector<uint32_t>(1,1);
  nodeS3s_.push_back(NodeS3(t, idsInternal));
  // NodeS3 3: 2 4 6 7
  t = TetraFromBox(box, 2, 4, 6, 7);
  idsInternal = std::vector<uint32_t>(1,2);
  nodeS3s_.push_back(NodeS3(t, idsInternal));
  // NodeS3 4: 3 5 6 7
  t = TetraFromBox(box, 3, 5, 6, 7);
  idsInternal = std::vector<uint32_t>(1,3);
  nodeS3s_.push_back(NodeS3(t, idsInternal));
  // NodeS3 5: 0 1 2 3
  t = TetraFromBox(box, 0, 1, 2, 3);
  idsInternal = std::vector<uint32_t>(1,4);
  nodeS3s_.push_back(NodeS3(t, idsInternal));
}

NodeAA::NodeAA(const NodeAA& node) 
  : BaseNode(node.GetIds(), node.GetLB(), node.GetUB()),
  nodeAA_(node.nodeAA_), nodeS3s_(node.nodeS3s_), q_lb_(node.q_lb_)
{ }

Tetrahedron4D NodeAA::TetraFromBox(const Box& box, uint32_t i0, uint32_t i1,
    uint32_t i2, uint32_t i3) {
  std::vector<Eigen::Vector3d> cs(4);
  Eigen::Matrix4d qs;
  box.GetCorner(i0, cs[0]); box.GetCorner(i1, cs[1]);
  box.GetCorner(i2, cs[2]); box.GetCorner(i3, cs[3]);
  for (uint32_t i=0; i<4; ++i) {
    double theta = cs[i].norm();
    if (theta > 1e-9) {
      qs(0,i) = cos(theta*0.5);
      qs.block<3,1>(1,i) = -cs[i]*sin(theta*0.5)/theta;
    } else {
      qs(0,i) = 1.;
      qs.block<3,1>(1,i).fill(0);
    }
  }
  return Tetrahedron4D(qs);
}

std::vector<NodeAA> NodeAA::Branch() const {
  std::vector<NodeAA> nodes;
  nodes.reserve(8);
  std::vector<NodeR3> boxs = nodeAA_.Branch();
  for (uint32_t i=0; i < boxs.size(); ++i) {
    std::vector<uint32_t> ids(this->ids_);
    ids.push_back(i);
    nodes.push_back(NodeAA(boxs[i].GetBox(), ids));
  }
  return nodes;
}

double NodeAA::GetVolume() const { 
  double V = 0.;
  for (uint32_t i=0; i<5; ++i)
    V += nodeS3s_[i].GetVolume();
  return V;
}

std::string NodeAA::ToString() const {
  std::stringstream ss;
  ss << " V=" << GetVolume()
    << " in lin space: " << nodeAA_.ToString() 
    << std::endl;
  for (const auto& node : nodeS3s_) 
    ss << "\t " << node.ToString() << std::endl;
  return ss.str();
};

std::string NodeAA::Serialize() const {
  return nodeAA_.Serialize();
};


}
