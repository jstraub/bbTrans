/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "bbTrans/node_AA.h"

namespace bb {

NodeAA::NodeAA(const Box& box, std::vector<uint32_t> ids) 
  : NodeLin(box,ids)
{
  Linearize(box);
}

NodeAA::NodeAA(const NodeAA& node) 
  : NodeLin(node)
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
//      qs(3,i) = cos(theta*0.5);
//      qs.block<3,1>(0,i) = -cs[i]*sin(theta*0.5)/theta;
    } else {
//      qs(3,i) = 1.;
//      qs.block<3,1>(0,i).fill(0);
      qs(0,i) = 1.;
      qs.block<3,1>(1,i).fill(0);
    }
  }
  return Tetrahedron4D(qs);
}

std::vector<NodeAA> NodeAA::Branch() const {
  std::vector<NodeAA> nodes;
  nodes.reserve(8);
  std::vector<NodeR3> boxs = nodeLin_.Branch();
  for (uint32_t i=0; i < boxs.size(); ++i) {
    std::vector<uint32_t> ids(this->ids_);
    ids.push_back(i);
    nodes.push_back(NodeAA(boxs[i].GetBox(), ids));
  }
  return nodes;
}

std::list<NodeAA> TessellateAA() {
  Eigen::Vector3d p_min(-M_PI,-M_PI,-M_PI);
  Eigen::Vector3d p_max(M_PI,M_PI,M_PI);
  NodeAA root(Box(p_min, p_max),std::vector<uint32_t>(0));
//  std::cout << root.ToString() << std::endl;
  std::vector<NodeAA> l1 = root.Branch();
  std::list<NodeAA> nodes;
  for (auto& node1 : l1) {
    std::vector<NodeAA> l2 = node1.Branch();
//    for (auto& node2 : l2) {
//      std::vector<NodeAA> l3 = node2.Branch();
//      for (auto& node3 : l3) {
//        std::vector<NodeAA> l4 = node3.Branch();
        nodes.insert(nodes.end(), l2.begin(), l2.end());
//      }
//    }
  }

}

}
