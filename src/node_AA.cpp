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

Eigen::Quaterniond NodeAA::Project(const Eigen::Vector3d& c) const {
  Eigen::Quaterniond q;
  double theta = c.norm();
  if (theta > 1e-9) {
    q.w() = cos(theta*0.5);
    q.vec() = -c*sin(theta*0.5)/theta;
  } else {
    q.w() = 1.;
    q.vec().fill(0.);
  }
  return q;
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
  return nodes;
}

}
