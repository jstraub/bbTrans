/* Copyright (c) 2016, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "bbTrans/node_Lin.h"

namespace bb {

NodeLin::NodeLin(const Box& box, std::vector<uint32_t> ids) 
  : BaseNode(ids), nodeLin_(box, ids)
{ }

NodeLin::NodeLin(const NodeLin& node) 
  : BaseNode(node.GetIds(), node.GetLB(), node.GetUB()),
  nodeLin_(node.nodeLin_), qs_(node.qs_), q_lb_(node.q_lb_)
{ }

//void NodeLin::Linearize(const Box& box) {
//  // subdivide box in Lin space into 4 tetrahedra
//  // https://www.ics.uci.edu/~eppstein/projects/tetra/
//  nodeS3s_.reserve(5);
//  // NodeS3 1: 0 4 5 7
//  Tetrahedron4D t = TetraFromBox(box, 0, 4, 5, 7);
//  std::vector<uint32_t> idsInternal(1,0);
//  nodeS3s_.push_back(NodeS3(t, idsInternal));
//  // NodeS3 2: 1 4 5 6
//  t = TetraFromBox(box, 1, 4, 5, 6);
//  idsInternal = std::vector<uint32_t>(1,1);
//  nodeS3s_.push_back(NodeS3(t, idsInternal));
//  // NodeS3 3: 2 4 6 7
//  t = TetraFromBox(box, 2, 4, 6, 7);
//  idsInternal = std::vector<uint32_t>(1,2);
//  nodeS3s_.push_back(NodeS3(t, idsInternal));
//  // NodeS3 4: 3 5 6 7
//  t = TetraFromBox(box, 3, 5, 6, 7);
//  idsInternal = std::vector<uint32_t>(1,3);
//  nodeS3s_.push_back(NodeS3(t, idsInternal));
//  // NodeS3 5: 0 1 2 3
//  t = TetraFromBox(box, 0, 1, 2, 3);
//  idsInternal = std::vector<uint32_t>(1,4);
//  nodeS3s_.push_back(NodeS3(t, idsInternal));
//}

void NodeLin::Linearize(const Box& box) {
  qs_.reserve(8);
  for (uint32_t i=0; i<8; ++i) {
    Eigen::Vector3d c;
    box.GetCorner(i,c);
    qs_.push_back(Project(c));
  }
}

Eigen::Quaterniond NodeLin::GetCenter() const {
  return Project(nodeLin_.GetBox().GetCenter());
}

double NodeLin::GetVolume() const { 
  // TODO would be better to use space on the surface of the sphere
  return nodeLin_.GetVolume();
}

std::string NodeLin::ToString() const {
  std::stringstream ss;
  ss << " V=" << GetVolume()
    << " in lin space: " << nodeLin_.ToString() 
    << std::endl;
  for (const auto& q : qs_) 
    ss << "\t " << q.coeffs().transpose() << std::endl;
  return ss.str();
};

std::string NodeLin::Serialize() const {
  return nodeLin_.Serialize();
};

}
