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
  nodeLin_(node.nodeLin_), nodeS3s_(node.nodeS3s_), q_lb_(node.q_lb_)
{ }

void NodeLin::Linearize(const Box& box) {
  // subdivide box in Lin space into 4 tetrahedra
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

double NodeLin::GetVolume() const { 
  double V = 0.;
  for (uint32_t i=0; i<5; ++i)
    V += nodeS3s_[i].GetVolume();
  return V;
}

std::string NodeLin::ToString() const {
  std::stringstream ss;
  ss << " V=" << GetVolume()
    << " in lin space: " << nodeLin_.ToString() 
    << std::endl;
  for (const auto& node : nodeS3s_) 
    ss << "\t " << node.ToString() << std::endl;
  return ss.str();
};

std::string NodeLin::Serialize() const {
  return nodeLin_.Serialize();
};

}
