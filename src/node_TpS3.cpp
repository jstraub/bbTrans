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

Tetrahedron4D NodeTpS3::TetraFromBox(const Box& box, uint32_t i0, uint32_t i1,
    uint32_t i2, uint32_t i3) {
  Eigen::Vector4d north;
  north << 0.,0.,0.,1.;
  S<double,4> s3(north);
  std::vector<Eigen::Vector3d> cs(4);
  Eigen::Matrix4d qs;
  box.GetCorner(i0, cs[0]); box.GetCorner(i1, cs[1]);
  box.GetCorner(i2, cs[2]); box.GetCorner(i3, cs[3]);

  for (uint32_t i=0; i<4; ++i)
    qs.col(i) = s3.Exp(s3.ToAmbient(cs[i])).vector();
  return Tetrahedron4D(qs);
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
