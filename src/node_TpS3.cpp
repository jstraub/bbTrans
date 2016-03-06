/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "bbTrans/node_TpS3.h"

namespace bb {

NodeTpS3::NodeTpS3(const Box& box, std::vector<uint32_t> ids) 
  : NodeAA(box, ids)
{ }

NodeTpS3::NodeTpS3(const NodeTpS3& node) 
  : BaseNode(node.GetIds(), node.GetLB(), node.GetUB()),
  nodeTpS3_(node.nodeTpS3_), nodeS3s_(node.nodeS3s_), q_lb_(node.q_lb_)
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


}
