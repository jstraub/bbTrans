/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include "optRot/upper_bound_log.h"

namespace OptRot {

UpperBoundLog::UpperBoundLog(const vMFMM<3>& vmf_mm_A, const vMFMM<3>&
    vmf_mm_B) : vmf_mm_A_(vmf_mm_A), vmf_mm_B_(vmf_mm_B)
{}

double UpperBoundLog::Evaluate(const Node& node) {
  std::vector<Eigen::Quaterniond> qs(4);
  for (uint32_t i=0; i<4; ++i)
    qs[i] = node.GetTetrahedron().GetVertex(i);
   
  Eigen::VectorXd ubElem(vmf_mm_A_.GetK()*vmf_mm_B_.GetK());
  for (std::size_t j=0; j < vmf_mm_A_.GetK(); ++j)
    for (std::size_t k=0; k < vmf_mm_B_.GetK(); ++k) {
      Eigen::Vector3d p_star = ClosestPointInTetrahedron(vmf_mm_A_.Get(j),
          vmf_mm_B_.Get(k), qs);
      ubElem(j*vmf_mm_B_.GetK() + k) =
        ComputeLogvMFtovMFcost<3>(vmf_mm_A_.Get(j), vmf_mm_B_.Get(k),
            p_star);
    }
  return SumExp(ubElem);
}

}
