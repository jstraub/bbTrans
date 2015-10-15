/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include "optRot/lower_bound_log.h"

namespace OptRot {

LowerBoundLog::LowerBoundLog(const vMFMM<3>& vmf_mm_A, const vMFMM<3>&
    vmf_mm_B) 
  : vmf_mm_A_(vmf_mm_A), vmf_mm_B_(vmf_mm_B)
{}

double LowerBoundLog::Evaluate(const Node& node) {
  Eigen::VectorXd lbs(5); 
  std::vector<Eigen::Quaterniond> qs(5);
  qs[0] = node.GetTetrahedron().GetCenterQuaternion();
  for (uint32_t i=0; i<4; ++i)
    qs[i+1] = node.GetTetrahedron().GetVertexQuaternion(i);
  for (uint32_t i=0; i<5; ++i) {
    Eigen::VectorXd lbElem(vmf_mm_A_.GetK()*vmf_mm_B_.GetK());
//    std::cout << qs[i].vec().transpose() << " " << qs[i].w() << std::endl;
    for (std::size_t j=0; j < vmf_mm_A_.GetK(); ++j) {
//      std::cout << "vMF " << vmf_mm_A_.Get(j).GetMu().transpose()
//        << " " << vmf_mm_A_.Get(j).GetTau() << " " << vmf_mm_A_.Get(j).GetPi()
//        << std::endl;
      for (std::size_t k=0; k < vmf_mm_B_.GetK(); ++k) {
//        std::cout << "vMF " << vmf_mm_B_.Get(k).GetMu().transpose()
//          << " " << vmf_mm_B_.Get(k).GetTau() << " " << vmf_mm_B_.Get(k).GetPi()
//          << std::endl;
        lbElem(j*vmf_mm_B_.GetK() + k) =
          ComputeLogvMFtovMFcost<3>(vmf_mm_A_.Get(j), vmf_mm_B_.Get(k),
              qs[i]._transformVector(vmf_mm_B_.Get(k).GetMu()));
//              qs[i].toRotationMatrix()*vmf_mm_B_.Get(k).GetMu());
//        std::cout << lbElem(j*vmf_mm_B_.GetK() + k) << std::endl;
//        std::cout << "rotated muB " << qs[i]._transformVector(vmf_mm_B_.Get(k).GetMu()).transpose()
//          << std::endl;
      }
    }
//    std::cout << lbElem.transpose() << std::endl;
    lbs(i) = SumExp(lbElem);
  }
  return lbs.maxCoeff();
}

}
