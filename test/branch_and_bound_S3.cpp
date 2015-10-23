/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include <vector>
#include "optRot/node.h"
#include "optRot/lower_bound_log.h"
#include "optRot/upper_bound_log.h"
#include "optRot/upper_bound_convexity_log.h"
#include "optRot/branch_and_bound.h"

using namespace OptRot;

int main(int argc, char** argv) {

  Eigen::Quaterniond q_true(1.,1.,1.,1.);
  q_true.normalize();
  std::cout << "true quaternion: " << q_true.coeffs().transpose() << std::endl;
  //std::cout << q_true.toRotationMatrix() << std::endl;
  
  Eigen::Vector3d muA1, muA2;
  muA1 << 1.,0.,0.;
  muA2 << 0.,1.,0.;
  Eigen::Vector3d muB1, muB2;
  muB1 = q_true.inverse()._transformVector(muA1);
  muB2 = q_true.inverse()._transformVector(muA2);

  std::vector<vMF<3>> vmfs_A;
  vmfs_A.push_back(vMF<3>(muA1,100.,0.3));
  vmfs_A.push_back(vMF<3>(muA2,10.,0.7));

  std::vector<vMF<3>> vmfs_B;
  vmfs_B.push_back(vMF<3>(muB1,100.,0.3));
  vmfs_B.push_back(vMF<3>(muB2,10.,0.7));
  
  vMFMM<3> vmf_mm_A(vmfs_A);
  vMFMM<3> vmf_mm_B(vmfs_B);

  std::list<NodeS3> nodes = GenerateNotesThatTessellateS3();
  LowerBoundLog lower_bound(vmf_mm_A, vmf_mm_B);
  UpperBoundLog upper_bound(vmf_mm_A, vmf_mm_B);
  UpperBoundConvexityLog upper_bound_convexity(vmf_mm_A, vmf_mm_B);
  
  double eps = 1.0e-5 * M_PI / 180.;
  uint32_t max_it = 1000;
  BranchAndBound<NodeS3> bb(lower_bound, upper_bound_convexity);
  NodeS3 node_star = bb.Compute(nodes, eps, max_it);

  std::cout << "optimum quaternion: " 
    << node_star.GetTetrahedron().GetCenter().transpose()
    << std::endl;
}

