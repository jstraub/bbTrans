/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include <vector>
#include "bbTrans/node.h"
#include "bbTrans/lower_bound_S3.h"
#include "bbTrans/upper_bound_indep_S3.h"
#include "bbTrans/upper_bound_convex_S3.h"
#include "bbTrans/branch_and_bound.h"
#include "jsCore/timer.hpp"

using namespace bb;

Eigen::Quaterniond RandomRotation() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> N(0,1);
  Eigen::Quaterniond q(N(gen), N(gen), N(gen), N(gen));
  q.normalize();
  return q;
}

int main(int argc, char** argv) {

  Eigen::Quaterniond q_true(1.,1.,1.,1.);
  q_true.normalize();
  q_true = RandomRotation();

  q_true = Eigen::Quaterniond(-0.285745, -0.234756, -0.690839, -0.621274);
//  std::cout << "true quaternion: " << q_true.coeffs().transpose() << std::endl;
  std::cout << "true quaternion: " << q_true.w() << " " << q_true.x() << " " << q_true.y() << " " << q_true.z() << std::endl;
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
  LowerBoundS3 lower_bound(vmf_mm_A, vmf_mm_B);
  UpperBoundIndepS3 upper_bound(vmf_mm_A, vmf_mm_B);
  UpperBoundConvexS3 upper_bound_convex(vmf_mm_A, vmf_mm_B);
  
  double eps = 1.0e-8;
  uint32_t max_it = 1000;
  uint32_t max_lvl = 15;
  BranchAndBound<NodeS3> bb(lower_bound, upper_bound_convex);

  jsc::Timer t0;
  NodeS3 node_star = bb.Compute(nodes, eps, max_lvl, max_it);
  t0.toctic("BB ");

  std::cout << "optimum quaternion: " 
    << node_star.GetLbArgument().coeffs().transpose()
    << " |.|=" << node_star.GetLbArgument().norm()
    << std::endl;
  std::cout << "true quaternion: " << q_true.coeffs().transpose()
    << " |.|=" << q_true.norm()
    << std::endl;

  std::cout << "Deviation from GT: "  <<
    q_true.angularDistance(node_star.GetLbArgument())
//    << 2.*acos(std::min(1.0, std::max(-1.,node_star.GetLbArgument().dot(q_true))))*180./M_PI
    << std::endl;
}

