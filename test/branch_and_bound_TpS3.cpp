/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include <vector>
#include "bbTrans/node.h"
#include "bbTrans/lower_bound_Lin.h"
#include "bbTrans/upper_bound_Lin.h"
#include "bbTrans/branch_and_bound.h"

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

  LowerBoundS3 lower_bound_S3(vmf_mm_A, vmf_mm_B);
  UpperBoundIndepS3 upper_bound_S3(vmf_mm_A, vmf_mm_B);
  UpperBoundConvexS3 upper_bound_convex_S3(vmf_mm_A, vmf_mm_B);

  LowerBoundTpS3 lower_bound(lower_bound_S3);
  UpperBoundIndepTpS3 upper_bound_indep(upper_bound_S3);
  UpperBoundConvexTpS3 upper_bound_convex(upper_bound_convex_S3);

  Eigen::Vector3d p_min(-M_PI*0.5,-M_PI*0.5,-M_PI*0.5);
  Eigen::Vector3d p_max( M_PI*0.5, M_PI*0.5, M_PI*0.5);
  NodeTpS3 root(Box(p_min, p_max),std::vector<uint32_t>(0));
//  std::cout << root.ToString() << std::endl;
  std::vector<NodeTpS3> l1 = root.Branch();
  std::list<NodeTpS3> nodes;
  for (auto& node1 : l1) {
    std::vector<NodeTpS3> l2 = node1.Branch();
//    for (auto& node2 : l2) {
//      std::vector<NodeTpS3> l3 = node2.Branch();
//      for (auto& node3 : l3) {
//        std::vector<NodeTpS3> l4 = node3.Branch();
        nodes.insert(nodes.end(), l2.begin(), l2.end());
//      }
//    }
  }
  std::cout << "initial # nodes: " << nodes.size() << std::endl;

//  for (const auto& node : nodes) {
//    std::cout << node.ToString() << std::endl;
//  }
  
  double eps = 1.0e-8;
  uint32_t max_it = 1000;
  uint32_t max_lvl = 20;
  BranchAndBound<NodeTpS3> bb(lower_bound, upper_bound_convex);
  NodeTpS3 node_star = bb.Compute(nodes, eps, max_lvl, max_it);

  std::cout << "optimum quaternion: " 
    << " w=" << node_star.GetLbArgument().w()
    << " " << node_star.GetLbArgument().vec().transpose()
    << "|.| " << node_star.GetLbArgument().norm()
    << std::endl;

  std::cout << "true quaternion: " << q_true.w() << " " << q_true.x()
    << " " << q_true.y() << " " << q_true.z() 
    << "|.| " << q_true.norm()
    << std::endl;

  std::cout << "Deviation from GT: " 
    << 2.*acos(std::min(1.0, std::max(-1.,node_star.GetLbArgument().dot(q_true))))*180./M_PI
    << std::endl;
}

