/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include <random>
#include <vector>
#include "optRot/node.h"
#include "optRot/lower_bound_log.h"
#include "optRot/upper_bound_log.h"
#include "optRot/upper_bound_convexity_log.h"

using namespace OptRot;

Eigen::Quaterniond RandomRotation() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> N(0,1);
  Eigen::Quaterniond q(N(gen), N(gen), N(gen), N(gen));
  q.normalize();
  return q;
}

int main(int argc, char ** argv) {

  Eigen::Quaterniond q_true = RandomRotation();
  Eigen::Vector3d muA1, muA2;
  muA1 << 1.,0.,0.;
  muA2 << 0.,1.,0.;
  Eigen::Vector3d muB1, muB2;
  muB1 = q_true.toRotationMatrix()*muA1;
  muB2 = q_true.toRotationMatrix()*muA2;

  std::vector<vMF<3>> vmfs_A;
  vmfs_A.push_back(vMF<3>(muA1,1.,0.3));
  vmfs_A.push_back(vMF<3>(muA2,100.,0.7));

  std::vector<vMF<3>> vmfs_B;
  vmfs_B.push_back(vMF<3>(muB1,1.,0.3));
  vmfs_B.push_back(vMF<3>(muB2,100.,0.7));
  
  vMFMM<3> vmf_mm_A(vmfs_A);
  vMFMM<3> vmf_mm_B(vmfs_B);

  std::vector<Node> nodes = GenerateNotesThatTessellateS3();
  LowerBoundLog lower_bound(vmf_mm_A, vmf_mm_B);
  UpperBoundLog upper_bound(vmf_mm_A, vmf_mm_B);
  UpperBoundConvexityLog upper_bound_convexity(vmf_mm_A, vmf_mm_B);
  Eigen::VectorXd lbs(nodes.size());
  Eigen::VectorXd ubs(nodes.size());
  Eigen::VectorXd ubCs(nodes.size());
  for (std::size_t i=0; i<nodes.size(); ++i) {
    lbs[i] = lower_bound.Evaluate(nodes[i]);
    ubs[i] = upper_bound.Evaluate(nodes[i]);
    ubCs[i] = upper_bound_convexity.Evaluate(nodes[i]);
  }
  std::cout << "Lower Bounds: " << std::endl;
  std::cout << lbs.transpose() << std::endl;
  std::cout << "Upper Bounds: " << std::endl;
  std::cout << ubs.transpose() << std::endl;
  std::cout << "Upper Bounds Convexity: " << std::endl;
  std::cout << ubCs.transpose() << std::endl;
  return 0;
}
