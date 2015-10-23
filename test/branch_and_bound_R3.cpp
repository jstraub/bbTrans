/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include <vector>
#include "optRot/normal.h"
#include "optRot/node.h"
#include "optRot/lower_bound_R3.h"
#include "optRot/upper_bound_indep_R3.h"
#include "optRot/branch_and_bound.h"

using namespace OptRot;

int main(int argc, char** argv) {

  Eigen::Vector3d t_true(1.,1.,1.);
  std::cout << "true translation: " << t_true.transpose() << std::endl;
  //std::cout << q_true.toRotationMatrix() << std::endl;
  
  Eigen::Vector3d muA1, muA2;
  muA1 << 1.,0.,0.;
  muA2 << 0.,1.,0.;
  Eigen::Vector3d muB1, muB2;
  muB1 = muA1 + t_true;
  muB2 = muA2 + t_true;

  Eigen::Matrix3d Sigma1 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Sigma2 = Eigen::Matrix3d::Identity();
  Sigma1 *= 0.001;
  Sigma2 *= 0.0001;

  std::vector<Normal<3>> gmmA;
  gmmA.push_back(Normal<3>(muA1,Sigma1,0.3));
  gmmA.push_back(Normal<3>(muA2,Sigma2,0.7));

  std::vector<Normal<3>> gmmB;
  gmmB.push_back(Normal<3>(muB1,Sigma1,0.3));
  gmmB.push_back(Normal<3>(muB2,Sigma2,0.7));
  
  Eigen::Vector3d min, max;
  min << 0.,0.,0.;
  max << 2.,2.,2.;
  std::list<NodeR3> nodes;
  nodes.push_back(NodeR3(Box(min, max), 0, std::vector<uint32_t>(1,0)));
  LowerBoundR3 lower_bound(gmmA, gmmB, Eigen::Quaterniond());
  UpperBoundIndepR3 upper_bound(gmmA, gmmB, Eigen::Quaterniond());
//  UpperBoundConvexityLog upper_bound_convexity(gmmA, gmmB);
  
  double eps = 1.e-4;
  uint32_t max_it = 1000;
  BranchAndBound<NodeR3> bb(lower_bound, upper_bound);
  NodeR3 node_star = bb.Compute(nodes, eps, max_it);

  std::cout << "optimum translation: " 
    << node_star.GetBox().GetCenter().transpose()
    << std::endl;
}

