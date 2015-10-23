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
  Eigen::Quaterniond q;
  q.setIdentity();
  //std::cout << q_true.toRotationMatrix() << std::endl;
  
  Eigen::Vector3d muA1, muA2;
  muA1 << 1.,0.,0.;
  muA2 << 0.,1.,0.;
  Eigen::Vector3d muB1, muB2;
  muB1 = muA1 + t_true;
  muB2 = muA2 + t_true;

  Eigen::Matrix3d Sigma1 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Sigma2 = Eigen::Matrix3d::Identity();
  Sigma1 *= 0.01;
  Sigma2 *= 0.001;

  std::vector<Normal<3>> gmmA;
  gmmA.push_back(Normal<3>(muA1,Sigma1,0.3));
  gmmA.push_back(Normal<3>(muA2,Sigma2,0.7));

  std::vector<Normal<3>> gmmB;
  gmmB.push_back(Normal<3>(muB1,Sigma1,0.3));
  gmmB.push_back(Normal<3>(muB2,Sigma2,0.7));
  
  Eigen::Vector3d min, max;
  min << 0.,0.,0.;
  max << 2.,2.,2.;

  NodeR3 node0(Box(min, max), 0, std::vector<uint32_t>(1,0));
  std::vector<std::vector<NodeR3>> node_tree;
  node_tree.push_back(std::vector<NodeR3>(1,node0));
  for (uint32_t lvl = 0; lvl < 4; ++lvl) {
    node_tree.push_back(node_tree[lvl][0].Branch());
    for (uint32_t i = 1; i < node_tree[lvl].size(); ++i) {
      std::vector<NodeR3> nodes_new = node_tree[lvl][i].Branch();
      for (auto& node: nodes_new) node_tree[lvl+1].push_back(node);
    }
    std::cout << "@" << lvl+1 << ": # " << node_tree[lvl+1].size() <<  std::endl;
  }

  uint32_t lvl = 4;
  std::list<NodeR3> nodes(node_tree[lvl].begin(), node_tree[lvl].end());

  LowerBoundR3 lower_bound(gmmA, gmmB, q);
  UpperBoundIndepR3 upper_bound(gmmA, gmmB, q); 
//  UpperBoundConvexityLog upper_bound_convexity(gmmA, gmmB);
  
  double eps = 1.e-4;
  uint32_t max_it = 1000;
  BranchAndBound<NodeR3> bb(lower_bound, upper_bound);
  NodeR3 node_star = bb.Compute(nodes, eps, max_it);

  std::cout << "optimum translation: " 
    << node_star.GetBox().GetCenter().transpose()
    << std::endl;
}

