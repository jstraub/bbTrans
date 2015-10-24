/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include <vector>
#include "optRot/normal.h"
#include "optRot/node.h"
#include "optRot/lower_bound_R3.h"
#include "optRot/upper_bound_indep_R3.h"
#include "optRot/upper_bound_convex_R3.h"
#include "optRot/branch_and_bound.h"

using namespace OptRot;

int main(int argc, char** argv) {

  Eigen::Vector3d t_true(1.,1.,1.);
  std::cout << "true translation: " << t_true.transpose() << std::endl;
  Eigen::Quaterniond q(1.,1.,1.,1.);
  q.normalize();
  Eigen::Matrix3d R = q.toRotationMatrix();
  std::cout << R << std::endl;
  
  Eigen::Vector3d muA1, muA2;
  muA1 << 1.,0.,0.;
  muA2 << 0.,1.,0.;
  Eigen::Vector3d muB1, muB2;
  muB1 = R*muA1 + t_true;
  muB2 = R*muA2 + t_true;

  Eigen::Matrix3d Sigma1 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Sigma2 = Eigen::Matrix3d::Identity();
  Sigma1 *= 0.01;
  Sigma2 *= 0.001;

  std::vector<Normal<3>> gmmA;
  gmmA.push_back(Normal<3>(muA1,Sigma1,0.3));
  gmmA.push_back(Normal<3>(muA2,Sigma2,0.7));

  std::vector<Normal<3>> gmmB;
  gmmB.push_back(Normal<3>(muB1,R*Sigma1*R.transpose(),0.3));
  gmmB.push_back(Normal<3>(muB2,R*Sigma2*R.transpose(),0.7));

  std::vector<Normal<3>> gmmT;
  ComputeGmmT(gmmA, gmmB, gmmT, q);

  std::cout << "GMM A" << std::endl;
  for (auto& g: gmmA) g.Print();
  std::cout << "GMM B" << std::endl;
  for (auto& g: gmmB) g.Print();
  std::cout << "GMM T" << std::endl;
  for (auto& g: gmmT) g.Print();
  
  Eigen::Vector3d min, max;
  min << 0.,0.,0.;
  max << 2.,2.,2.;
  std::list<NodeR3> nodes = GenerateNotesThatTessellateR3(min, max, 0.3);

  LowerBoundR3 lower_bound(gmmA, gmmB, q);
  UpperBoundIndepR3 upper_bound(gmmA, gmmB, q); 
  UpperBoundConvexR3 upper_bound_convex(gmmA, gmmB, q);
  
  double eps = 1.e-4;
  uint32_t max_it = 1000;
  BranchAndBound<NodeR3> bb1(lower_bound, upper_bound);
  NodeR3 node_star = bb1.Compute(nodes, eps, max_it);

  std::cout << "optimum translation (indep UB): " 
    << node_star.GetBox().GetCenter().transpose()
    << std::endl;

  nodes = GenerateNotesThatTessellateR3(min, max, 0.3);
  BranchAndBound<NodeR3> bb2(lower_bound, upper_bound_convex);
  node_star = bb2.Compute(nodes, eps, max_it);

  std::cout << "optimum translation (convex UB): " 
    << node_star.GetBox().GetCenter().transpose()
    << std::endl;
}

