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
#include "optRot/lower_bound_log.h"
#include "optRot/upper_bound_log.h"
#include "optRot/upper_bound_convexity_log.h"
#include "optRot/branch_and_bound.h"

using namespace OptRot;

Eigen::Quaterniond RandomRotation() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> N(0,1);
  Eigen::Quaterniond q(N(gen), N(gen), N(gen), N(gen));
  q.normalize();
  return q;
}

int main(int argc, char** argv) {

  Eigen::Vector3d t_true(1.,1.,1.);
  std::cout << "true translation: " << t_true.transpose() << std::endl;
  Eigen::Quaterniond q_true(1.,1.,1.,1.);
//  q_true.normalize();
  q_true = RandomRotation();
  std::cout << "true quaternion: " << q_true.coeffs().transpose() << std::endl;
  Eigen::Matrix3d R = q_true.toRotationMatrix();
  std::cout << R << std::endl;
  
  // Setup GMMs
  Eigen::Vector3d mA1, mA2;
  mA1 << 1.,0.,0.;
  mA2 << 0.,1.,0.;
  Eigen::Vector3d mB1, mB2;
  mB1 = R*mA1 + t_true;
  mB2 = R*mA2 + t_true;

  Eigen::Matrix3d Sigma1 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Sigma2 = Eigen::Matrix3d::Identity();
  Sigma1 *= 0.01;
  Sigma2 *= 0.001;

  std::vector<Normal<3>> gmmA;
  gmmA.push_back(Normal<3>(mA1,Sigma1,0.3));
  gmmA.push_back(Normal<3>(mA2,Sigma2,0.7));
  std::vector<Normal<3>> gmmB;
  gmmB.push_back(Normal<3>(mB1,R*Sigma1*R.transpose(),0.3));
  gmmB.push_back(Normal<3>(mB2,R*Sigma2*R.transpose(),0.7));

  // Setup vMF MMs
  Eigen::Vector3d muA1, muA2;
  muA1 << 1.,0.,0.;
  muA2 << 0.,1.,0.;
  Eigen::Vector3d muB1, muB2;
  muB1 = q_true._transformVector(muA1);
  muB2 = q_true._transformVector(muA2);

  std::vector<vMF<3>> vmfs_A;
  vmfs_A.push_back(vMF<3>(muA1,100.,0.3));
  vmfs_A.push_back(vMF<3>(muA2,10.,0.7));
  std::vector<vMF<3>> vmfs_B;
  vmfs_B.push_back(vMF<3>(muB1,100.,0.3));
  vmfs_B.push_back(vMF<3>(muB2,10.,0.7));
  vMFMM<3> vmf_mm_A(vmfs_A);
  vMFMM<3> vmf_mm_B(vmfs_B);

  // Compute rotation.
  std::list<NodeS3> nodes = GenerateNotesThatTessellateS3();
  LowerBoundLog lower_bound(vmf_mm_A, vmf_mm_B);
  UpperBoundLog upper_bound(vmf_mm_A, vmf_mm_B);
  UpperBoundConvexityLog upper_bound_convexity(vmf_mm_A, vmf_mm_B);
  
  double eps = 1.0e-8;
  uint32_t max_it = 1000;
  BranchAndBound<NodeS3> bb(lower_bound, upper_bound);
//  BranchAndBound<NodeS3> bb(lower_bound, upper_bound_convexity);
  NodeS3 node_star = bb.Compute(nodes, eps, max_it);
  Eigen::Quaterniond q = node_star.GetTetrahedron().GetCenterQuaternion();
  q = q.inverse();

  std::cout << "optimum quaternion: " 
    << q.coeffs().transpose()
    << std::endl;

  // Compute translation.
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
  std::list<NodeR3> nodesR3 = GenerateNotesThatTessellateR3(min, max, 0.5);

  LowerBoundR3 lower_boundR3(gmmA, gmmB, q);
  UpperBoundIndepR3 upper_boundR3(gmmA, gmmB, q); 
  UpperBoundConvexR3 upper_bound_convex_R3(gmmA, gmmB, q);
  
  eps = 1.e-30;
  max_it = 5000;
  BranchAndBound<NodeR3> bbR3(lower_boundR3, upper_bound_convex_R3);
  NodeR3 node_starR3 = bbR3.Compute(nodesR3, eps, max_it);
  Eigen::Vector3d t = node_starR3.GetBox().GetCenter();

  std::cout << "optimum translation: " << t.transpose() << std::endl;
  R = q.toRotationMatrix();

  std::cout << "muAs: " << std::endl
    << muA1.transpose() << std::endl
    << muA2.transpose() << std::endl
    << "muBs: " << std::endl
    << muB1.transpose() << std::endl
    << muB2.transpose() << std::endl
    << "muBs after applying transformation" << std::endl
    << (R.transpose()*(muB1)).transpose() << std::endl
    << (R.transpose()*(muB2)).transpose() << std::endl;

  if (((muA1.array() - (R.transpose()*(muB1)).array()).abs() < 1e-6).all()) {
    std::cout << " Rotation is correct according to surface normal distributions!."
      << std::endl;
  } else {
    std::cout << " Rotation is NOT correct according to surface normal distributions!."
      << std::endl;
  }

  std::cout << "mAs: " << std::endl
    << mA1.transpose() << std::endl
    << mA2.transpose() << std::endl
    << "mBs: " << std::endl
    << mB1.transpose() << std::endl
    << mB2.transpose() << std::endl
    << "mBs after applying transformation" << std::endl
    << (R.transpose()*(mB1-t)).transpose() << std::endl
    << (R.transpose()*(mB2-t)).transpose() << std::endl;

  if (((mA1.array() - (R.transpose()*(mB1-t)).array()).abs() < 1e-6).all()) {
    std::cout << " Rotation and translation is correct according to point distributions!."
      << std::endl;
  } else {
    std::cout << " Rotation or translation is NOT correct according to point distributions!."
      << std::endl;
  }
}

