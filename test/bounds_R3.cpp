/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "optRot/node.h"
#include "optRot/lower_bound_R3.h"
#include "optRot/upper_bound_indep_R3.h"
#include "optRot/upper_bound_convex_R3.h"

using namespace OptRot;

Eigen::Quaterniond RandomRotation() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> N(0,1);
  Eigen::Quaterniond q(N(gen), N(gen), N(gen), N(gen));
  q.normalize();
  return q;
}

Eigen::Vector3d RandomTranslation() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> N(0,1);
  Eigen::Vector3d t(N(gen), N(gen), N(gen));
  return t*3.;
}

int main(int argc, char ** argv) {
  Eigen::Vector3d t_true;
  t_true << 1.,1.,1.;
  Eigen::Quaterniond q;
  q = RandomRotation();
//  q.setIdentity();
  Eigen::Matrix3d R = q.toRotationMatrix();

  Eigen::Vector3d muA1, muA2;
  muA1 << 1.,0.,0.;
  muA2 << 0.,1.,0.;
  Eigen::Vector3d muB1, muB2;
  muB1 = R*muA1 + t_true;
  muB2 = R*muA2 + t_true;

  Eigen::Matrix3d Sigma1 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Sigma2 = Eigen::Matrix3d::Identity();
  Sigma1 *= 0.1;
  Sigma2 *= 0.01;

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
  min << -2.,-2.,-2.;
  max << 2.,2.,2.;
  std::list<NodeR3> nodes = GenerateNotesThatTessellateR3(min, max, 0.2);

  LowerBoundR3 lower_bound(gmmA, gmmB, q);
  UpperBoundIndepR3 upper_bound(gmmA, gmmB, q);
  UpperBoundConvexR3 upper_bound_convex(gmmA, gmmB, q);
  
  Eigen::VectorXd lbs(nodes.size());
  Eigen::VectorXd ubs(nodes.size());
  Eigen::VectorXd ubCs(nodes.size());
  std::size_t i=0;
  for (NodeR3& node : nodes) {
    std::cout << "---- box " << std::endl;
    std::cout << node.GetBox().GetCenter().transpose() << std::endl;

    std::cout << " ------ " << std::endl;
    lbs[i] = lower_bound.Evaluate(node);
    std::cout << "lower bound: " << lbs[i] << std::endl;
    std::cout << " ------ " << std::endl;
    ubCs[i] = upper_bound_convex.Evaluate(node);
    std::cout << "upper bound C: " << ubCs[i] << std::endl;
    if (ubCs[i] - lbs[i] > 1e3) 
      std::cout << " !! large deviation !!" << std::endl;
    std::cout << " ------ " << std::endl;
    ubs[i] = upper_bound.Evaluate(node);
    std::cout << "upper bound: " << ubs[i] << std::endl;
    ++i;
  }
  std::cout << "Lower Bounds: " << std::endl;
  std::cout << lbs.transpose() << std::endl;
  std::cout << "Upper Bounds: " << std::endl;
  std::cout << ubs.transpose() << std::endl;
  std::cout << "# upper < lower " << (ubs.array() < lbs.array()).transpose()
    << std::endl;
  std::cout << "Upper Bounds Convexity: " << std::endl;
  std::cout << ubCs.transpose() << std::endl;
  std::cout << "# lower < upperC " << (ubCs.array() > lbs.array()).transpose()
    << std::endl;
  std::cout << "# upperC < upper " << (ubCs.array() < ubs.array()).transpose()
    << std::endl;
//  std::cout << "# upperC < lower " << (ubCs.array() < lbs.array()-1e-6).transpose()
//    << std::endl;
//  std::cout << "# upperC - lower " << (ubCs.array() - lbs.array()).transpose()
//    << std::endl;
  std::vector<NodeR3> nodes_v(nodes.begin(), nodes.end());
  std::vector<size_t> idx(nodes.size());
  for (size_t i=0; i<idx.size(); ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&nodes_v,&t_true](size_t i1, size_t i2) {
      return (nodes_v[i1].GetBox().GetCenter()-t_true).norm()
      < (nodes_v[i2].GetBox().GetCenter()-t_true).norm();});

  std::ofstream out("./testBound.csv");
  for (uint32_t i=0; i<lbs.size(); ++i) {
    out << lbs(idx[i]) << " " << ubs(idx[i]) << " " << ubCs(idx[i])
      << " " 
      << (nodes_v[idx[i]].GetBox().GetCenter()-t_true).norm()
      << std::endl;
  }
  out.close();

  return 0;
}
