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

int main(int argc, char ** argv) {
  Eigen::Vector3d min, max;
  min << 0.,0.,0.;
  max << 1.,1.,1.;
  NodeR3 node(Box(min,max), 0, std::vector<uint32_t>(1,0));

  std::cout << " Testing InSide() " << std::endl;
  Eigen::Vector3d p;
  p << 0.3,0.3,0.3;
  if (!node.GetBox().IsInside(p)) std::cout << "- error: " <<
    p.transpose() << std::endl; 
  p << 0.3,0.3,1.;
  if (!node.GetBox().IsInside(p)) std::cout << "- error: " <<
    p.transpose() << std::endl; 
  p << 1.0,0.3,1.;
  if (!node.GetBox().IsInside(p)) std::cout << "- error: " <<
    p.transpose() << std::endl; 
  p << 1.0,0.,1.;
  if (!node.GetBox().IsInside(p)) std::cout << "- error: " <<
    p.transpose() << std::endl; 
  p << 1.1,0.,1.;
  if (node.GetBox().IsInside(p)) std::cout << "- error: " <<
    p.transpose() << std::endl; 
  p << 1.1,999.,1.;
  if (node.GetBox().IsInside(p)) std::cout << "- error: " <<
    p.transpose() << std::endl;
  p << 1.,.5,-1.e-6;
  if (node.GetBox().IsInside(p)) std::cout << "- error: " <<
    p.transpose() << std::endl; 

  std::cout << " Testing FindMinTranslationInNode() " << std::endl;
  // min inside
  Normal<3> g(Eigen::Vector3d(0.5, 0.5, 0.5), Eigen::Matrix3d::Identity(), 1.0);
  Eigen::Vector3d t = FindMinTranslationInNode(g.GetOmega(), g.GetXi(), node);
  if (!((t.array() - g.GetMu().array()).abs() < 1e-6).all())
    std::cout << "- error: " << t.transpose() << std::endl;
  else
    std::cout << "- ok: " << t.transpose() << std::endl;
  t = FindMaxTranslationInNode(g.GetOmega(), g.GetXi(), node);
  std::cout << "- check (should be any of the box corners): " <<
    t.transpose() << std::endl;
  // min on a side
  g = Normal<3>(Eigen::Vector3d(1.5, 0.5, 0.5), Eigen::Matrix3d::Identity(), 1.0);
  t = FindMinTranslationInNode(g.GetOmega(), g.GetXi(), node);
  if (!(((t.array() - Eigen::Vector3d(1.,0.5,0.5).array())).abs() < 1e-6).all())
    std::cout << "- error: " << t.transpose() << std::endl;
  else
    std::cout << "- ok: " << t.transpose() << std::endl;
  t = FindMaxTranslationInNode(g.GetOmega(), g.GetXi(), node);
  if (!fabs(t(0)) < 1e-6)
    std::cout << "- error: " << t.transpose() << std::endl;
  else
    std::cout << "- ok: " << t.transpose() << std::endl;
  // min on an edge
  g = Normal<3>(Eigen::Vector3d(1.5, 1.5, 0.5), Eigen::Matrix3d::Identity(), 1.0);
  t = FindMinTranslationInNode(g.GetOmega(), g.GetXi(), node);
  if (!(((t.array() - Eigen::Vector3d(1.,1.,0.5).array())).abs() < 1e-6).all())
    std::cout << "- error: " << t.transpose() << std::endl;
  else
    std::cout << "- ok: " << t.transpose() << std::endl;
  t = FindMaxTranslationInNode(g.GetOmega(), g.GetXi(), node);
  if (!(fabs(t(0)) < 1e-6 && fabs(t(1)) < 1e-6))
    std::cout << "- error: " << t.transpose() << std::endl;
  else
    std::cout << "- ok: " << t.transpose() << std::endl;
  // min on a corner
  g = Normal<3>(Eigen::Vector3d(1.5, 1.5, 1.5), Eigen::Matrix3d::Identity(), 1.0);
  t = FindMinTranslationInNode(g.GetOmega(), g.GetXi(), node);
  if (!(((t.array() - Eigen::Vector3d(1.,1.,1.).array())).abs() < 1e-6).all())
    std::cout << "- error: " << t.transpose() << std::endl;
  else
    std::cout << "- ok: " << t.transpose() << std::endl;
  t = FindMaxTranslationInNode(g.GetOmega(), g.GetXi(), node);
  if (!(((t.array() - Eigen::Vector3d(0.,0.,0.).array())).abs() < 1e-6).all())
    std::cout << "- error: " << t.transpose() << std::endl;
  else
    std::cout << "- ok: " << t.transpose() << std::endl;
}
