/* Copyright (c) 2016, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include "bbTrans/node_TpS3.h"
#include "bbTrans/node_AA.h"
#include "bbTrans/s3_tessellation.h"
using namespace bb;

int main(int argc, char**argv) {

  std::list<NodeAA> nodeAAs = TessellateAA();
  std::cout << "initial # nodeAAs: " << nodeAAs.size() << std::endl;

  double V = 0.;
  for (auto& node : nodeAAs) {
    V += node.GetVolume();
  }
  std::cout << "Tessellation AA volume: " << V << std::endl;
  std::cout << "true volume of half sphere: " << M_PI*M_PI << std::endl;

  std::vector<Tetrahedron4D> tetrahedra;
  for (const auto& node : nodeAAs) {
    for (size_t i=0; i<5; ++i) {
      tetrahedra.push_back(node.GetNodeS3(i).GetTetrahedron());
    }
  }
  std::cout << tetrahedra.size() << std::endl;

  uint32_t N = 0;
  uint32_t Nall = 100000;
  for (uint32_t i=0; i<Nall; ++i) {
    S4d q = S4d::Random();
//    if (i%100 == 0) 
//      std::cout << q.vector().transpose() << std::endl;
    for (const auto& tetra : tetrahedra) 
      if (tetra.Intersects(q.vector())) {
        ++N;
//        break;
      }
  }
  std::cout << "fraction all over sphere intersected with tessellation: "
    << static_cast<double>(N)/static_cast<double>(Nall)
    << std::endl;

  N = 0.;
  for (uint32_t i=0; i<Nall; ++i) {
    S4d q = S4d::Random();
    q.vector()(0) = q.vector()(0) < 0. ? -q.vector()(0) : q.vector()(0);
    for (const auto& tetra : tetrahedra) 
      if (tetra.Intersects(q.vector())) {
        ++N;
        break;
      }
  }
  std::cout << "fraction on top half-sphere intersected with tessellation: "
    << static_cast<double>(N)/static_cast<double>(Nall)
    << std::endl;

  N = 0.;
  for (uint32_t i=0; i<Nall; ++i) {
    S4d q = S4d::Random();
    q.vector()(0) = q.vector()(0) < 0. ? q.vector()(0) : -q.vector()(0);
    for (const auto& tetra : tetrahedra) 
      if (tetra.Intersects(q.vector())) {
        ++N;
        break;
      }
  }
  std::cout << "fraction on bottom half-sphere intersected with tessellation: "
    << static_cast<double>(N)/static_cast<double>(Nall)
    << std::endl;

//  TessellationTest(tetrahedra, 50000);

  std::list<NodeTpS3> nodes = TessellateTpS3();
  std::cout << "initial # nodes: " << nodes.size() << std::endl;

  V = 0.;
  for (auto& node : nodes) {
    V += node.GetVolume();
  }
  std::cout << "Tessellation TpS3 volume: " << V << std::endl;
//  std::vector<Tetrahedron4D> tetrahedra;
//  for (const auto& node : nodes) {
//    for (uint32_t i=0; i<5; ++i) 
//      tetrahedra.push_back(node.GetNodeS3(i).GetTetrahedron());
//  }
//  std::cout << tetrahedra.size() << std::endl;
//
//  TessellationTest(tetrahedra, 50000);


}
