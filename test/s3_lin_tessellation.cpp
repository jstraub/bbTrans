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
  for (const auto& node : nodeAAs) {
    V += node.GetVolume();
  }
  std::cout << "Tessellation AA volume: " << V << std::endl;

//  tetrahedra.clear();
//  for (const auto& node : nodeAAs) {
//    for (uint32_t i=0; i<5; ++i) 
//      tetrahedra.push_back(node.GetNodeS3(i).GetTetrahedron());
//  }
//  std::cout << tetrahedra.size() << std::endl;
//
//  TessellationTest(tetrahedra, 50000);


  std::list<NodeTpS3> nodes = TessellateTpS3();
  std::cout << "initial # nodes: " << nodes.size() << std::endl;

  V = 0.;
  for (const auto& node : nodes) {
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
