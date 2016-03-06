/* Copyright (c) 2016, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include "bbTrans/node_TpS3.h"
#include "bbTrans/node_AA.h"
#include "bbTrans/s3_tessellation.h"
using namespace bb;

int main(int argc, char**argv) {

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
  std::vector<Tetrahedron4D> tetrahedra;
  for (const auto& node : nodes) {
    for (uint32_t i=0; i<5; ++i) 
      tetrahedra.push_back(node.GetNodeS3(i).GetTetrahedron());
  }
  std::cout << tetrahedra.size() << std::endl;

  TessellationTest(tetrahedra, 50000);

  p_min = Eigen::Vector3d(-M_PI,-M_PI,-M_PI);
  p_max = Eigen::Vector3d( M_PI, M_PI, M_PI);
  NodeAA rootAA(Box(p_min, p_max),std::vector<uint32_t>(0));
//  std::cout << rootAA.ToString() << std::endl;
  std::vector<NodeAA> l1AA = rootAA.Branch();
  std::list<NodeAA> nodeAAs;
  for (auto& node1 : l1AA) {
    std::vector<NodeAA> l2 = node1.Branch();
//    for (auto& node2 : l2) {
//      std::vector<NodeAA> l3 = node2.Branch();
//      for (auto& node3 : l3) {
//        std::vector<NodeAA> l4 = node3.Branch();
        nodeAAs.insert(nodeAAs.end(), l2.begin(), l2.end());
//      }
//    }
  }
  std::cout << "initial # nodeAAs: " << nodeAAs.size() << std::endl;
  tetrahedra.clear();
  for (const auto& node : nodeAAs) {
    for (uint32_t i=0; i<5; ++i) 
      tetrahedra.push_back(node.GetNodeS3(i).GetTetrahedron());
  }
  std::cout << tetrahedra.size() << std::endl;

  TessellationTest(tetrahedra, 50000);


}
