/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "optRot/node.h"

namespace OptRot {

BaseNode::BaseNode(uint32_t lvl, std::vector<uint32_t> ids) :
  lvl_(lvl), ids_(ids), lb_(-1e12), ub_(1e12)
{}

BaseNode::BaseNode(uint32_t lvl, std::vector<uint32_t> ids, double lb,
    double ub) : lvl_(lvl), ids_(ids), lb_(lb), ub_(ub) {
}

BaseNode::BaseNode(const BaseNode& node) : lvl_(node.GetLevel()),
  ids_(node.GetIds()), lb_(node.GetLB()), ub_(node.GetUB()) {
}

NodeS3::NodeS3(const Tetrahedron4D& tetrahedron, uint32_t lvl,
    std::vector<uint32_t> ids) : BaseNode(lvl, ids),
  tetrahedron_(tetrahedron) {
}

NodeS3::NodeS3(const NodeS3& node) : BaseNode(node.GetLevel(), node.GetIds(),
    node.GetLB(), node.GetUB()), tetrahedron_(node.GetTetrahedron()) {
}

NodeR3::NodeR3(const Box& box, uint32_t lvl,
    std::vector<uint32_t> ids) : BaseNode(lvl, ids), box_(box) {
}

NodeR3::NodeR3(const NodeR3& node) : BaseNode(node.GetLevel(), node.GetIds(),
    node.GetLB(), node.GetUB()), box_(node.GetBox()) {
}

std::vector<NodeS3> NodeS3::Branch() const {
  std::vector<NodeS3> nodes;
  nodes.reserve(8);
  std::vector<Tetrahedron4D> tetrahedra = tetrahedron_.Subdivide();
  for (uint32_t i=0; i < tetrahedra.size(); ++i) {
    std::vector<uint32_t> ids(this->ids_);
    ids.push_back(i);
    nodes.push_back(NodeS3(tetrahedra[i], lvl_+1, ids));
  }
  return nodes;
}

std::vector<NodeR3> NodeR3::Branch() const {
  std::vector<NodeR3> nodes;
  nodes.reserve(8);
  std::vector<Box> boxs = box_.Subdivide();
  for (uint32_t i=0; i < boxs.size(); ++i) {
    std::vector<uint32_t> ids(this->ids_);
    ids.push_back(i);
    nodes.push_back(NodeR3(boxs[i], lvl_+1, ids));
  }
  return nodes;
}

std::list<NodeS3> GenerateNotesThatTessellateS3() {
  std::vector<Tetrahedron4D> tetrahedra = TessellateS3();
  std::list<NodeS3> nodes; 
//  nodes.reserve(tetrahedra.size());
  for (uint32_t i=0; i<tetrahedra.size(); ++i) {
    nodes.push_back(NodeS3(tetrahedra[i], 0, std::vector<uint32_t>(1,i)));
  }
  return nodes;
}

std::list<NodeR3> GenerateNotesThatTessellateR3(const Eigen::Vector3d&
    min, const Eigen::Vector3d& max, double max_side_len) {
  NodeR3 node0(Box(min, max), 0, std::vector<uint32_t>(1,0));
  std::vector<std::vector<NodeR3>> node_tree;
  node_tree.push_back(std::vector<NodeR3>(1,node0));
  for (uint32_t lvl = 0; lvl < 20; ++lvl) {
    node_tree.push_back(node_tree[lvl][0].Branch());
    for (uint32_t i = 1; i < node_tree[lvl].size(); ++i) {
      std::vector<NodeR3> nodes_new = node_tree[lvl][i].Branch();
      for (auto& node: nodes_new) node_tree[lvl+1].push_back(node);
    }
    std::cout << "@" << lvl+1 << ": # " << node_tree[lvl+1].size() <<  std::endl;
    if (node_tree[lvl+1][0].GetBox().GetSideLengths().maxCoeff() < max_side_len)
      break;
  }
  uint32_t lvl = node_tree.size() -1;
  return std::list<NodeR3>(node_tree[lvl].begin(), node_tree[lvl].end());
}

}
