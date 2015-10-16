/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#include "optRot/node.h"

namespace OptRot {

Node::Node(const Tetrahedron4D& tetrahedron, uint32_t lvl,
    std::vector<uint32_t> ids) : tetrahedron_(tetrahedron), lvl_(lvl),
  ids_(ids), lb_(-1e12), ub_(1e12)
{}

std::vector<Node> Node::Branch() const {
  std::vector<Node> nodes;
  nodes.reserve(8);
  std::vector<Tetrahedron4D> tetrahedra = tetrahedron_.Subdivide();
  for (uint32_t i=0; i < tetrahedra.size(); ++i) {
    std::vector<uint32_t> ids(ids_.begin(), ids_.end());
    ids.push_back(i);
    nodes.push_back(Node(tetrahedra[i], lvl_+1, ids));
  }
  return nodes;
}

std::vector<Node> GenerateNotesThatTessellateS3() {
  std::vector<Tetrahedron4D> tetrahedra = TessellateS3();
  std::vector<Node> nodes; 
  nodes.reserve(tetrahedra.size());
  for (uint32_t i=0; i<tetrahedra.size(); ++i) {
    nodes.push_back(Node(tetrahedra[i], 0, std::vector<uint32_t>(1,i)));
  }
  return nodes;
}
}
