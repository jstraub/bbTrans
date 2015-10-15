/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <stdint.h>
#include <vector>
#include "optRot/tetrahedron.h"
#include "optRot/s3_tessellation.h"

namespace OptRot {

class Node {
 public:
  Node(const Tetrahedron4D& tetrahedron, uint32_t lvl,
      std::vector<uint32_t> ids);
  ~Node() = default;
  std::vector<Node> Branch() const;
  const Tetrahedron4D& GetTetrahedron() const { return tetrahedron_;}
  uint32_t GetLevel() const {return lvl_;}
  std::vector<uint32_t> GetIds() const {return ids_;}
  double GetUB() const { return ub_;}
  double GetLB() const { return lb_;}
  void SetUB(double ub) { ub_ = ub;}
  void SetLB(double lb) { lb_ = lb;}
  double GetBoundGap() const {return ub_-lb_;}
 private:
  Tetrahedron4D tetrahedron_;
  uint32_t lvl_;
  std::vector<uint32_t> ids_;
  double lb_;
  double ub_;
};

// For use with std::forward_list::remove_if
class IsPrunableNode {
 public:
  IsPrunableNode(double lb) : lb_(lb) {}
  bool operator() (const Node& node) {return node.GetUB() < lb_;}
 private:
  double lb_;
};

struct LessThanNodeUB {
  bool operator() (const Node& node_a, const Node& node_b) { return node_a.GetUB() < node_b.GetUB();}
};

struct LessThanNodeLB {
  bool operator() (const Node& node_a, const Node& node_b) { return node_a.GetLB() < node_b.GetLB();}
};

std::vector<Node> GenerateNotesThatTessellateS3();

}
