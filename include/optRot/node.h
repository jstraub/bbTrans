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
 private:
  Tetrahedron4D tetrahedron_;
  uint32_t lvl_;
  std::vector<uint32_t> ids_;
};

std::vector<Node> GenerateNotesThatTessellateS3();

}
