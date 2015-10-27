/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <stdint.h>
#include <vector>
#include <list>
#include "optRot/node.h"
#include "optRot/tetrahedron.h"
#include "optRot/s3_tessellation.h"

namespace OptRot {

class NodeS3 : public BaseNode {
 public:
  NodeS3(const Tetrahedron4D& tetrahedron, uint32_t lvl,
      std::vector<uint32_t> ids);
  NodeS3(const NodeS3& node);
  virtual ~NodeS3() = default;
  virtual std::vector<NodeS3> Branch() const;
  const Tetrahedron4D& GetTetrahedron() const { return tetrahedron_;}
 protected:
  Tetrahedron4D tetrahedron_;
};

std::list<NodeS3> GenerateNotesThatTessellateS3();
}
