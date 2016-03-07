/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include <stdint.h>
#include <vector>
#include <list>
#include <sstream>
#include <string>

#include "manifold/S.h"
#include "bbTrans/node_Lin.h"
#include "bbTrans/node_R3.h"
#include "bbTrans/node_S3.h"
#include "bbTrans/box.h"
#include "bbTrans/tetrahedron.h"

namespace bb {

class NodeAA : public NodeLin {
 public:
  NodeAA(const Box& box, std::vector<uint32_t> ids);
  NodeAA(const NodeAA& node);
  virtual ~NodeAA() = default;
  virtual std::vector<NodeAA> Branch() const;

  std::string GetSpace() const { return "AA"; }
 protected:
  virtual Tetrahedron4D TetraFromBox(const Box& box, uint32_t i0, uint32_t i1,
    uint32_t i2, uint32_t i3);
};
std::list<NodeAA> TessellateAA();
}
