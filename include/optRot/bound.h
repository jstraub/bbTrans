/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
#pragma once

#include "optRot/node.h"

namespace OptRot {

class Bound {
 public:
  Bound() = default;
  virtual ~Bound() = default;
  virtual double Evaluate(const Node& node) = 0; 
 private:
};
}
