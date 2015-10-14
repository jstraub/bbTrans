/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include "optRot/combinations.h"

void Combinations(std::vector<uint32_t>& comb_current, uint32_t start,
    uint32_t k_current, uint32_t n, uint32_t k, std::vector<uint32_t>*
    combinations) {
  // If one combination is finished save it to the return list of
  // combinations and return;
  if (k_current >= k) {
    for (uint32_t elem : comb_current) {
      combinations->push_back(elem);
    }
    return;
  }
  for (uint32_t i = start; i < n; ++i) {
    comb_current[k_current] = i;
    Combinations(comb_current, i+1, k_current+1, n, k, combinations);
  }
}
