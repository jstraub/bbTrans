/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include "optRot/combinations.h"

using namespace OptRot;

int main(int argc, char** argv) {
  
  uint32_t k = 4;
  uint32_t n = 120;
  std::vector<uint32_t> comb_work(k);
  std::vector<uint32_t> combinations;

  Combinations(comb_work, 0, 0, n, k, &combinations);
  
  std::cout << combinations.size() << "combinations for " << n << "C"
    << k << std::endl;

  for (uint32_t i = 0; i < combinations.size()/k; ++i) {
    for (uint32_t j = 0; j < k; ++j) {
      std::cout << combinations[i*k+j] << "\t";
    }
    std::cout << std::endl;
  }
}
