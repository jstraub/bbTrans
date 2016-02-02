/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <iostream>
#include "bbTrans/combinations.h"

using namespace bb;

int main(int argc, char** argv) {
  
  uint32_t k = 4;
  uint32_t n = 120;
  Combinations combinations(n, k);
  
  std::cout << combinations.Get().size() << "combinations for " << n << "C"
    << k << std::endl;

  for (uint32_t i = 0; i < combinations.Get().size(); ++i) {
    for (uint32_t j = 0; j < k; ++j) {
      std::cout << combinations.Get()[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}
