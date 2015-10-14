/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

#include <stdint.h>
#include <vector>

/// Compute all combinations of n numbers when taking k at a time
/// (nCk).  Inspired by
/// http://www.cs.utexas.edu/users/djimenez/utsa/cs3343/lecture25.html 
void Combinations(std::vector<uint32_t>& comb_current, uint32_t start,
    uint32_t k_current, uint32_t n, uint32_t k, std::vector<uint32_t>*
    combinations);
