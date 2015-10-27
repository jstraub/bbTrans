/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */
namespace OptRot {

template<class Node>
std::vector<uint32_t> CountBranchesInTree(const std::list<Node>& nodes) {
  uint32_t max_lvl = 0;
  for (auto& node : nodes) if (node.GetLevel() > max_lvl) 
    max_lvl = node.GetLevel();
  std::cout << "Max depth in nodes: " << max_lvl << std::endl;
  std::vector<std::map<uint64_t, uint32_t>> counts(max_lvl);
  for (auto& node : nodes) {
    for (uint32_t lvl = 0; lvl < node.GetLevel(); ++lvl) {
      ++ counts[lvl][node.GetIdAtLevel(lvl)];
    }
  }
  std::vector<uint32_t> countsPerLevel(max_lvl,0);
  for (uint32_t lvl = 0; lvl < max_lvl; ++lvl) {
    std::cout << "@ lvl " << lvl << std::endl;
    for (auto& count : counts[lvl]) {
      std::cout << "  id: " << count.first << "\t#:" << count.second << std::endl;
      ++ countsPerLevel[lvl];
    }
  }
  return countsPerLevel;
}

}
