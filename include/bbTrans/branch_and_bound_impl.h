/* Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu> Licensed
 * under the MIT license. See the license file LICENSE.
 */

namespace bb {

template <class Node>
BranchAndBound<Node>::BranchAndBound(Bound<Node>& lower_bound,
    Bound<Node>& upper_bound)
  : lower_bound_(lower_bound), upper_bound_(upper_bound)
{}

template <class Node>
uint32_t BranchAndBound<Node>::BoundAndPrune(std::list<Node>& nodes, double& lb,
    double& ub, double eps) {
//  lb = -1.; 
  ub = -1.; // only reset UB but LB should only be able to go up
  for (auto& node : nodes) {
    // Because of numerics...
//    if (node.GetUB() < node.GetLB()) {
//      std::cout << " ub < lb: " << node.GetUB() << " < " <<
//        node.GetLB() << " - " << node.GetUB() - node.GetLB()
//        << " lvl " << node.GetLevel()
//        << std::endl;
//      node.SetUB(node.GetLB()+eps); 
//    }
    lb = std::max(lb, node.GetLB());
    ub = std::max(ub, node.GetUB());
  }
  // Prune
  nodes.remove_if(IsPrunableNode<Node>(lb));
  return std::distance(nodes.begin(), nodes.end());
}

template <class Node>
void BranchAndBound<Node>::WriteStats(std::ofstream& out,
    std::list<Node>& nodes, double lb, double ub, double dt,
    typename
    std::list<Node>::iterator& node_star) {
    double V = 0.; // Volume that the nodes explain
    for (auto& node : nodes) V += node.GetVolume(); 
    uint32_t n_nodes = std::distance(nodes.begin(), nodes.end());
    out << lb << " " << ub << " " << n_nodes << " " << V 
      << " " << dt
      << " " << node_star->GetLB()
      << " " << node_star->GetUB()
      << " " << node_star->GetLevel() << std::endl;
};

template <class Node>
void BranchAndBound<Node>::WriteNodes(std::ofstream& out, std::list<Node>& nodes, double lb, double ub) {
    double V = 0.; // Volume that the nodes explain
    for (auto& node : nodes) V += node.GetVolume(); 
    uint32_t n_nodes = std::distance(nodes.begin(), nodes.end());
    out << lb << " " << ub << " " << n_nodes << " " << V << std::endl;
    for (auto& node : nodes) {
      out << node.GetLB() << " " << node.GetUB() << " " << node.GetLevel()
        << " " << node.GetVolume() << std::endl;
      out << node.Serialize();
    }
};

template <class Node>
Node BranchAndBound<Node>::Compute(std::list<Node>& nodes, double eps,
    uint32_t max_lvl, uint32_t max_it, double lb0, double ub0) {

  // Prepare output
  bool write_stats = false;
  std::ofstream out, outNodes;
  if (write_stats) {
    std::stringstream ss;
    ss << "./bb_iteration_stats_" << nodes.begin()->GetSpace() << ".csv";
    out.open(ss.str());
    std::stringstream ssNodes;
    ssNodes << "./bb_nodes_per_iteration_" << nodes.begin()->GetSpace() << ".csv";
    outNodes.open(ssNodes.str());
  }
  jsc::Timer t0;
  // Get initial upper and lower bounds
  double lb = lb0;
  double ub = ub0;
  Node node0 = nodes.front();
  // ugly but works
  std::vector<Node*> parForNodes(nodes.size(), nullptr);
  size_t i=0;
  for (auto& node : nodes) 
    parForNodes[i++] = &node;
#pragma omp parallel for
  for (size_t i=0; i<parForNodes.size(); ++i) {
    lower_bound_.EvaluateAndSet(*parForNodes[i]);
    upper_bound_.EvaluateAndSet(*parForNodes[i]);
    // Because of numerics in S3 case...
//    if (node.GetUB() < node.GetLB()) {
//      node.SetUB(node.GetLB()+eps); 
//    }
//    }
  }
//  for (const auto& node : nodes) 
//    std::cout << node.GetLB() << " " << node.GetUB() << std::endl;
  typename std::list<Node>::iterator node_star = FindBestNode(nodes, eps);
  if (write_stats) WriteStats(out, nodes, lb, ub, t0.toc(), node_star);
  if (write_stats) WriteNodes(outNodes, nodes, lb, ub);
  // Prune
  nodes.remove_if(IsPrunableNode<Node>(lb));
  uint32_t n_nodes  = std::distance(nodes.begin(), nodes.end());
  uint32_t it = 0;
  do  {
    t0.tic();
    // Find node with the biggest upper bound (the most promising node)
    // to explore further.
    auto node_i = FindBestNodeToExplore(nodes, eps);
    // Find the current best estimate as the maximizer of the lower
    // bounds
    node_star = FindBestNode(nodes, eps);
    if (node_star->GetLevel() >= max_lvl) break;
    // Find the node with the biggest lower bound (the most
    // conservative node to return).
    if (it%(max_it/10) == 0)
      std::cout << "@" << it << " # " << n_nodes 
        << ": exploring " << node_i->GetLB() << " < " << node_i->GetUB() 
        << " lvl " << node_i->GetLevel()
        << ": best " << node_star->GetLB() << " < " << node_star->GetUB() 
        << " lvl " << node_star->GetLevel()
        << "\t global "
        << lb << " < " << ub << " " << " |.| " << (ub - lb)/lb
        << std::endl;
      // Branch and check resulting nodes.
      std::vector<Node> new_nodes = node_i->Branch();
//      for (auto& node : new_nodes) {
#pragma omp parallel for
      for (size_t i=0; i<new_nodes.size(); ++i) {
        upper_bound_.EvaluateAndSet(new_nodes[i]);
        lower_bound_.EvaluateAndSet(new_nodes[i]);
        // Because of numerics in S3 case...
//        if (node.GetUB() < node.GetLB()) {
//          if (node.GetLevel() < 8) 
//            std::cout << " ub < lb: " << node.GetUB() << " < " <<
//              node.GetLB() << " - " << node.GetUB() - node.GetLB() 
//              << " lvl " << node.GetLevel() << "@" << it << " # " <<
//              n_nodes << ": cur " << node_i->GetLB() << " < " <<
//              node_i->GetUB() << " lvl " << node_i->GetLevel() 
//              << "\t global LB " << lb << " UB " << ub << " " << " |.| " <<
//              (ub - lb)/lb << std::endl;
//          node.SetUB(node.GetLB()+10*eps); 
//        }
      }
      for (auto& node : new_nodes) {
        if (node.GetUB() > lb) {
          // Remember this node since we cannot prune it.
          nodes.push_back(node);
          ++n_nodes;
        }
      }
    // Pop the examined node since we are branching here
    nodes.erase(node_i); 
    n_nodes = BoundAndPrune(nodes, lb, ub, eps);
    if (write_stats) WriteStats(out, nodes, lb, ub, t0.toc(), node_star);
    if (write_stats) WriteNodes(outNodes, nodes, lb, ub);

//    for (auto& node : new_nodes) 
//      std::cout << " (" << node.GetLB() << " " << node.GetUB() << ") ";
//    std::cout << std::endl;
  } while (it++ < max_it && (ub - lb)/lb > eps && n_nodes >= 1);

  if (write_stats) out.close();
  if (write_stats) outNodes.close();

  if (n_nodes > 0) {
    node_star = FindBestNode(nodes, eps);
    std::cout << "@" << it << " # " << n_nodes << ": global " 
      << lb << " < " << ub << " " << " |.| " << fabs(ub - lb)/lb << "\t selected "
      << node_star->GetLB() << " < " << node_star->GetUB() << " lvl "
      << node_star->GetLevel() << " " //<< node_star.GetIds()[0] << "\t" 
      << std::endl;
    return *node_star;
  } else {
    std::cout << "@" << it << " pruned all nodes " 
      << lb << " < " << ub << " " << " |.| " << fabs(ub - lb)/lb << std::endl;
    node0.SetLB(-1);
    node0.SetUB(-1);
    return node0;
  }

}

template<class Node>
typename std::list<Node>::iterator BranchAndBound<Node>::FindBestNodeToExplore(std::list<Node>& nodes, double eps) {
  auto node_i = std::max_element(nodes.begin(), nodes.end(),
      LessThanNodeUB<Node>());
  uint32_t Neq = 0;
  auto node_istar = node_i;
  for (auto it=nodes.begin(); it!=nodes.end(); it++) 
    if (node_i->GetUB() - it->GetUB() < eps 
        && it->GetLB() > node_istar->GetLB()) {
      ++ Neq;
      node_istar = it;
      std::cout << "angle between equal nodes: " 
        << node_i->DistanceTo(*node_istar)
        << std::endl;
    }
  if (Neq > 1) {
    std::cout << " equal to explore node updates: " << Neq << std::endl;
  }
  return node_istar;
}

template<class Node>
typename std::list<Node>::iterator BranchAndBound<Node>::FindBestNode(std::list<Node>& nodes, double eps) {
  auto node_star = std::max_element(nodes.begin(), nodes.end(),
      LessThanNodeLB<Node>());
  uint32_t Neq = 0;
  auto node_starstar = node_star;
  for (auto it=nodes.begin(); it!=nodes.end(); it++) 
    if (node_star->GetLB() - it->GetLB() < eps && it->GetUB() < node_starstar->GetUB()) {
      node_starstar = it;
      ++ Neq;
    }
  if (Neq > 1) {
    std::cout << " equal best nodes: " << Neq << std::endl;
  }
  return node_starstar;
}

}
