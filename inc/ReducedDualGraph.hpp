#ifndef BLUNTIFIER_DUAL_GRAPH_HPP
#define BLUNTIFIER_DUAL_GRAPH_HPP

/**
 * \file ReducedDualGraph.hpp
 *
 * Defines an algorithm to compute biclique cover using a dual graph
 */

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <algorithm>
#include <bitset>
#include <cmath>

#include "handlegraph/types.hpp"
#include "BipartiteGraph.hpp"
#include "VertexColoring.hpp"
#include "utility.hpp"

namespace bluntifier {

using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::function;
using std::pair;

/*
 * Implicitly represents the dual graph described by Ene, et al. (2008)
 */
class ReducedDualGraph {
public:
    
    ReducedDualGraph(const BipartiteGraph& graph);
    ReducedDualGraph() = default;
    ~ReducedDualGraph() = default;
    
    // computes the biclique cover, after which this object
    // can no longer be meaningfully queried
    vector<bipartition> biclique_cover(bool& is_exact_out);
    
protected:
    
    
    
    // the number of dual nodes left the reduced graph
    size_t reduced_size() const;
    
    // apply reductions until convergence
    void reduce();
    
    // remove a node from the dual graph and maintain the implicit
    // structures`
    void do_reduction(size_t i, size_t dominatee);
    
    // execute a function on the neighborhood of dual node i, with i
    // included as well. requires a bank of bools for the left, right
    // and dual nodes. while executing the function, the dual neighbor
    // bank will be accurate and complete. the left/right neighbor bank
    // may not be accurate complete. all banks should be all false
    // going into this method, and they will be reset before leaving it.
    void dual_neighborhood_do(size_t i, vector<bool>& is_left_neighbor,
                              vector<bool>& is_right_neighbor,
                              vector<bool>& is_dual_neighbor,
                              const function<void(const vector<size_t>&)>& lambda);
    
    vector<size_t> reduced_clique_partition(bool& is_exact_out);
    
    // undo the dual graph reduction operation while adding the dual nodes
    // into the partition
    // note: dual graph is no longer queryable after this
    void reverse_reduction(vector<size_t>& reduced_partition);
    
    // convert a clique partition of the dual graph into a biclique cover
    vector<bipartition> convert_to_biclique_cover(vector<size_t>& clique_partition) const;
    
    // apply Ene's, et al. (2008) postprocessing heuristic
    void lattice_polish(vector<bipartition>& cover) const;
    
    vector<vector<size_t>> left_edges;
    vector<unordered_map<size_t, size_t>> left_edge_index;
    vector<vector<size_t>> right_edges;
    vector<unordered_map<size_t, size_t>> right_edge_index;
    
    vector<pair<size_t, size_t>> dual_nodes;
    unordered_map<pair<size_t, size_t>, size_t> edge_to_dual_node;
    
    // records of (reduced index, dominatee index/max (if no dominatee));
    vector<pair<size_t, size_t>> reductions;
    
    const BipartiteGraph* graph;
};

}

#endif /* BLUNTIFIER_DUAL_GRAPH_HPP */
