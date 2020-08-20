#ifndef BLUNTIFIER_BIPARTITE_GRAPH_HPP
#define BLUNTIFIER_BIPARTITE_GRAPH_HPP

/**
 * \file BipartiteGraph.hpp
 *
 * Defines an interface for a bipartite graph or subgraph
 */

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <functional>

#include "handlegraph/handle_graph.hpp"
#include "handlegraph/types.hpp"
#include "handlegraph/util.hpp"

namespace bluntifier {

using std::pair;
using std::vector;
using std::unordered_set;
using std::unordered_map;
using std::function;
using handlegraph::handle_t;
using handlegraph::HandleGraph;


/*
 * Represents the division of the two sides of a bipartite graph
 */
typedef pair<unordered_set<handle_t>, unordered_set<handle_t>> bipartition;
typedef pair<vector<handle_t>, vector<handle_t>> ordered_bipartition;

/*
 * A bipartite subgraph of a HandleGraph
 */
class BipartiteGraph {
public:
    using const_iterator = vector<handle_t>::const_iterator;
    
    // note: constructor does not validate bipartiteness
    BipartiteGraph(const HandleGraph& graph,
                   const bipartition& partition);
    
    ~BipartiteGraph();
    
    size_t get_degree(handle_t node) const;
    const_iterator left_begin() const;
    const_iterator left_end() const;
    const_iterator left_iterator(const handle_t node) const;
    size_t left_size() const;
    const_iterator right_begin() const;
    const_iterator right_end() const;
    const_iterator right_iterator(const handle_t node) const;
    size_t right_size() const;
    bool is_left_side(const handle_t node) const;
    
    // Amilhastre, et al. 1998, algorithm 1. returns a simplified version of this bipartite graph
    // and records all handles along with their successors for each step of the simplification
    // in order
    BipartiteGraph simplify(vector<pair<handle_t, vector<handle_t>>>& simplifications) const;
    
    void for_each_adjacent_side(const handle_t& side,
                                const function<void(handle_t)>& lambda) const;
        
    const HandleGraph& get_graph() const;
private:
    
    BipartiteGraph(const HandleGraph& graph,
                   const ordered_bipartition& partition);
    
    // Amilhastre algorithm
    void simplify_side(const vector<handle_t>& simplifying_partition,
                       const vector<handle_t>& opposite_partition,
                       vector<vector<size_t>>& simplifying_edges,
                       vector<vector<size_t>>& opposite_edges,
                       vector<pair<handle_t, vector<handle_t>>>& simplifications) const;
    
    const HandleGraph* graph;
    ordered_bipartition _partition;
    vector<vector<size_t>> left_edges;
    vector<vector<size_t>> right_edges;
    unordered_map<handle_t, size_t> left_partition_index;
    unordered_map<handle_t, size_t> right_partition_index;
};
}

#endif /* BLUNTIFIER_BIPARTITE_GRAPH_HPP */
