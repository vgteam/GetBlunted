#ifndef BLUNTIFIER_BIPARTITE_GRAPH_HPP
#define BLUNTIFIER_BIPARTITE_GRAPH_HPP

/**
 * \file bipartiteGraph.hpp
 *
 * Defines an interface for a bipartite graph or subgraph
 */

#include "handlegraph/handle_graph.hpp"
#include "handlegraph/types.hpp"

namespace bluntifier {

/*
 * represents the division of the two sides of a bipartite graph
 */
typedef pair<unordered_set<handle_t>, unordered_set<handle_t>> bipartition;
typedef pair<vector<handle_t>, vector<handle_t>> ordered_bipartition;

class BipartiteGraph {
public:
    using const_iterator = vector<handle_t>::const_iterator;
    
    BipartiteGraph(const HandleGraph& graph,
                   const bipartition& partition);
    BipartiteGraph(const HandleGraph& graph,
                   const ordered_bipartition& partition);
    ~BipartiteGraph();
    
    const_iterator left_begin() const;
    const_iterator left_end() const;
    const_iterator left_iterator(const handle_t node) const;
    size_t left_size() const;
    const_iterator right_begin() const;
    const_iterator right_end() const;
    const_iterator right_iterator(const handle_t node) const;
    size_t right_size() const;
    
    // lambda returns true if iteration should continue. function returns
    // true if iteration was not stopped early by lambda.
    bool for_each_adjacent_side(const handle_t& side,
                                const function<bool(handle_t)>& lambda) const;
    
    const ordered_bipartition& bipartition() const;
    
private:
    
    const HandleGraph* graph;
    ordered_bipartition _partition
    unordered_map<handle_t, size_t> left_partition_index;
    unordered_map<handle_t, size_t> right_partition_index;
}


}

#endif /* BLUNTIFIER_BIPARTITE_GRAPH_HPP */
