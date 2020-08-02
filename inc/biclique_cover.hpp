#ifndef BLUNTIFIER_BICLIQUE_COVER_HPP
#define BLUNTIFIER_BICLIQUE_COVER_HPP

/**
 * \file biclique_cover.hpp
 *
 * Defines algorithm for computing the biclique cover of a bipartite graph.
 */

#include <vector>
#include <functional>
#include <unordered_map>
#include <unordered_set>

#include "handlegraph/handle_graph.hpp"
#include "handlegraph/types.hpp"
#include "adjacency_components.hpp"
#include "subtractive_graph.hpp"
#include "utility.hpp"

namespace bluntifier {

using std::vector;
using std::unordered_map;
using std::function;
using handlegraph::HandleGraph;

class BicliqueCover {
public:
    
    // initialize with a graph and partition of node sides. the
    // subgraph induced by the partition must be bipartite to be
    // valie (this is not checked)
    BicliqueCover(const HandleGraph& graph,
                  const bipartition& partition);
    ~BicliqueCover();
    
    // compute and return a biclique cover of the partition, where
    // bicliques are each represented by a bipartition of some
    // subset of the nodes
    vector<bipartition> get() const;
    
    
    bool is_domino_free() const;
    
private:
    
    // TODO: redundant with adjacency component
    
    // lambda returns true if iteration should continue. function returns
    // true if iteration was not stopped early by lambda.
    bool for_each_adjacent_side(const handle_t& side,
                                const function<bool(handle_t)>& lambda) const;
    
    void simplify_side(const vector<handle_t>& simplifying_side,
                       SubtractiveHandleGraph& simplifying) const;
    
    
    SubtractiveHandleGraph simplify() const;
    
    vector<bipartition> domino_free_cover(const HandleGraph& simple_graph) const;
    
    vector<bipartition> heuristic_cover() const;
    
    const HandleGraph& graph;
    vector<handle_t> left_partition;
    vector<handle_t> right_partition;
    unordered_map<handle_t, size_t> left_partition_index;
    unordered_map<handle_t, size_t> right_partition_index;
};


}

#endif /* BLUNTIFIER_BICLIQUE_COVER_HPP */
