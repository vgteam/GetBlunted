#ifndef BLUNTIFIER_BICLIQUE_COVER_HPP
#define BLUNTIFIER_BICLIQUE_COVER_HPP

/**
 * \file biclique_cover.hpp
 *
 * Defines algorithm for computing the biclique cover of a bipartite graph.
 */

#include <vector>
#include <functional>

#include "handlegraph/handle_graph.hpp"
#include "handlegraph/types.hpp"
#include "adjacency_components.hpp"
#include "subtractive_graph.hpp"
#include "utility.hpp"

namespace bluntifier {

using std::vector;
using std::function;
using handlegraph::HandleGraph;

class BicliqueCover {
public:
    
    BicliqueCover(const HandleGraph& graph,
                  const bipartition& partition);
    ~BicliqueCover();
    
    // compute and return a biclique cover of the graph, where
    // bicliques are each represented by a bipartition of some
    // subset of the nodes
    vector<bipartition> get() const;
    
    
    bool is_domino_free() const;
    
    SubtractiveHandleGraph simplify() const;
    
    vector<bipartition> domino_free_cover() const;
    
    vector<bipartition> heuristic_cover() const;
    
private:
    
    // TODO: redundant with adjacency component
    
    // lambda returns true if iteration should continue. function returns
    // true if iteration was not stopped early by lambda.
    bool for_each_adjacent_side(const handle_t& side,
                                const function<bool(handle_t)>& lambda) const;
    
    const HandleGraph& graph;
    const bipartition& partition;
};


}

#endif /* BLUNTIFIER_BICLIQUE_COVER_HPP */
