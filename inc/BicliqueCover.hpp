#ifndef BLUNTIFIER_BICLIQUE_COVER_HPP
#define BLUNTIFIER_BICLIQUE_COVER_HPP

/**
 * \file BicliqueCover.hpp
 *
 * Defines algorithm for computing the biclique cover of a bipartite graph.
 */

#include <vector>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <deque>
#include <queue>
#include <tuple>


#include "handlegraph/handle_graph.hpp"
#include "handlegraph/types.hpp"
#include "BipartiteGraph.hpp"
#include "GaloisLattice.hpp"
#include "ReducedDualGraph.hpp"

namespace bluntifier {

using std::vector;
using std::unordered_map;
using std::function;
using handlegraph::HandleGraph;

/*
 * Represents an instance of the minimum biclique cover problem for a
 * bipartite subgraph of a larger graph
 */
class BicliqueCover {
public:
    
    // initialize with a graph and partition of node sides. the
    // subgraph induced by the partition must be bipartite to be
    // valid (this is not checked)
    BicliqueCover(const BipartiteGraph& graph);
    ~BicliqueCover();
    
    // compute and return a biclique cover of the partition, where
    // bicliques are each represented by a bipartition of some
    // subset of the nodes
    vector<bipartition> get() const;
    
    
    
    
    
    // use Ene's, et al. (2008) fast heuristic
    vector<bipartition> biclique_cover_apx() const;
    
    // use Ene's, et al. (2008) lattice-based post-processing for an
    // approximate solution
    void lattice_polish(vector<bipartition>& biclique_cover) const;
    
private:
    
    // convert biclique cover for a simplified graph into cover
    // for the original graph (using procecure described in proof
    // of property 4.2 in Amilhastre, et al. 1998)
    void unsimplify(vector<bipartition>& simplified_cover,
                    const vector<pair<handle_t, vector<handle_t>>>& simplifications) const;
    
    const BipartiteGraph& graph;
    
};

}

#endif /* BLUNTIFIER_BICLIQUE_COVER_HPP */
