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
#include <algorithm>

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

/*
 * Represents an instance of the minimum biclique cover problem for a
 * bipartite subgraph of a larger graph
 */
class BicliqueCover {
public:
    
    // initialize with a graph and partition of node sides. the
    // subgraph induced by the partition must be bipartite to be
    // valid (this is not checked)
    BicliqueCover(const HandleGraph& graph,
                  const bipartition& partition);
    ~BicliqueCover();
    
    // compute and return a biclique cover of the partition, where
    // bicliques are each represented by a bipartition of some
    // subset of the nodes
    vector<bipartition> get() const;
    
private:
    
    struct GaloisLattice {
        GaloisLattice();
        ~GaloisLattice();
        
        void separator();
    }
    
    // lambda returns true if iteration should continue. function returns
    // true if iteration was not stopped early by lambda.
    // TODO: redundant with adjacency component
    bool for_each_adjacent_side(const handle_t& side,
                                const function<bool(handle_t)>& lambda) const;
    
    void simplify_side(const vector<handle_t>& simplifying_side,
                       SubtractiveHandleGraph& simplifying) const;
    
    // Amilhastre, et al. 1998, algorithm 2. remove edges without affecting
    // the biclique cover
    SubtractiveHandleGraph simplify() const;
    
    // apply local checking (Amilhastre, et al. 1998, algorithm 3) to all
    // nodes in order to check if the bipartite graph is domino-free
    bool is_domino_free() const;
    
    vector<bipartition> domino_free_cover() const;
    
    vector<bipartition> heuristic_cover() const;
    
    GaloisLattice get_galois_lattice(const HandleGraph& simple_graph) const;
    
    const HandleGraph& graph;
    vector<handle_t> left_partition;
    vector<handle_t> right_partition;
    unordered_map<handle_t, size_t> left_partition_index;
    unordered_map<handle_t, size_t> right_partition_index;
    
    friend class QuotientBall;
};

/*
 * Represents the quotient graph of two-hop subgraph starting at a center
 * node over the equivalence relationship of having the same neighborhood
 * (applied to nodes on the same side as the center)
 */
class QuotientBall {
public:
    QuotientBall(const BicliqueCover& parent, handle_t center);
    QuotientBall() = delete;
    ~QuotientBall() = default;
    
    // Amilhastre, et al (1998) algorithm 3
    bool check_neighbor_ordering_property() const;
    
private:
    vector<vector<size_t>> left_edges, right_edges;
};


}

#endif /* BLUNTIFIER_BICLIQUE_COVER_HPP */
