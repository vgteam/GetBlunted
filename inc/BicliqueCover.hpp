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

#include "handlegraph/handle_graph.hpp"
#include "handlegraph/types.hpp"
#include "AdjacencyComponent.hpp"
#include "SubtractiveHandleGraph.hpp"
#include "BipartiteGraph.hpp"
#include "utility.hpp"

namespace bluntifier {

using std::vector;
using std::unordered_map;
using std::function;
using handlegraph::HandleGraph;

class GaloisLattice;

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
    
    // lambda returns true if iteration should continue. function returns
    // true if iteration was not stopped early by lambda.
    // TODO: redundant with adjacency component
    bool for_each_adjacent_side(const handle_t& side,
                                const function<bool(handle_t)>& lambda) const;
    
private:
    
    
    void simplify_side(const vector<handle_t>& simplifying_side,
                       SubtractiveHandleGraph& simplifying) const;
    
    // Amilhastre, et al. 1998, algorithm 2. remove edges without affecting
    // the biclique cover
    SubtractiveHandleGraph simplify() const;
    
    vector<bipartition> domino_free_cover() const;
    
    vector<bipartition> heuristic_cover() const;
        
    const BipartiteGraph& graph;
    
};


/*
 * Represents the quotient graph of two-hop subgraph starting at a center
 * node over the equivalence relationship of having the same neighborhood
 * (applied to nodes on the same side as the center)
 */
class CenteredGaloisTree {
public:
    
    class edge_iterator;
    
    CenteredGaloisTree(const BipartiteGraph& cover, handle_t center);
    CenteredGaloisTree() = delete;
    ~CenteredGaloisTree() = default;
    
    // Amilhastre, et al (1998) algorithm 3
    bool has_neighbor_ordering_property() const;
    
    // returns the number of maximal bicliques
    size_t size() const;
    
    // the immediate predecessors in the Hasse diagram of the maximal
    // bicliques
    const vector<size_t>& predecessors(size_t i) const;
    
    // the successor of a biclique in the Hasse diagram, or
    // numeric_limits<size_t>::max() if there is none
    size_t successor(size_t i) const;
    
    // get the biclique that corresponds to the root of the Galois tree
    size_t central_equivalence_class() const;
    
    // returns the size of the right partition part of the i-th biclique
    size_t right_size(size_t i) const;
    
    // get the i-th biclique in this tree
    bipartition biclique(size_t i) const;
    
    // iterate over the edges of the i-th equivalence class (note: not biclique)
    edge_iterator edge_begin(size_t i) const;
    edge_iterator edge_end(size_t i) const;
    
    class edge_iterator {
    public:
        
        edge_iterator& operator=(const edge_iterator& other) = default;
        edge_iterator& operator++();
        pair<handle_t, handle_t> operator*() const;
        bool operator==(const edge_iterator& other) const;
        bool operator!=(const edge_iterator& other) const;
    private:
        edge_iterator() = delete;
        edge_iterator(size_t left, size_t right, size_t eq_class,
                      const CenteredGaloisTree* iteratee);
        size_t left, right, eq_class;
        const CenteredGaloisTree* iteratee;
        friend class CenteredGaloisTree;
    };
    
private:
    
    // clear internal structures to indicate that the N.O.P. doesn't hold
    void clear();
    
    // the equivalanece classes of right nodes that have the same neighborhoods
    vector<vector<handle_t>> equiv_classes;
    // the neighborhood of each equivalence class
    vector<vector<handle_t>> neighborhoods;
    // the immediate successor of the corresponding biclique for each
    // equivalence class
    vector<size_t> successors;
    // the immediate predecessors of the corresponding biclique for each
    // equivalence class
    vector<vector<size_t>> equiv_class_predecessors;
    
    friend class edge_iterator;
};

/*
 * The lattice representing the containment ordering of
 */
class GaloisLattice {
public:
    // construct with algorithm 4 and 5 from Amilhastre, et al. (1998)
    GaloisLattice(const BipartiteGraph& graph);
    ~GaloisLattice();
    
    bool is_domino_free() const;
    
    vector<bipartition> biclique_cover() const;
    
private:
    
    void clear();
    
    //void separator();
    
    vector<CenteredGaloisTree> galois_trees;
    unordered_map<pair<size_t, size_t>, vector<pair<size_t, size_t>>> lattice;
};


}

#endif /* BLUNTIFIER_BICLIQUE_COVER_HPP */
