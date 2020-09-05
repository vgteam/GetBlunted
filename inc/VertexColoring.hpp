#ifndef BLUNTIFIER_VERTEX_COLORING_HPP
#define BLUNTIFIER_VERTEX_COLORING_HPP

/**
 * \file VertexColoring.hpp
 *
 * Defines a vertex coloring algorithm
 */

#include <vector>
#include <cstdint>
#include <limits>
#include <utility>
#include <iostream>
#include <cmath>

namespace bluntifier {

using std::vector;

/*
 * Represents an instance of the vertex coloring problem
 */
class VertexColoring {
public:
    
    VertexColoring(const vector<vector<size_t>>& graph);
    VertexColoring() = default;
    ~VertexColoring() = default;
    
    // compute and return a vertex coloring. the paired bool flag indicates whether
    // this vertex coloring has been verified to be a minimum coloring.
    vector<size_t> get(bool& is_exact) const;
    
protected:
    
    
    // use Lawler's (1976) algorithm, implementation is valid for no more than 15 nodes
    vector<size_t> lawlers_algorithm() const;
    
    // use greedy coloring
    vector<size_t> greedy_coloring(const vector<size_t>& order) const;
    
    // Tsukiyama's, et al. (1977) algorithm to list maximal independent
    // sets, returns bitsets to represent each set (only for small graphs)
    // note: requires adjacency lists to be ordered by index
    vector<uint16_t> maximal_independent_sets(const vector<vector<size_t>>& subgraph,
                                              const vector<size_t>& subgraph_trans) const;
    
    // Tsukiyama's, et al. (1977) backtrack procedure
    void maximal_independent_sets_internal(const vector<vector<size_t>>& subgraph,
                                           const vector<size_t>& subgraph_trans,
                                           size_t max_node,
                                           vector<uint16_t>& intersection_size,
                                           vector<vector<uint16_t>>& reset_buckets,
                                           vector<uint16_t>& maximal_ind_sets) const;
    
    
    // TODO: more orderings and heuristic algorithms?
    // TODO: include Mehrota & Trick's (1995) ILP formulation?
    
    vector<size_t> degree_ordering() const;
    
    vector<vector<size_t>> graph;
    
};

}

#endif /* BLUNTIFIER_BICLIQUE_COVER_HPP */
