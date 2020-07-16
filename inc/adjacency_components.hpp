#ifndef BLUNTIFIER_ADJACENCY_COMPONENTS_HPP
#define BLUNTIFIER_ADJACENCY_COMPONENTS_HPP

/**
 * \file adjacency_components.hpp
 *
 * Defines algorithms for identifying and manipulating components
 * of node side adjacencies.
 */

#include <functional>
#include <random>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <pair>
#include <limits>
#include <cmath>

#include "bdsg/internal/packed_structs.hpp"
#include "handlegraph/handle_graph.hpp"
#include "handlegraph/util.hpp"

namespace bluntifier {

using std::vector;
using std::unordered_set;
using std::function;
using std::pair;
using std::numeric_limits;
using handlegraph::HandleGraph;
using handlegraph::handle_t;
using std::sort;

class AdjacencyComponent {
public:
    
    using bipartition = pair<unordered_set<handle_t>, unordered_set<handle_t>>;
    using const_iterator = vector<handle_t>::const_iterator;
    
    template<typename SideIter>
    AdjacencyComponent(const HandleGraph& graph,
                       SideIter begin, SideIter end);
    
    AdjacencyComponent(AdjacencyComponent&&) = default;
    
    // lambda returns true if iteration should continue. function returns
    // true if iteration was not stopped early by lambda.
    bool for_each_adjacent_side(const handle_t& side,
                                const function<bool(handle_t)>& lambda) const;
    
    const_iterator begin() const;
    
    const_iterator end() const;
    
    size_t size() const;
    
    bool empty() const;
    
    bool is_bipartite() const;
    
    // returns empty sets if the adjacency component is not actually bipartite
    bipartition bipartite_partition() const;
    
    // return the maximum bipartite partition computed in O(Delta * 2^(n-1)) time
    bipartition exhaustive_maximum_bipartite_partition() const;
    
    // return a 1/2-approximation of the maximum bipartite partition
    bipartition maximum_bipartite_partition_apx_1_2(uint64_t seed = 8477176661834875934ull) const;
    
    // greedily modify a bipartite partition until it is locally optimal.
    // uses Bylka, Idzik, & Tuza's (1999) iteration scheme to guarantee Edwards-Erdos
    // bound that #edges across partition is >= 1/2 * m + 1/8 * (sqrt(8 * m + 1) - 1)
    void refine_apx_partition(bipartition& partition,
                              size_t max_opt_steps = numeric_limits<size_t>::max()) const;
    
    // iterate through pairs of subgraphs and bipartitions that have the following properties:
    // - each edge in this adjacency component occurs in exactly one subgraph
    // - each subgraphs's edges are a subset of the parent graph's edges
    // - every subgraph is bipartite with respect to the paired bipartition
    void decompose_into_bipartite_blocks(const function<void(const HandleGraph&,const bipartition&)>& lambda) const;
    
private:
    
    // use a recursively defined Gray code to iterate over bipartitions
    uint64_t recursive_gray_code(bipartition& partition, uint64_t score, size_t index, size_t to_flip,
                                 uint64_t& best_score, uint64_t& best_index)
    
    vector<handle_t> component;
    
    const HandleGraph* graph;
};
    
vector<AdjacencyComponent> adjacency_components(const HandleGraph& graph);

void for_each_adjacency_component(const HandleGraph& graph,
                                  const function<void(AdjacencyComponent&)>& lambda);







/// Template implementations

template<typename SideIter>
AdjacencyComponent::AdjacencyComponent(const HandleGraph& graph,
                                       SideIter begin, SideIter end)
    : graph(&graph), component(begin, end)
{
    // to remove any system dependent behavior due to ordering
    sort(component.begin(), component.end());
}

}

#endif /* BLUNTIFIER_ADJACENCY_COMPONENTS_HPP */
