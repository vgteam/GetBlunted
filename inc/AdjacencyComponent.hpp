#ifndef BLUNTIFIER_ADJACENCY_COMPONENTS_HPP
#define BLUNTIFIER_ADJACENCY_COMPONENTS_HPP

/**
 * \file adjacenceComponent.hpp
 *
 * Defines algorithms for identifying and manipulating components
 * of node side adjacencies.
 */

#include <functional>
#include <random>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <limits>
#include <cmath>

#include "bdsg/internal/packed_structs.hpp"
#include "handlegraph/handle_graph.hpp"
#include "handlegraph/util.hpp"
#include "handlegraph/types.hpp"
#include "SubtractiveHandleGraph.hpp"
#include "BipartiteGraph.hpp"
#include "utility.hpp"

namespace bluntifier {

using std::vector;
using std::unordered_set;
using std::function;
using std::pair;
using std::numeric_limits;
using handlegraph::HandleGraph;
using handlegraph::handle_t;
using std::sort;

// forward declaration
class AdjacencyComponent;

// iterate over the adjacency components of a graph
void for_each_adjacency_component(const HandleGraph& graph,
                                  const function<void(AdjacencyComponent&)>& lambda);

// get a list of all of the adjacency components in a graph
vector<AdjacencyComponent> adjacency_components(const HandleGraph& graph);

/*
 * Class that represents a collection of node sides that are connected by
 * edges (without crossing any nodes)
 */
class AdjacencyComponent {
public:
    
    // iterator over the node sides of the component
    using const_iterator = vector<handle_t>::const_iterator;
    
    // Initialize using a container of node sides (represented by handle_t's)
    // Does not check that the sides provided constitute an adjacency component.
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
    
    // iterate through pairs of subgraphs and bipartitions that have the following properties:
    // - each edge in this adjacency component occurs in exactly one subgraph
    // - every subgraph is bipartite
    void decompose_into_bipartite_blocks(const function<void(const BipartiteGraph&)>& lambda) const;
    
    
    
    bool is_bipartite() const;
    
    // returns empty sets if the adjacency component is not actually bipartite
    bipartition bipartite_partition() const;
    
    // return the maximum bipartite partition computed in O(Delta * 2^(n-1)) time
    bipartition exhaustive_maximum_bipartite_partition() const;
    
    // return a expected 1/2-approximation of the maximum bipartite partition
    bipartition maximum_bipartite_partition_apx_1_2(uint64_t seed = 8477176661834875934ull) const;
    
    // greedily modify a bipartite partition until it is locally optimal.
    // uses Bylka, Idzik, & Tuza's (1999) iteration scheme to guarantee Edwards-Erdos
    // bound that #edges across partition is >= 1/2 * m + 1/8 * (sqrt(8 * m + 1) - 1)
    void refine_apx_partition(bipartition& partition,
                              size_t max_opt_steps = numeric_limits<size_t>::max()) const;
    
    // TODO: include Goemans-Williamson SDP algorithm with .88 approx ratio?
    
private:
    
    // use a recursively defined Gray code to iterate over all bipartitions
    uint64_t recursive_gray_code(bipartition& partition, uint64_t score, size_t index, size_t to_flip,
                                 uint64_t& best_score, uint64_t& best_index) const;
    
    vector<handle_t> component;
    
    const HandleGraph* graph;
};







/// Template implementations

template<typename SideIter>
AdjacencyComponent::AdjacencyComponent(const HandleGraph& graph,
                                       SideIter begin, SideIter end)
    : component(begin, end), graph(&graph)
{
    // to remove any system dependent behavior due to ordering
    sort(component.begin(), component.end());
}

}

#endif /* BLUNTIFIER_ADJACENCY_COMPONENTS_HPP */
