#ifndef BLUNTIFIER_ADJACENCY_COMPONENTS_HPP
#define BLUNTIFIER_ADJACENCY_COMPONENTS_HPP

/**
 * \file adjacency_components.hpp
 *
 * Defines algorithms for identifying and manipulating components
 * of node side adjacencies.
 */

#include <functional>
#include <vector>
#include <unordered_set>
#include <pair>

#include "bdsg/internal/packed_structs.hpp"
#include "handlegraph/handle_graph.hpp"
#include "handlegraph/util.hpp"

namespace bluntifier {

using std::vector;
using std::unordered_set;
using std::function;
using std::pair;
using handlegraph::HandleGraph;
using handlegraph::handle_t;
    
vector<unordered_set<handle_t>> adjacency_components(const HandleGraph& graph);

void for_each_adjacency_component(const HandleGraph& graph,
                                  const function<void(unordered_set<handle_t>&)> lambda);

bool adjacency_component_is_bipartite(const HandleGraph& graph,
                                      const unordered_set<handle_t>& component);

// returns empty sets if the adjacency component is not actually bipartite
pair<unordered_set<handle_t>, unordered_set<handle_t>>
bipartite_partition(const HandleGraph& graph, const unordered_set<handle_t>& component);


}

#endif /* BLUNTIFIER_ADJACENCY_COMPONENTS_HPP */
