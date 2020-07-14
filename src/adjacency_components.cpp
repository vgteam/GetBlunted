/**
 * \file adjacency_components.cpp
 *
 * Implements algorithms for identifying and manipulating components
 * of node side adjacencies.
 */

#include "adjacency_components.hpp"

namespace bluntifier {

using bdsg::PackedSet;
using handlegraph::as_integer;
using std::move;
using std::make_pair;
using std::tie;

vector<unordered_set<handle_t>> adjacency_components(const HandleGraph& graph) {
    
    vector<unordered_set<handle_t>> return_val;
    
    for_each_adjacency_component(graph, [&](unordered_set<handle_t>& component) {
        return_val.emplace_back(move(component));
    });
    
    return return_val;
}

void for_each_adjacency_component(const HandleGraph& graph,
                                  const function<void(const unordered_set<handle_t>&)> lambda) {
    PackedSet sides_seen;
    
    graph.for_each_handle([&](const handle_t& handle) {
        for (handle_t side : {handle, graph.flip(handle)}) {
            if (handles_seen.find(as_integer(side))) {
                // this node side hasn't been visited yet, start
                // a new adjacency component
                unordered_set<handle_t> component;
                
                // init the stack and trackers
                component.insert(side);
                sides_seen.insert(side);
                vector<handle_t> stack(1, side);
                while (!stack.empty()) {
                    
                    auto side_here = stack.back();
                    stack.pop_back();
                    
                    graph.follow_edges(side_here, false, [&](const handle_t& neighbor) {
                        handle_t adjacent_side = graph.flip(neighbor);
                        if (!component.count(adjacent_side)) {
                            // we found a new side to add to this component
                            component.insert(adjacent_side);
                            sides_seen.insert(as_integer(adjacent_side));
                            stack.push_back(adjacent_side);
                        }
                    });
                }
                
                lambda(component);
            }
        }
    });
}

bool adjacency_component_is_bipartite(const HandleGraph& graph,
                                      const unordered_set<handle_t>& component) {
    auto partition = bipartite_partition(component);
    return partition.first.size() + partition.second.size() < component.size();
}

pair<unordered_set<handle_t>, unordered_set<handle_t>>
bipartite_partition(const HandleGraph& graph, const unordered_set<handle_t>& component) {
    
    pair<unordered_set<handle_t>, unordered_set<handle_t>> return_val;
    
    if (!component.empty()) {
        
        handle_t start_side = *component.begin();
        
        vector<pair<handle_t, bool>> stack(1, make_pair(start_side, false));
        return_val.first.insert(start_side);
        
        while (!stack.empty()) {
            bool going_left;
            handle_t side_here;
            tie(side_here, going_left) = stack.back();
            stack.pop_back();
            
            auto& partition_across = going_left ? return_val.second : return_val.first;
            auto& partition_here = going_left ? return_val.first : return_val.second;
            
            bool still_bipartite = graph.follow_edges(side_here, false,
                                                      [&](const handle_t& neighbor) {
                handle_t adjacent_side = graph.flip(neighbor);
                if (partition_here.count(adjacent_side)) {
                    // this side was seen with both even and odd parities, the
                    // adjacency component is not bipartite
                    return false;
                }
                else if (!partition_across.count(adjacent_side)) {
                    // add this side to the partition and prepare a search from it
                    // with the opposite parity
                    partition_across.count(adjacent_side);
                    stack.emplace_back(adjacent_side, !going_left);
                }
                return true;
            });
            
            if (!still_bipartite) {
                // prepare to return a sentinel and stop searching
                return_val.first.clear();
                return_val.second.clear();
                break;
            }
        }
        
    }
    
    return return_val;
}

}
