/**
 * \file adjacency_components.cpp
 *
 * Implements algorithms for identifying and manipulating components
 * of node side adjacencies.
 */
#include <iostream>
#include "adjacency_components.hpp"

//#define debug_is_bipartite
//#define debug_refinement

namespace bluntifier {

using bdsg::PackedSet;
using handlegraph::as_integer;
using handlegraph::edge_t;
using std::move;
using std::make_pair;
using std::tie;
using std::linear_congruential_engine;
using std::uniform_int_distribution;
using std::cerr;
using std::endl;

bool AdjacencyComponent::for_each_adjacent_side(const handle_t& side,
                                           const function<bool(handle_t)>& lambda) const {
    return graph->follow_edges(side, false, [&](const handle_t& neighbor) {
        return lambda(graph->flip(neighbor));
    });
}


AdjacencyComponent::const_iterator AdjacencyComponent::begin() const {
    return component.begin();
}

AdjacencyComponent::const_iterator AdjacencyComponent::end() const {
    return component.end();
}

size_t AdjacencyComponent::size() const {
    return component.size();
}

bool AdjacencyComponent::empty() const {
    return component.empty();
}

vector<AdjacencyComponent> adjacency_components(const HandleGraph& graph) {
    
    vector<AdjacencyComponent> return_val;
    
    for_each_adjacency_component(graph, [&](AdjacencyComponent& component) {
        return_val.emplace_back(move(component));
    });
    
    return return_val;
}

void for_each_adjacency_component(const HandleGraph& graph,
                                  const function<void(AdjacencyComponent&)>& lambda) {
    PackedSet sides_seen;
    
    graph.for_each_handle([&](const handle_t& handle) {
        for (handle_t side : {handle, graph.flip(handle)}) {
            if (!sides_seen.find(as_integer(side))) {
                // this node side hasn't been visited yet, start
                // a new adjacency component
                unordered_set<handle_t> component;
                
                // init the stack and trackers
                component.insert(side);
                sides_seen.insert(as_integer(side));
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
                
                AdjacencyComponent adj_component(graph, component.begin(), component.end());
                lambda(adj_component);
            }
        }
    });
}

bool AdjacencyComponent::is_bipartite() const {
    auto partition = bipartite_partition();
    return partition.first.size() + partition.second.size() == component.size();
}

AdjacencyComponent::bipartition AdjacencyComponent::bipartite_partition() const {
    
    bipartition return_val;
    
    if (!component.empty()) {
        
        handle_t start_side = *component.begin();
        
        vector<pair<handle_t, bool>> stack(1, make_pair(start_side, false));
        return_val.first.insert(start_side);
        
        while (!stack.empty()) {
            bool going_left;
            handle_t side_here;
            tie(side_here, going_left) = stack.back();
            stack.pop_back();
            
#ifdef debug_is_bipartite
            cerr << "traversal at " << graph->get_id(side_here) << " " << graph->get_is_reverse(side_here) << ", going left? " << going_left << endl;
#endif
            
            auto& partition_across = going_left ? return_val.first : return_val.second;
            auto& partition_here = going_left ? return_val.second : return_val.first;
            
            bool still_bipartite = for_each_adjacent_side(side_here,
                                                          [&](handle_t adjacent_side) {
#ifdef debug_is_bipartite
                cerr << "\tcheck adjacent " <<  graph->get_id(adjacent_side) << " " << graph->get_is_reverse(adjacent_side) << endl;
#endif
                if (partition_here.count(adjacent_side)) {
                    // this side was seen with both even and odd parities, the
                    // adjacency component is not bipartite
#ifdef debug_is_bipartite
                    cerr << "\t\twrong parity, component is not bipartite" << endl;
#endif
                    return false;
                }
                else if (!partition_across.count(adjacent_side)) {
                    // add this side to the partition and prepare a search from it
                    // with the opposite parity
                    partition_across.insert(adjacent_side);
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

AdjacencyComponent::bipartition AdjacencyComponent::maximum_bipartite_partition_apx_1_2(uint64_t seed) const {
    
    // m = largest prime less than 2^64-1
    // a, b were generated by uniform random variables in [0, m)
    linear_congruential_engine<uint64_t,
                               12842602976330790676ull,
                               2553097010709119664ull,
                               18446744073709551557ull> gen(seed);
    uniform_int_distribution<int> coin_flip(0, 1);
    
    bipartition return_val;
    
    for (const handle_t& side : component) {
        if (coin_flip(gen)) {
            return_val.first.insert(side);
        }
        else {
            return_val.second.insert(side);
        }
    }
    
    // ensure that both sides of the partition are non-empty8
    if (return_val.first.empty() && !return_val.second.empty()) {
        return_val.first.insert(component.front());
        return_val.second.erase(component.front());
    }
    else if (return_val.second.empty() && !return_val.first.empty()) {
        return_val.second.insert(component.front());
        return_val.first.erase(component.front());
    }
    
    return return_val;
}

void AdjacencyComponent::refine_apx_partition(bipartition& partition, size_t max_opt_steps) const {
    
#ifdef debug_refinement
    cerr << "doing local search to refine an approximate partition:" << endl;
    cerr << "left" << endl;
    for (auto h : partition.first) {
        cerr << "\t" << graph->get_id(h) << " " << graph->get_is_reverse(h) << endl;
    }
    cerr << "right" << endl;
    for (auto h : partition.second) {
        cerr << "\t" << graph->get_id(h) << " " << graph->get_is_reverse(h) << endl;
    }
#endif
    
    bool no_optimizations_left = false;
    size_t opt_steps = 0;
    while (!no_optimizations_left && opt_steps < max_opt_steps) {
        // keep looping through the sides until we make an entire loop without
        // finding a greedy improvement (or until hitting max iters)
        bool new_loop = true;
        for (size_t i = 0, end = 0; (i != end || new_loop) && opt_steps < max_opt_steps;
             i = (i + 1) % component.size()) {
            
            new_loop = false;
            
            // we will check if we can improve the partition by switching
            // this side
            handle_t side = component[i];
            
            // count the edges across the partition and within it
            bool on_left = partition.first.count(side);
            if (on_left ? partition.first.size() <= 1 : partition.second.size() <= 1) {
                // we
                continue;
            }
            int edges_across = 0, edges_within = 0;
            for_each_adjacent_side(side, [&](handle_t adj_side) {
                bool nbr_on_left = partition.first.count(adj_side);
                if (on_left == nbr_on_left) {
                    ++edges_within;
                }
                else {
                    ++edges_across;
                }
                return true;
            });
            // move if the node is in Bylka, Idzik, & Tuza's S_1
            // or S_3, third component
            if (edges_within > edges_across || (edges_within == edges_across && on_left)) {
                // move it to the other side of the partition
                new_loop = true;
                end = i;
                ++opt_steps;
                if (on_left) {
                    partition.second.insert(side);
                    partition.first.erase(side);
                }
                else {
                    partition.first.insert(side);
                    partition.second.erase(side);
                }
#ifdef debug_refinement
                cerr << "found a local move by swapping node " << graph->get_id(side) << " " << graph->get_is_reverse(side) << ", partition is now:" << endl;
                cerr << "left" << endl;
                for (auto h : partition.first) {
                    cerr << "\t" << graph->get_id(h) << " " << graph->get_is_reverse(h) << endl;
                }
                cerr << "right" << endl;
                for (auto h : partition.second) {
                    cerr << "\t" << graph->get_id(h) << " " << graph->get_is_reverse(h) << endl;
                }
#endif
            }
        }
        
        // keep track of whether we find any optimization with our
        // sweep over the edges
        no_optimizations_left = true;
        for (size_t i = 0; i < component.size() && opt_steps < max_opt_steps; ++i) {
            handle_t side = component[i];
            if (partition.first.count(side)) {
                // look at swapping along edges of this side
                
                // count the edges within and across for this side
                int edges_across = 0, edges_within = 0;
                for_each_adjacent_side(side, [&](handle_t adj_side) {
                    if (partition.second.count(adj_side)) {
                        ++edges_across;
                    }
                    else {
                        ++edges_within;
                    }
                    return true;
                });
                // check all of its edges across the partition
                for_each_adjacent_side(side, [&](handle_t adj_side) {
                    
                    if (!partition.second.count(adj_side)) {
                        // this edge is not across the partition, skip over it
                        return true;
                    }
                    
                    // count up the adjacent side's edges across and within
                    int adj_edges_across = 0, adj_edges_within = 0;
                    for_each_adjacent_side(adj_side, [&](handle_t adj_adj_side) {
                        if (partition.first.count(adj_adj_side)) {
                            ++adj_edges_across;
                        }
                        else {
                            ++adj_edges_within;
                        }
                        return true;
                    });
                    
                    // move if edge is in  Bylka, Idzik, & Tuza's S_3,
                    // second component
                    if (edges_within + adj_edges_within - edges_across - adj_edges_across >= -1) {
                        // swapping the assignments of both ends of this edge will
                        // improve the partition
                        no_optimizations_left = false;
                        partition.second.insert(side);
                        partition.first.erase(side);
                        partition.first.insert(adj_side);
                        partition.second.erase(adj_side);
                        ++opt_steps;
                        
#ifdef debug_refinement
                        cerr << "found a local move by swapping edge " << graph->get_id(side) << " " << graph->get_is_reverse(side) << " -- " << graph->get_id(adj_side) << " " << graph->get_is_reverse(adj_side) << ", partition is now:" << endl;
                        cerr << "left" << endl;
                        for (auto h : partition.first) {
                            cerr << "\t" << graph->get_id(h) << " " << graph->get_is_reverse(h) << endl;
                        }
                        cerr << "right" << endl;
                        for (auto h : partition.second) {
                            cerr << "\t" << graph->get_id(h) << " " << graph->get_is_reverse(h) << endl;
                        }
#endif
                        
                        // we can't keep iterating on this node's edges because it's
                        // on the other side of the partition now
                        return false;
                    }
                    // keep iterating
                    return true;
                });
            }
        }
    }
}

uint64_t AdjacencyComponent::recursive_gray_code(bipartition& partition, uint64_t score,
                                                 size_t index, size_t to_flip, uint64_t& best_score,
                                                 uint64_t& best_index) const {
    // recursively handle the indexes to the left
    if (to_flip) {
        score = recursive_gray_code(partition, score, index - (1 << (to_flip - 1)),
                                    to_flip - 1, best_score, best_index);
    }
    
    // count how many edges are currently across or within the partition
    bool on_left = partition.first.count(component[to_flip]);
    uint64_t edges_within = 0, edges_across = 0;
    for_each_adjacent_side(component[to_flip], [&](handle_t adj_side) {
        bool adj_on_left = partition.first.count(adj_side);
        if (on_left == adj_on_left) {
            ++edges_within;
        }
        else {
            ++edges_across;
        }
        return true;
    });
    
    // switch the flipping side to the other part of the partition
    if (on_left) {
        partition.first.erase(component[to_flip]);
        partition.second.insert(component[to_flip]);
    }
    else {
        partition.second.erase(component[to_flip]);
        partition.first.insert(component[to_flip]);
    }
    // update the score for this side's edges
    // TODO: this will break if there are reversing self-loops
    score += edges_within - edges_across;
    
    if (score > best_score) {
        // this is the best score we've seen so far
        best_score = score;
        best_index = index;
    }
    
    // recursively handle the indexes to the right
    if (to_flip) {
        score = recursive_gray_code(partition, score, index + (1 << (to_flip - 1)),
                                    to_flip - 1, best_score, best_index);
    }
    
    return score;
}

AdjacencyComponent::bipartition AdjacencyComponent::exhaustive_maximum_bipartite_partition() const {
    
    // a partition that we will maintain throughout iteration
    bipartition partition;
    partition.first.reserve(component.size());
    partition.first.insert(component.begin(), component.end());
    partition.second.reserve(component.size() - 1);
    
    uint64_t best_score = 0;
    uint64_t best_index = 0;
    
    if (component.size() > 1) {
        // recursively iterate in order through a reflected binary Gray code (so that
        // only one side be swapped in the partition per iteration)
        recursive_gray_code(partition, 0, 1 << (component.size() - 2),
                            component.size() - 2, best_score, best_index);
    }
    
    // reconstruct the partition that gave rise to the best score (conversion algorithm
    // taken from https://en.wikipedia.org/wiki/Gray_code#Converting_to_and_from_Gray_code)
    uint64_t best_gray_code = best_index ^ (best_index >> 1);
    bipartition best_partition;
    for (size_t i = 0; i < component.size(); ++i) {
        if (best_gray_code & (1 << i)) {
            best_partition.second.insert(component[i]);
        }
        else {
            best_partition.first.insert(component[i]);
        }
    }
    return best_partition;
}

void AdjacencyComponent::decompose_into_bipartite_blocks(const function<void(const BipartiteGraph&)>& lambda) const {
    
    bipartition partition = bipartite_partition();
    if (partition.first.size() + partition.second.size() == component.size()) {
        // the whole component is bipartite, so there is only need for the one bipartite block
        lambda(BipartiteGraph(*graph, partition));
    }
    else {
        // TODO: magic constants
        if (component.size() < 8) {
            // the case is small enough to solve with brute force
            partition = exhaustive_maximum_bipartite_partition();
        }
        else {
            // start off with an approximate bipartition
            partition = maximum_bipartite_partition_apx_1_2();
            
            // first iteration bound heuristically derived from Kaul & West (2008)
            size_t max_opt_iters = min<size_t>(ceil(0.5 * pow(component.size(), 1.5)),
                                               10 * component.size());
            
            // use greedy local search to find a locally optimal partition (or bail early)
            refine_apx_partition(partition, max_opt_iters);
        }
        
        // divvy up the nodes and edges based on whether the edges are across
        // the partition or not
        SubtractiveHandleGraph partition_graph(*graph), remainder_graph(*graph);
        unordered_set<handle_t> partition_sides, remainder_sides;
        for (auto it = begin(), e = end(); it != e; ++it) {
            handle_t side = *it;
            bool on_left = partition.first.count(side);
            for_each_adjacent_side(side, [&](handle_t adj_side) {
                bool adj_on_left = partition.first.count(adj_side);
                if (on_left == adj_on_left) {
                    // this is a within partition edge
                    remainder_sides.insert(side);
                    remainder_sides.insert(adj_side);
                    partition_graph.subtract_edge(side, graph->flip(adj_side));
                }
                else {
                    // this is an across partition edge
                    partition_sides.insert(side);
                    partition_sides.insert(adj_side);
                    remainder_graph.subtract_edge(side, graph->flip(adj_side));
                }
                return true;
            });
        }
        
        // remove any sides from the partition that don't have any edges
        // that cross the partition
        vector<handle_t> to_remove;
        to_remove.reserve(component.size() - partition_sides.size());
        for (auto it = begin(), e = end(); it != e; ++it) {
            if (!partition_sides.count(*it)) {
                to_remove.push_back(*it);
            }
        }
        for (auto side : to_remove) {
            if (partition.first.count(side)) {
                partition.first.erase(side);
            }
            else {
                partition.second.erase(side);
            }
        }
        // execute on the part of the component that we bipartitioned
        lambda(BipartiteGraph(partition_graph, partition));
        
        // repeat the whole procedure again on the part of the component that we didn't
        // manage to bipartition
        AdjacencyComponent remainder_component(remainder_graph,
                                               remainder_sides.begin(), remainder_sides.end());
        remainder_component.decompose_into_bipartite_blocks(lambda);
    }
}

}
