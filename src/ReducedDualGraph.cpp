/**
 * \file ReducedDualGraph.cpp
 *
 * Implements an algorithm to compute biclique cover using a dual graph
 */
#include "ReducedDualGraph.hpp"

//#define debug_dual_graph
//#define debug_independent_sets
//#define debug_vertex_cover

namespace bluntifier {

using std::numeric_limits;
using std::make_pair;
using std::max_element;
using std::remove_if;
using std::bitset;
using std::endl;
using std::cerr;

ReducedDualGraph::ReducedDualGraph(const BipartiteGraph& graph) : graph(&graph) {
    
#ifdef debug_dual_graph
    cerr << "constructing dual graph" << endl;
#endif
    
    // make a local, mutable copy of the bipartite graph and init
    // the dual nodes (which correspond to edges)
    left_edges.resize(graph.left_size());
    left_edge_index.resize(graph.left_size());
    right_edges.resize(graph.right_size());
    right_edge_index.resize(graph.right_size());
    for (size_t i = 0; i < left_edges.size(); ++i) {
        left_edges[i].reserve(graph.get_degree(*(graph.left_begin() + i)));
        graph.for_each_adjacent_side(*(graph.left_begin() + i),
                                     [&](handle_t adj) {
            size_t j = graph.right_iterator(adj) - graph.right_begin();
            
            left_edge_index[i][j] = left_edges[i].size();
            left_edges[i].push_back(j);
            right_edge_index[j][i] = right_edges[j].size();
            right_edges[j].push_back(i);
            
            edge_to_dual_node[make_pair(i, j)] = dual_nodes.size();
            dual_nodes.emplace_back(i, j);
#ifdef debug_dual_graph
            cerr << "\tdual node " << dual_nodes.size() - 1 << " for edge " << i << " " << j << " (" << graph.get_graph().get_id(*(graph.left_begin() + i)) << " " << graph.get_graph().get_id(adj) << ")" << endl;
#endif
        });
    }
    
    // do the reduction
    reduce();
}

void ReducedDualGraph::dual_neighborhood_do(size_t i, vector<bool>& is_left_neighbor,
                                            vector<bool>& is_right_neighbor,
                                            vector<bool>& is_dual_neighbor,
                                            const function<void(const vector<size_t>&)>& lambda) {
    auto& dual_node = dual_nodes[i];
    
    vector<size_t> dual_neighborhood;
    
    // we'll switch off between two methods of finding the neighborhood
    if (left_edges[dual_node.first].size() * right_edges[dual_node.second].size() < reduced_size()) {
        // it is faster to iterate over pairs of nodes in the non-dual neighborhood
        
#ifdef debug_dual_graph
        cerr << "finding neighborhood of " << i << " using non-dual method" << endl;
#endif
        
        for (size_t j : left_edges[dual_node.first]) {
            size_t k = edge_to_dual_node[make_pair(dual_node.first, j)];
            is_dual_neighbor[k] = true;
            dual_neighborhood.push_back(k);
        }
        for (size_t j : right_edges[dual_node.second]) {
            size_t k = edge_to_dual_node[make_pair(j, dual_node.second)];
            if (!is_dual_neighbor[k]) {
                is_dual_neighbor[k] = true;
                dual_neighborhood.push_back(k);
            }
        }
        
        for (size_t j : left_edges[dual_node.first]) {
            for (size_t k : right_edges[dual_node.second]) {
                if (right_edge_index[j].count(k) && left_edge_index[k].count(j)) {
                    // there is an edge between these two neighbors
                    size_t l = edge_to_dual_node[make_pair(dual_node.first, j)];
                    if (!is_dual_neighbor[l]) {
                        is_dual_neighbor[l] = true;
                        dual_neighborhood.push_back(l);
                    }
                }
            }
        }
    }
    else {
        // it is faster to iterate over dual nodes
        
#ifdef debug_dual_graph
        cerr << "finding neighborhood of " << i << " using dual method" << endl;
#endif
        
        // record which nodes are adjacent to the endpoints of the edge
        for (size_t j : left_edges[dual_node.first]) {
            is_right_neighbor[j] = true;
        }
        for (size_t j : right_edges[dual_node.second]) {
            is_left_neighbor[j] = true;
        }
        
        // any dual node that is adjacent on both sides in the bipartite
        // graph is adjacent in the dual graph
        for (size_t j = 0; j < reduced_size(); ++j) {
            auto& other_dual_node = dual_nodes[j];
            if (is_left_neighbor[other_dual_node.first]
                && is_right_neighbor[other_dual_node.second]) {
                is_dual_neighbor[j] = true;
                dual_neighborhood.push_back(j);
            }
        }
        
        // reset the neighbor status
        for (size_t j : left_edges[dual_node.first]) {
            is_right_neighbor[j] = false;
        }
        for (size_t j : right_edges[dual_node.second]) {
            is_left_neighbor[j] = false;
        }
    }
    
    // execute on the neighborhood
    lambda(dual_neighborhood);
    
    // reset the dual neighbor status
    for (size_t j : dual_neighborhood) {
        is_dual_neighbor[j] = false;
    }
}

void ReducedDualGraph::do_reduction(size_t i, size_t dominatee) {
    
    auto& dual_node = dual_nodes[i];
    
    // remove the corresponding edge in the left adjacency lists
    auto& left_adj = left_edges[dual_node.first];
    auto& left_adj_idx = left_edge_index[dual_node.first];
    left_adj[left_adj_idx[dual_node.second]] = left_adj.back();
    left_adj.pop_back();
    left_adj_idx.erase(dual_node.second);
    
    // remove the corresponding edge in the right adjacency lists
    auto& right_adj = right_edges[dual_node.second];
    auto& right_adj_idx = right_edge_index[dual_node.second];
    right_adj[right_adj_idx[dual_node.first]] = right_adj.back();
    right_adj.pop_back();
    right_adj_idx.erase(dual_node.first);
    
    // move the dual node to the back
    size_t unreduced_back = reduced_size() - 1;
    swap(dual_node, dual_nodes[unreduced_back]);
    edge_to_dual_node[dual_nodes[i]] = i;
        
    // record the reduction
    reductions.emplace_back(i, dominatee);
    
#ifdef debug_dual_graph
    cerr << "reduced away " << i << " (" << dual_node.first << " " << dual_node.second << ") ";
    if (dominatee == numeric_limits<size_t>::max()) {
        cerr << "as isolated node";
    }
    else {
        cerr << "as dominator of " << dominatee;
    }
    cerr << ", dual node " << unreduced_back << " (" << dual_nodes[i].first << " " << dual_nodes[i].second << ") moved to position " << i << endl;
    cerr << "remaining dual nodes:" << endl;
    for (size_t i = 0; i < reduced_size(); ++i) {
        cerr << "\t" << i << ": " << dual_nodes[i].first << " " << dual_nodes[i].second << endl;
    }
#endif
}

void ReducedDualGraph::reduce() {
    
#ifdef debug_dual_graph
    cerr << "reducing dual graph" << endl;
#endif
    
    // TODO: implement a bail out to avoid worst case O(VE^3) run time?
    
    // create the banks of bool values that we will use to make membership
    // queries
    vector<bool> is_left_neighbor(left_edges.size(), false);
    vector<bool> is_right_neighbor(right_edges.size(), false);
    vector<bool> is_dual_neighbor(dual_nodes.size(), false);
    vector<bool> is_other_left_neighbor(left_edges.size(), false);
    vector<bool> is_other_right_neighbor(right_edges.size(), false);
    vector<bool> is_other_dual_neighbor(dual_nodes.size(), false);
    
    bool fully_reduced = false;
    while (!fully_reduced) {
        
        fully_reduced = true;
        
        // check for reductions involving dominated nodes
        for (size_t i = 0; i < reduced_size();) {
            bool removed_i = false;
            dual_neighborhood_do(i, is_left_neighbor, is_right_neighbor, is_dual_neighbor,
                                 [&](const vector<size_t>& neighborhood) {
                
                size_t neighborhood_size = neighborhood.size();
                for (size_t j : neighborhood) {
                    // a node can't dominate itself
                    if (i == j) {
                        continue;
                    }
                    
                    // check if either of these two dual nodes dominate each other
                    size_t num_shared = 0;
                    size_t other_size = 0;
                    dual_neighborhood_do(j, is_other_left_neighbor, is_other_right_neighbor,
                                         is_other_dual_neighbor,
                                         [&](const vector<size_t>& other_neighborhood) {
                        for (size_t k : other_neighborhood) {
                            if (is_dual_neighbor[k]) {
                                ++num_shared;
                            }
                        }
                        other_size = other_neighborhood.size();
                    });
                    
                    
                    if (num_shared == other_size) {
                        // i dominates j, reduce by removing i
                        do_reduction(i, j);
                        fully_reduced = false;
                        removed_i = true;
                        // we can't reduce with i anymore, it's gone
                        return;
                    }
                    else if (num_shared == neighborhood_size) {
                        // j dominates i, reduce by removing j
                        do_reduction(j, i);
                        fully_reduced = false;
                        // to keep looking through i's neighbors we have to
                        // maintain these tracking variables
                        is_dual_neighbor[j] = false;
                        --neighborhood_size;
                    }
                }
            });
            
            if (!removed_i) {
                ++i;
            }
        }
        
        // check for reductions involving isolated nodes
        for (size_t i = 0; i < reduced_size();) {
            
            auto& dual_node = dual_nodes[i];
            
            if (left_edges[dual_node.first].size() == 1 &&
                right_edges[dual_node.second].size() == 1) {
                // this edge's dual node is isolated,
                do_reduction(i, numeric_limits<size_t>::max());
                fully_reduced = false;
            }
            else {
                ++i;
            }
        }
    }
#ifdef debug_dual_graph
    cerr << "finished reducing graph" << endl;
#endif
}

size_t ReducedDualGraph::reduced_size() const {
    return dual_nodes.size() - reductions.size();
}

vector<size_t> ReducedDualGraph::reduced_clique_partition(bool& is_exact_out) {
    
    vector<size_t> complement_coloring;
    if (reduced_size() == 0) {
#ifdef debug_dual_graph
        cerr << "graph is completely reduced, return trivial coloring" << endl;
#endif
        is_exact_out = true;
        return complement_coloring;
    }
    
#ifdef debug_dual_graph
    cerr << "making complement graph" << endl;
#endif
    
    // construct the complement graph
    vector<bool> is_left_neighbor(left_edges.size(), false);
    vector<bool> is_right_neighbor(right_edges.size(), false);
    vector<bool> is_dual_neighbor(reduced_size(), false);
    
    size_t total_degree = 0;
    vector<vector<size_t>> complement_graph(reduced_size());
    for (size_t i = 0; i < reduced_size(); ++i) {
        dual_neighborhood_do(i, is_left_neighbor, is_right_neighbor, is_dual_neighbor,
                             [&](const vector<size_t>& dual_neighborhood) {
            auto& compl_adj = complement_graph[i];
            compl_adj.reserve(reduced_size() - dual_neighborhood.size());
            for (size_t j = 0; j < reduced_size(); ++j) {
                if (!is_dual_neighbor[j]) {
                    complement_graph[i].push_back(j);
                }
            }
        });
        total_degree += complement_graph[i].size();
    }
    
#ifdef debug_dual_graph
    cerr << "complement graph:" << endl;
    for (size_t i = 0; i < complement_graph.size(); ++i) {
        cerr << "\t" << i << ":";
        for (auto j : complement_graph[i]) {
            cerr << " " << j;
        }
        cerr << endl;
    }
#endif
    
    // estimate the time usage of the exact algorithm (constant is 1 + 3^(1/3), see Lawler 1976)
    size_t num_edges = total_degree / 2;
    size_t num_nodes = complement_graph.size();
    size_t cost = num_nodes * num_edges * size_t(ceil(pow(2.44225, num_nodes)));
    
    // TODO: magic number
    // compute a vertex coloring of the complement graph (use a fairly generous bound
    // because in practice the optimizations seem to speed it up quite a bit)
    if (cost <= (1 << 20)) {
#ifdef debug_dual_graph
        cerr << "computing exact vertex coloring" << endl;
#endif
        complement_coloring = vertex_coloring_exact(complement_graph);
        is_exact_out = true;
    }
    else {
#ifdef debug_dual_graph
        cerr << "computing approximate vertex coloring" << endl;
#endif
        complement_coloring = vertex_coloring_apx(complement_graph);
        is_exact_out = false;
    }
    
    return complement_coloring;
}

vector<uint16_t> ReducedDualGraph::maximal_independent_sets(const vector<vector<size_t>>& graph,
                                                            const vector<size_t>& subgraph_trans) const {
    
    // IS in Tsukiyama
    vector<uint16_t> intersection_size(graph.size(), 0);
    // Bucket in Tsukiyama
    vector<vector<uint16_t>> reset_buckets(graph.size());
    vector<uint16_t> maximal_ind_sets;
    maximal_independent_sets_internal(graph, subgraph_trans, 0, intersection_size,
                                      reset_buckets, maximal_ind_sets);
    return maximal_ind_sets;
}

void ReducedDualGraph::maximal_independent_sets_internal(const vector<vector<size_t>>& graph,
                                                         const vector<size_t>& subgraph_trans,
                                                         size_t max_node,
                                                         vector<uint16_t>& intersection_size,
                                                         vector<vector<uint16_t>>& reset_buckets,
                                                         vector<uint16_t>& maximal_ind_sets) const {
    
    if (max_node + 1 == graph.size()) {
        // we've reached the end, the nodes whose intersection doesn't have any
        // nodes in common with the maximal indpendent set are its members
        uint16_t ind_set = 0;
        for (size_t i = 0; i < intersection_size.size(); ++i) {
            if (intersection_size[i] == 0) {
                ind_set |= (1 << subgraph_trans[i]);
            }
        }
        maximal_ind_sets.push_back(ind_set);
#ifdef debug_independent_sets
        cerr << "at a full MIS, with current IS" << endl;
        for (size_t i = 0; i < intersection_size.size(); ++i) {
            cerr << "\t" << i << ": " << intersection_size[i] << endl;
        }
        cerr << "yielding result " << bitset<16>(ind_set) << endl;
#endif
        return;
    }
    
    // x in Tsukiyama
    size_t next_node = max_node + 1;
    
#ifdef debug_independent_sets
    cerr << "MIS iteration with next node to be added " << next_node << endl;
    cerr << "current IS" << endl;
    for (size_t i = 0; i < intersection_size.size(); ++i) {
        cerr << "\t" << i << ": " << intersection_size[i] << endl;
    }
#endif
    
    // c in Tsukiyama
    size_t nbr_count = 0;
    auto& next_adj = graph[next_node];
    for (size_t i = 0; i < next_adj.size() && next_adj[i] <= max_node; ++i) {
        if (intersection_size[next_adj[i]] == 0) {
            nbr_count += 1;
        }
    }
    if (nbr_count == 0) {
#ifdef debug_independent_sets
        cerr << "this IS doesn't overlap with neighborhood of " << next_node << endl;
#endif
        for (size_t i = 0; i < next_adj.size() && next_adj[i] <= max_node; ++i) {
            intersection_size[next_adj[i]] += 1;
        }
        // procedure B1
        maximal_independent_sets_internal(graph, subgraph_trans, next_node, intersection_size,
                                          reset_buckets, maximal_ind_sets);
        for (size_t i = 0; i < next_adj.size() && next_adj[i] <= max_node; ++i) {
            intersection_size[next_adj[i]] -= 1;
        }
    }
    else {
#ifdef debug_independent_sets
        cerr << "this IS satisfies condition b" << endl;
#endif
        intersection_size[next_node] = nbr_count;
        maximal_independent_sets_internal(graph, subgraph_trans, next_node, intersection_size,
                                          reset_buckets, maximal_ind_sets);
        intersection_size[next_node] = 0;
        
#ifdef debug_independent_sets
        cerr << "checking for condition b" << endl;
#endif
        auto& bucket = reset_buckets[next_node];
        
        // f in Tsukiyama, called condition B in exposition
        bool condition_b = true;
        for (size_t i = 0; i < next_adj.size() && next_adj[i] <= max_node; ++i) {
            // y in Tsukiyama
            auto nbr = next_adj[i];
            if (intersection_size[nbr] == 0) {
                bucket.push_back(nbr);
                auto& nbr_adj = graph[nbr];
                for (size_t j = 0; j < nbr_adj.size() && nbr_adj[j] <= max_node; ++j) {
                    // z in Tsukiyama
                    auto nbr_nbr = nbr_adj[j];
                    intersection_size[nbr_nbr] -= 1;
                    condition_b &= (intersection_size[nbr_nbr] != 0);
                }
            }
            intersection_size[nbr] += 1;
        }
#ifdef debug_independent_sets
        cerr << "done checking condition b, reset bucket for " << next_node << " contains" << endl;
        for (auto j : bucket) {
            cerr << "\t" << j << endl;
        }
#endif
        if (condition_b) {
#ifdef debug_independent_sets
            cerr << "this IS satisfies condition b" << endl;
#endif
            maximal_independent_sets_internal(graph, subgraph_trans, next_node, intersection_size,
                                              reset_buckets, maximal_ind_sets);
        }
        for (size_t i = 0; i < next_adj.size() && next_adj[i] <= max_node; ++i) {
            intersection_size[next_adj[i]] -= 1;
        }
        while (!bucket.empty()) {
            // y in Tsukiyama
            auto nbr = bucket.back();
            bucket.pop_back();
            auto& nbr_adj = graph[nbr];
            for (size_t i = 0; i < nbr_adj.size() && nbr_adj[i] <= max_node; ++i) {
                // z in Tsukiyama
                intersection_size[nbr_adj[i]] += 1;
            }
        }
#ifdef debug_independent_sets
        cerr << "IS after resetting with bucket" << endl;
        for (size_t i = 0; i < intersection_size.size(); ++i) {
            cerr << "\t" << i << ": " << intersection_size[i] << endl;
        }
#endif
    }
}

vector<size_t> ReducedDualGraph::vertex_coloring_exact(const vector<vector<size_t>>& complement_graph) const {
    
    if (complement_graph.size() > 15) {
        cerr << "only can compute exact vertex coloring for <= 15 vertices" << endl;
        exit(1);
    }
    
    // we need to ensure that adjacency lists are ordered by index for the MIS algorithm
    vector<vector<size_t>> ordered_graph(complement_graph.size());
    for (size_t i = 0; i < complement_graph.size(); ++i) {
        for (size_t j : complement_graph[i]) {
            ordered_graph[j].push_back(i);
        }
    }
    
    vector<size_t> subset_graph_trans;
    vector<vector<size_t>> subset_graph;
    subset_graph_trans.reserve(ordered_graph.size());
    subset_graph.reserve(ordered_graph.size());
    
    // fill subset_graph and subset_graph_trans with the induced subgraph corresponding
    // to the bitset
    auto load_subset_graph = [&](uint16_t subset) {
        // prepare to load the new subgraph
        subset_graph_trans.clear();
        subset_graph.clear();
        
        // figure out how to compress the indexes
        vector<size_t> removed_before(ordered_graph.size(), 0);
        for (size_t i = 0; i + 1 < ordered_graph.size(); ++i) {
            if (subset & (1 << i)) {
                removed_before[i + 1] = removed_before[i];
            }
            else {
                removed_before[i + 1] = removed_before[i] + 1;
            }
        }
        for (size_t i = 0; i < removed_before.size(); ++i) {
        }
        
        // copy into the subgraph
        for (size_t i = 0; i < ordered_graph.size(); ++i) {
            if ((subset >> i) & 1) {
                subset_graph.emplace_back();
                subset_graph_trans.emplace_back(i);
                auto& edges = subset_graph.back();
                for (size_t adj : ordered_graph[i]) {
                    if ((subset >> adj) & 1) {
                        edges.emplace_back(adj - removed_before[adj]);
                    }
                }
            }
        }
#ifdef debug_vertex_cover
        cerr << "loaded subset graph:" << endl;
        for (size_t i = 0; i < subset_graph.size(); ++i) {
            cerr << i << " (" << subset_graph_trans[i] << "):";
            for (size_t j : subset_graph[i]) {
                cerr << " " << j << "(" << subset_graph_trans[j] << ")";
            }
            cerr << endl;
        }
#endif
    };
    
    // returns a bipartite partition as bitsets if one exists, otherwise
    // a pair of sentinels
    auto bipartite_partition = [&]() {
        pair<uint16_t, uint16_t> return_val(numeric_limits<uint16_t>::max(),
                                            numeric_limits<uint16_t>::max());
        
        static uint8_t odd_parity = 0, even_parity = 1, unvisited = 2;
        
        vector<uint8_t> parity(subset_graph.size(), unvisited);
        
        for (size_t i = 0; i < subset_graph.size(); ++i) {
            // TODO: this also finds connected components, which can be used
            // to find further independent sets
            if (parity[i] != unvisited) {
                continue;
            }
            
            // DFS, looking for odd cycles
            vector<size_t> stack(1, i);
            parity[i] = even_parity;
            while (!stack.empty()) {
                size_t here = stack.back();
                stack.pop_back();
                for (size_t j : subset_graph[here]) {
                    if (parity[j] == unvisited) {
                        parity[j] = parity[here] ^ 1;
                        stack.push_back(j);
                    }
                    else if (parity[j] != (parity[here] ^ 1)) {
                        // found an odd cycle, graph is not bipartite
                        return return_val;
                    }
                    else {
                    }
                }
            }
        }
        
        // the DFS didn't find any odd cycles, convert the parity into
        // two independent sets in the full graph
        return_val.first = 0;
        return_val.second = 0;
        for (size_t i = 0; i < subset_graph.size(); ++i) {
            if (parity[i] == even_parity) {
                return_val.first |= (1 << subset_graph_trans[i]);
            }
            else {
                return_val.second |= (1 << subset_graph_trans[i]);
            }
        }
        return return_val;
    };
    
    // each DP record consists of (chromatic number, backpointer)
    // indexed by the bitset corresponding to the nodes in a subgraph
    size_t num_subsets = 1 << ordered_graph.size();
    vector<pair<uint16_t, uint16_t>> dp(num_subsets, make_pair(numeric_limits<uint16_t>::max(),
                                                               numeric_limits<uint16_t>::max()));
    
    // base case
    dp[0].first = 0;
    
    // stack records of (subset, parent set), we use the stack to explicitly do
    // recursion rather than filling out the whole DP table in hopes that the space
    // is somewhat sparse in the middle
    vector<pair<uint16_t, uint16_t>> stack;
    stack.emplace_back(num_subsets - 1, numeric_limits<uint16_t>::max());
    
    while (!stack.empty()) {
        
        uint16_t set_here = stack.back().first;
        uint16_t parent = stack.back().second;
        
#ifdef debug_vertex_cover
        cerr << "destack " << bitset<16>(set_here) << " with parent " << bitset<16>(parent) << endl;
#endif
        
        if (dp[set_here].first != numeric_limits<uint16_t>::max()) {
            // we've computed the DP value for this node subset already
#ifdef debug_vertex_cover
            cerr << "\tDP has been computed" << endl;
#endif
            if (parent != numeric_limits<uint16_t>::max() && dp[parent].first > dp[set_here].first + 1) {
#ifdef debug_vertex_cover
                cerr << "\t\tupdate parent to " << dp[set_here].first << endl;
#endif
                // adding this independent set improves on the parent's chromatic number
                dp[parent] = pair<uint16_t, uint16_t>(dp[set_here].first + 1, set_here);
            }
            stack.pop_back();
            continue;
        }
        
        // make the corresponding graph
        load_subset_graph(set_here);
        
        // Lawler's third optimization
        auto bipartite_ind_sets = bipartite_partition();
        if (bipartite_ind_sets.first != numeric_limits<uint16_t>::max()) {
            // the graph is bipartite, yielding a trivial 2-coloring
#ifdef debug_vertex_cover
            cerr << "\tsubgraph is bipartite with partial sets " << bitset<16>(bipartite_ind_sets.first) << ", " << bitset<16>(bipartite_ind_sets.second) << endl;
#endif
            
            dp[set_here] = pair<uint16_t, uint16_t>(2, bipartite_ind_sets.first);
            // moreover the components are independent sets
            dp[bipartite_ind_sets.first] = pair<uint16_t, uint16_t>(1, 0);
            dp[bipartite_ind_sets.second] = pair<uint16_t, uint16_t>(1, 0);
            
            if (parent != numeric_limits<uint16_t>::max() && dp[parent].first > 3) {
                // adding this independent set improves on the parent's chromatic number
                dp[parent] = pair<uint16_t, uint16_t>(3, set_here);
            }
            stack.pop_back();
            continue;
        }
        
        size_t max_degree_node = 0;
        for (size_t i = 0; i < subset_graph.size(); ++i) {
            if (subset_graph[i].size() > subset_graph[max_degree_node].size()) {
                max_degree_node = i;
            }
            // TODO: Lawler's fourth opimzation?
            // i think do actually do this correctly you have to backtrace
            // down the DP to figure out which independent set it should be added to, which
            // gets hard because using the first optimization means that the corresponding
            // DP problem after removing that independent set might not be computed, which
            // means you won't be able to backtrace in the final step (except if we maybe
            // restrict to degree 2 nodes, because then the bipartite routine guarantees it)
        }
        
        uint16_t mask_filter = 1 << subset_graph_trans[max_degree_node];
        
#ifdef debug_vertex_cover
        cerr << "\tchoose node " << subset_graph_trans[max_degree_node] << " with subset degree " << subset_graph[max_degree_node].size() << " for filter: " << bitset<16>(mask_filter) << endl;
#endif
        
        // queue up the next DP problems
        for (uint16_t ind_set : maximal_independent_sets(subset_graph, subset_graph_trans)) {
#ifdef debug_vertex_cover
            cerr << "\trecurse with independent set " << bitset<16>(ind_set) << endl;
#endif
            
            // Lawler's first optimization, only consider sets containing
            // one arbitrary node
            if (mask_filter & ind_set) {
                stack.emplace_back(set_here ^ ind_set, set_here);
            }
#ifdef debug_vertex_cover
            else {
                cerr << "\t\tset doesn't pass filter" << endl;
            }
#endif
            
            // Lawler's second optimization
            dp[ind_set] = pair<uint16_t, uint16_t>(1, 0);
        }
    }
    
#ifdef debug_vertex_cover
    cerr << "obtain final chromatic number " << dp.back().first << endl;
#endif
    
    // backtrace the independent sets of the optimal solution
    vector<uint16_t> ind_sets;
    ind_sets.reserve(dp.back().first);
    uint16_t set_here = num_subsets - 1;
    while (set_here != 0) {
        uint16_t next_set = dp[set_here].second;
        ind_sets.emplace_back(set_here ^ next_set);
        set_here = next_set;
    }
    
    // unpack the independent sets into a vertex coloring
    vector<size_t> return_val(ordered_graph.size(), 0);
    for (size_t i = 1; i < ind_sets.size(); ++i) {
        for (size_t j = 0; j < ordered_graph.size(); ++j) {
            if (ind_sets[i] & (1 << j)) {
                return_val[j] = i;
            }
        }
    }
    return return_val;
}

vector<size_t> ReducedDualGraph::vertex_coloring_apx(const vector<vector<size_t>>& complement_graph) const {
    
    // aggregate by degree
    vector<vector<size_t>> degree_sets;
    for (size_t i = 0; i < complement_graph.size(); ++i) {
        size_t degree = complement_graph[i].size();
        while (degree_sets.size() <= degree) {
            degree_sets.emplace_back();
        }
        degree_sets[degree].push_back(i);
    }
    
    // convert the degree sets into an ordering and an index lookup
    // descending by degree
    vector<size_t> ordering(complement_graph.size());
    vector<size_t> index(complement_graph.size());
    for (size_t ordering_idx = 0, i = degree_sets.size(); i > 0; --i) {
        auto& degree_set = degree_sets[i - 1];
        for (size_t j = 0; j < degree_set.size(); ++j, ++ordering_idx) {
            ordering[ordering_idx] = degree_set[j];
            index[degree_set[j]] = ordering_idx;
        }
    }
    
#ifdef debug_dual_graph
    cerr << "degree ordering:" << endl;
    for (auto i : ordering) {
        cerr << "\t" << i << endl;
    }
#endif
    
    vector<size_t> coloring(complement_graph.size());
    vector<bool> color_used;
    for (size_t i : ordering) {
        
#ifdef debug_dual_graph
        cerr << "assigning color to " << i << endl;
#endif
        
        // mark the colors of the predecessor neighbors used
        for (auto j : complement_graph[i]) {
            if (index[j] < index[i]) {
#ifdef debug_dual_graph
                cerr << "\tpredecessor " << j << " has color " << coloring[j] << endl;
#endif
                color_used[coloring[j]] = true;
            }
        }
        
        // find the lowest valued unused color
        size_t color = numeric_limits<size_t>::max();
        for (size_t j = 0; j < color_used.size(); ++j) {
            if (!color_used[j]) {
                color = j;
#ifdef debug_dual_graph
                cerr << "assigning previous color " << color << endl;
#endif
                break;
            }
        }
        
        if (color == numeric_limits<size_t>::max()) {
            // all of the colors are already used, we need a new color
            color = color_used.size();
            color_used.emplace_back(false);
#ifdef debug_dual_graph
            cerr << "assigning new color " << color << endl;
#endif
        }
        
        coloring[i] = color;
        
        // reset the used status for the neighbor's colors
        for (auto j : complement_graph[i]) {
            if (index[j] < index[i]) {
                color_used[j] = false;
            }
        }
    }
    
    return coloring;
}

void ReducedDualGraph::reverse_reduction(vector<size_t>& reduced_partition) {
    
#ifdef debug_dual_graph
    cerr << "reversing reductions" << endl;
#endif
    
    size_t next_clique_id = 0;
    if (!reduced_partition.empty()) {
        next_clique_id = 1 + *max_element(reduced_partition.begin(), reduced_partition.end());
    }
    
    reduced_partition.resize(dual_nodes.size());
    // iterate over reductions in reverse order
    while (!reductions.empty()) {
        
        auto reduction = reductions.back();
        reductions.pop_back();
        size_t unreduced_back = reduced_size() - 1;
        
#ifdef debug_dual_graph
        cerr << "reversing reduction " << reduction.first << " " << reduction.second << " swapping positions of " << endl;
        cerr << "\t" << reduction.first << ": " << dual_nodes[reduction.first].first << " " << dual_nodes[reduction.first].second << endl;
        cerr << "\t" << unreduced_back << ": " << dual_nodes[unreduced_back].first << " " << dual_nodes[unreduced_back].second << endl;
#endif
        
        // swap the reduced node out from the back position
        swap(dual_nodes[reduction.first], dual_nodes[unreduced_back]);
        reduced_partition[unreduced_back] = reduced_partition[reduction.first];
        if (reduction.second == numeric_limits<size_t>::max()) {
            // this was an isolated dual node, it gets its own new clique
            reduced_partition[reduction.first] = next_clique_id;
            ++next_clique_id;
#ifdef debug_dual_graph
            cerr << "assigned " << reduction.first << " to new clique " << next_clique_id - 1 << endl;
#endif
        }
        else {
            // this was a dominator, it gets assigned to the clique of the dominatee
            reduced_partition[reduction.first] = reduced_partition[reduction.second];
#ifdef debug_dual_graph
            cerr << "assigned " << reduction.first << " to dominatee " << reduction.second << "'s clique " << reduced_partition[reduction.second] << endl;
#endif
        }
    }
}

vector<bipartition> ReducedDualGraph::convert_to_biclique_cover(vector<size_t>& clique_partition) const {
    size_t num_bicliques = 1 + *max_element(clique_partition.begin(), clique_partition.end());
#ifdef debug_dual_graph
    cerr << "computed dual clique cover of size " << num_bicliques << ", converting to bicliques" << endl;
#endif
    vector<bipartition> return_val(num_bicliques);
    for (size_t i = 0; i < clique_partition.size(); ++i) {
        auto& biclique = return_val[clique_partition[i]];
        auto& dual_node = dual_nodes[i];
        biclique.first.insert(*(graph->left_begin() + dual_node.first));
        biclique.second.insert(*(graph->right_begin() + dual_node.second));
    }
    return return_val;
}

vector<bipartition> ReducedDualGraph::biclique_cover(bool& is_exact_out) {
    // get a clique partition of the reduced dual graph
    vector<size_t> partition = reduced_clique_partition(is_exact_out);
    // unwind the reduction
    reverse_reduction(partition);
    return convert_to_biclique_cover(partition);
}

}
