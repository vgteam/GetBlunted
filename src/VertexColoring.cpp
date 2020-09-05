/**
 * \file VertexColoring.cpp
 *
 * Implements algorithm for computing vertex coloring
 */
#include "VertexColoring.hpp"
 
namespace bluntifier {

using std::numeric_limits;
using std::cerr;
using std::endl;
using std::pair;
using std::make_pair;

VertexColoring::VertexColoring(const vector<vector<size_t>>& graph) {
    // create a local copy of the graph that is ordered by index in the
    // adjacency lists
    this->graph.resize(graph.size());
    for (size_t i = 0; i < graph.size(); ++i) {
        this->graph[i].reserve(graph[i].size());
    }
    for (size_t i = 0; i < graph.size(); ++i) {
        for (size_t j : graph[i]) {
            this->graph[j].emplace_back(i);
        }
    }
}

vector<size_t> VertexColoring::get(bool& is_exact) const {
    
    vector<size_t> coloring;
    
    size_t total_degree = 0;
    for (size_t i = 0; i < graph.size(); ++i) {
        total_degree += graph[i].size();
    }
    // estimate the time usage of the exact algorithm (constant is 1 + 3^(1/3), see Lawler 1976)
    size_t num_edges = total_degree / 2;
    size_t num_nodes = graph.size();
    size_t max_cost = num_nodes * num_edges * size_t(ceil(pow(2.44225, num_nodes)));
    
    // TODO: magic number
    // compute a vertex coloring of the complement graph (use a fairly generous bound
    // because in practice the optimizations seem to speed it up quite a bit)
    if (max_cost <= (1 << 20) && num_nodes < 16) {
#ifdef debug_vertex_coloring
        cerr << "computing exact vertex coloring" << endl;
#endif
        coloring = lawlers_algorithm();
        is_exact = true;
    }
    else {
#ifdef debug_vertex_coloring
        cerr << "computing approximate vertex coloring" << endl;
#endif
        auto order = degree_ordering();
        coloring = greedy_coloring(order);
        is_exact = false;
    }
    
    return coloring;
}

vector<size_t> VertexColoring::degree_ordering() const {
    // aggregate by degree
    vector<vector<size_t>> degree_sets;
    for (size_t i = 0; i < graph.size(); ++i) {
        size_t degree = graph[i].size();
        while (degree_sets.size() <= degree) {
            degree_sets.emplace_back();
        }
        degree_sets[degree].push_back(i);
    }
    
    // convert the degree sets into an ordering and an index lookup
    // descending by degree
    vector<size_t> ordering(graph.size());
    for (size_t ordering_idx = 0, i = degree_sets.size(); i > 0; --i) {
        auto& degree_set = degree_sets[i - 1];
        for (size_t j = 0; j < degree_set.size(); ++j, ++ordering_idx) {
            ordering[ordering_idx] = degree_set[j];
        }
    }
    
    return ordering;
}

vector<size_t> VertexColoring::lawlers_algorithm() const {
    if (graph.size() > 15) {
        cerr << "Lawler's algorithm implementation only valid for <= 15 vertices" << endl;
        exit(1);
    }
    
    vector<size_t> subset_graph_trans;
    vector<vector<size_t>> subset_graph;
    subset_graph_trans.reserve(graph.size());
    subset_graph.reserve(graph.size());
    
    // fill subset_graph and subset_graph_trans with the induced subgraph corresponding
    // to the bitset
    auto load_subset_graph = [&](uint16_t subset) {
        // prepare to load the new subgraph
        subset_graph_trans.clear();
        subset_graph.clear();
        
        // figure out how to compress the indexes
        vector<size_t> removed_before(graph.size(), 0);
        for (size_t i = 0; i + 1 < graph.size(); ++i) {
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
        for (size_t i = 0; i < graph.size(); ++i) {
            if ((subset >> i) & 1) {
                subset_graph.emplace_back();
                subset_graph_trans.emplace_back(i);
                auto& edges = subset_graph.back();
                for (size_t adj : graph[i]) {
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
    size_t num_subsets = 1 << graph.size();
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
    vector<size_t> return_val(graph.size(), 0);
    for (size_t i = 1; i < ind_sets.size(); ++i) {
        for (size_t j = 0; j < graph.size(); ++j) {
            if (ind_sets[i] & (1 << j)) {
                return_val[j] = i;
            }
        }
    }
    return return_val;
}

vector<size_t> VertexColoring::greedy_coloring(const vector<size_t>& order) const {
    
    // convert the degree sets into an ordering and an index lookup
    // descending by degree
    vector<size_t> index(graph.size());
    for (size_t i = 0; i < graph.size(); ++ i) {
        index[order[i]] = i;
    }
    
#ifdef debug_vertex_coloring
    cerr << "degree ordering:" << endl;
    for (auto i : order) {
        cerr << "\t" << i << endl;
    }
#endif
    
    vector<size_t> coloring(graph.size());
    vector<bool> color_used;
    for (size_t i : order) {
        
#ifdef debug_vertex_coloring
        cerr << "assigning color to " << i << endl;
#endif
        
        // mark the colors of the predecessor neighbors used
        for (auto j : graph[i]) {
            if (index[j] < index[i]) {
#ifdef debug_vertex_coloring
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
#ifdef debug_vertex_coloring
                cerr << "assigning previous color " << color << endl;
#endif
                break;
            }
        }
        
        if (color == numeric_limits<size_t>::max()) {
            // all of the colors are already used, we need a new color
            color = color_used.size();
            color_used.emplace_back(false);
#ifdef debug_vertex_coloring
            cerr << "assigning new color " << color << endl;
#endif
        }
        
        coloring[i] = color;
        
        // reset the used status for the neighbor's colors
        for (auto j : graph[i]) {
            if (index[j] < index[i]) {
                color_used[coloring[j]] = false;
            }
        }
    }
    
    return coloring;
}

vector<uint16_t> VertexColoring::maximal_independent_sets(const vector<vector<size_t>>& subgraph,
                                                          const vector<size_t>& subgraph_trans) const {
    // IS in Tsukiyama
    vector<uint16_t> intersection_size(subgraph.size(), 0);
    // Bucket in Tsukiyama
    vector<vector<uint16_t>> reset_buckets(subgraph.size());
    vector<uint16_t> maximal_ind_sets;
    maximal_independent_sets_internal(subgraph, subgraph_trans, 0, intersection_size,
                                      reset_buckets, maximal_ind_sets);
    return maximal_ind_sets;
    
}

void VertexColoring::maximal_independent_sets_internal(const vector<vector<size_t>>& subgraph,
                                                       const vector<size_t>& subgraph_trans,
                                                       size_t max_node,
                                                       vector<uint16_t>& intersection_size,
                                                       vector<vector<uint16_t>>& reset_buckets,
                                                       vector<uint16_t>& maximal_ind_sets) const {
    
    if (max_node + 1 == subgraph.size()) {
        // we've reached the end, the nodes whose neighborhood doesn't have any
        // nodes in common with the maximal independent set are its members
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
    auto& next_adj = subgraph[next_node];
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
        maximal_independent_sets_internal(subgraph, subgraph_trans, next_node, intersection_size,
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
        maximal_independent_sets_internal(subgraph, subgraph_trans, next_node, intersection_size,
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
                auto& nbr_adj = subgraph[nbr];
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
            maximal_independent_sets_internal(subgraph, subgraph_trans, next_node, intersection_size,
                                              reset_buckets, maximal_ind_sets);
        }
        for (size_t i = 0; i < next_adj.size() && next_adj[i] <= max_node; ++i) {
            intersection_size[next_adj[i]] -= 1;
        }
        while (!bucket.empty()) {
            // y in Tsukiyama
            auto nbr = bucket.back();
            bucket.pop_back();
            auto& nbr_adj = subgraph[nbr];
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

}
