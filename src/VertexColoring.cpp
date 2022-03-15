/**
 * \file VertexColoring.cpp
 *
 * Implements algorithm for computing vertex coloring
 */
#include "VertexColoring.hpp"
#include <algorithm>
#include <bitset>
#include <cassert>

//#define debug_interchange
//#define debug_independent_sets
//#define debug_vertex_cover
//#define debug_vertex_coloring

namespace bluntifier {

using std::numeric_limits;
using std::cerr;
using std::endl;
using std::pair;
using std::make_pair;
using std::move;
using std::swap;
using std::deque;
using std::max_element;
using std::bitset;

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
    size_t m = total_degree / 2;
    size_t n = graph.size();
    size_t max_lawler_cost = m * n * size_t(round(pow(2.44225, n)));
    size_t max_interchange_cost = m * n;

    // TODO: magic number
    // compute a vertex coloring of the complement graph (use a fairly generous bound
    // because in practice the optimizations seem to speed it up quite a bit)
    if (max_lawler_cost <= (1 << 20) && graph.size() < 16) {
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

        size_t max_color = numeric_limits<size_t>::max();
        
        // use greedy and (maybe) interchange coloring algorithms on
        auto attempt_greedy_coloring = [&](const vector<size_t>& order) {
            
#ifdef debug_vertex_coloring
            cerr << "use greedy coloring algorithm" << endl;
#endif
            auto candidate_coloring = greedy_coloring(order);
            auto candidate_max_color = *max_element(candidate_coloring.begin(),
                                                    candidate_coloring.end());
            if (candidate_max_color < max_color) {
                max_color = candidate_max_color;
                coloring = move(candidate_coloring);
#ifdef debug_vertex_coloring
                cerr << "found new coloring minimum " << max_color << endl;
                for (size_t i = 0; i < coloring.size(); ++i) {
                    cerr << "\t" << i << ": " << coloring[i] << endl;
                }
#endif
            }
            
            // TODO: magic number
            if (max_interchange_cost < (1 << 15)) {
                
#ifdef debug_vertex_coloring
                cerr << "use interchange greedy coloring" << endl;
#endif
                // the size is small enough that we can use the interchange algorithm
                // without too much cost
                auto ix_candidate_coloring = interchange_greedy_coloring(order);
                auto ix_candidate_max_color = *max_element(ix_candidate_coloring.begin(),
                                                           ix_candidate_coloring.end());
                if (ix_candidate_max_color < max_color) {
                    max_color = ix_candidate_max_color;
                    coloring = move(ix_candidate_coloring);
#ifdef debug_vertex_coloring
                    cerr << "found new coloring minimum " << max_color << endl;
                    for (size_t i = 0; i < coloring.size(); ++i) {
                        cerr << "\t" << i << ": " << coloring[i] << endl;
                    }
#endif
                }
            }
        };
        
        // get a lower bound so we can check for optimality
        size_t lower_bnd = lower_bound();
        
#ifdef debug_vertex_coloring
        cerr << "using least first ordering" << endl;
#endif
        
        // greedy coloring by least first ordering
        attempt_greedy_coloring(least_first_ordering());
        
        // greedy coloring by degree first ordering
        if (max_color != lower_bnd) {
#ifdef debug_vertex_coloring
            cerr << "using degree ordering" << endl;
#endif
            attempt_greedy_coloring(degree_ordering());
        }
        
        // greedy coloring by two different connected sequence orderings
        if (max_color != lower_bnd) {
            
#ifdef debug_vertex_coloring
            cerr << "using depth first ordering" << endl;
#endif
            attempt_greedy_coloring(depth_first_order());
        }
        if (max_color != lower_bnd) {
            
#ifdef debug_vertex_coloring
            cerr << "using bredth first ordering" << endl;
#endif
            attempt_greedy_coloring(breadth_first_order());
        }
        
        // greedy coloring by log(n) random orderings
        size_t logn = ceil(log(graph.size()));
        uint64_t seed = 14847024944434445584ull;
        for (size_t i = 0; i < logn && max_color != lower_bnd; ++i) {
            
#ifdef debug_vertex_coloring
            cerr << "using random ordering " << i << endl;
#endif
            attempt_greedy_coloring(random_ordering(seed));
        }
        
        is_exact = (max_color == lower_bnd);
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

vector<size_t> VertexColoring::least_first_ordering() const {
    
    // put the nodes into levels according to their degree
    vector<vector<size_t>> degree_queue;
    vector<size_t> degree_remaining(graph.size());
    for (size_t i = 0; i < graph.size(); ++i) {
        auto degree = graph[i].size();
        degree_remaining[i] = degree;
        while (degree_queue.size() <= degree) {
            degree_queue.emplace_back();
        }
        degree_queue[degree].push_back(i);
    }
    
    vector<size_t> ordering(graph.size());
    size_t ordering_idx = graph.size();
    vector<bool> ordered(graph.size(), false);
    
    size_t i = 0;
    while (ordering_idx > 0) {
        if (degree_queue[i].empty()) {
            // this degree set is empty, go to higher degrees
            ++i;
        }
        else {
            // choose an arbitrary elemment of this degree set
            auto& degree_level = degree_queue[i];
            auto j = degree_level.back();
            degree_level.pop_back();
            if (!ordered[j]) {
                // we haven't added this to the ordering yet, put it in
                // the last available spot
                ordering[--ordering_idx] = j;
                ordered[j] = true;
                // all of this nodes neighbors degree value reduces by 1
                // (we don't worry about removing the other record because
                // it will be filtered out)
                for (auto k : graph[j]) {
                    if (!ordered[k]) {
                        degree_queue[--degree_remaining[k]].emplace_back(k);
                    }
                }
                // the smallest non empty degree set may now be one spot lower
                --i;
            }
        }
    }
    return ordering;
}

vector<size_t> VertexColoring::random_ordering(uint64_t& seed) const {
    // m = largest prime less than 2^64-1
    // a, b were generated by uniform random variables in [0, m)
    static const uint64_t a = 14494735964225224815ull;
    static const uint64_t b = 11973402385035416413ull;
    static const uint64_t m = 18446744073709551557ull;
    
    vector<size_t> ordering(graph.size());
    for (size_t i = 0; i < graph.size(); ++i) {
        ordering[i] = i;
        seed = (a * seed + b) % m;
        std::swap(ordering[i], ordering[seed % (i + 1)]);
    }
    return ordering;
}

vector<size_t> VertexColoring::depth_first_order() const {
    vector<bool> enqueued(graph.size(), false);
    vector<size_t> order(graph.size());
    size_t ordering_idx = 0;
    for (size_t i = 0; i < graph.size(); ++i) {
        if (enqueued[i]) {
            continue;
        }
        vector<size_t> stack(1, i);
        enqueued[i] = true;
        while (!stack.empty()) {
            auto here = stack.back();
            stack.pop_back();
            order[ordering_idx] = here;
            ++ordering_idx;
            for (auto j : graph[here]) {
                if (!enqueued[j]) {
                    enqueued[j] = true;
                    stack.push_back(j);
                }
            }
        }
    }
    return order;
}

vector<size_t> VertexColoring::breadth_first_order() const {
    vector<bool> enqueued(graph.size(), false);
    vector<size_t> order(graph.size());
    size_t ordering_idx = 0;
    for (size_t i = 0; i < graph.size(); ++i) {
        if (enqueued[i]) {
            continue;
        }
        deque<size_t> queue(1, i);
        enqueued[i] = true;
        while (!queue.empty()) {
            auto here = queue.front();
            queue.pop_front();
            order[ordering_idx] = here;
            ++ordering_idx;
            for (auto j : graph[here]) {
                if (!enqueued[j]) {
                    enqueued[j] = true;
                    queue.push_back(j);
                }
            }
        }
    }
    return order;
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
    
    // TODO: reverse the ordering (DP starting with full set as 0) so that it can be
    // used for branch and bound as we go
    
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
    cerr << "ordering:" << endl;
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

vector<size_t> VertexColoring::interchange_greedy_coloring(const vector<size_t>& order) const {
        
#ifdef debug_interchange
    cerr << "graph:" << endl;
    for (size_t i = 0; i < graph.size(); ++i) {
        cerr << i << ":";
        for (auto j : graph[i]) {
            cerr << " " << j;
        }
        cerr << endl;
    }
    cerr << "ordering:" << endl;
    for (auto i : order) {
        cerr << "\t" << i << endl;
    }
#endif
    
    // TODO: reimplement using bank of records in a vector for better cache efficiency?
    
    // set up the linked list representation of the graph that we'll be manipulating
    vector<InterchangeEdge*> interx_adj_list(graph.size(), nullptr);
    for (size_t i = 0; i < graph.size(); ++i) {
        const auto& adj = graph[i];
        for (size_t j = 0; j < adj.size() && adj[j] < i; ++j) {
            
            auto dest = adj[j];
            auto edge_to = new InterchangeEdge();
            auto edge_from = new InterchangeEdge();
            
#ifdef debug_interchange
            cerr << "edge pair between " << i << " " << dest << endl;
            cerr << "\t" << edge_to << endl;
            cerr << "\t" << edge_from << endl;
#endif
            
            edge_from->vertex = dest;
            edge_from->next = interx_adj_list[i];
            interx_adj_list[i] = edge_from;
            edge_from->mate = edge_to;
            
            edge_to->vertex = i;
            edge_to->next = interx_adj_list[dest];
            interx_adj_list[dest] = edge_to;
            edge_to->mate = edge_from;
        }
    }
    
    vector<vector<InterchangeEdge*>> color_lists(graph.size());
    
    // convert the degree sets into an ordering and an index lookup
    vector<size_t> index(graph.size());
    for (size_t i = 0; i < graph.size(); ++i) {
        index[order[i]] = i;
    }
    
    vector<size_t> coloring(graph.size());
    for (size_t i : order) {
#ifdef debug_interchange
        cerr << "coloring node " << i << endl;
#endif
        
        auto& color_list = color_lists[i];
        
        // identify the edges that this node has in this iteration by color
        for (auto edge = interx_adj_list[i]; edge != nullptr; edge = edge->next) {
            if (index[edge->vertex] > index[i]) {
                // we're only looking at nodes that are earlier in the ordering
                continue;
            }
            auto nbr_color = coloring[edge->vertex];
            
            edge->color_next = color_list[nbr_color];
            color_list[nbr_color] = edge;
            // we don't add this edge to the other color lists yet, because
            // we haven't given i a color
        }
        
#ifdef debug_interchange
        cerr << "current color topology:" << endl;
        for (size_t j = 0; j < graph.size(); ++j) {
            if (index[j] <= index[i]) {
                auto& color_list = color_lists[j];
                cerr << j;
                if (j != i) {
                    cerr << " (" << coloring[j] << ")";
                    
                }
                for (size_t k = 0; k < color_list.size(); ++k) {
                    cerr << "\t" << k << ":";
                    for (auto e = color_list[k]; e != nullptr; e = e->color_next) {
                        cerr << " " << e->vertex;
                    }
                    cerr << endl;
                }
                if (color_list.size() == 0) {
                    cerr << endl;
                }
            }
        }
#endif
        
        // find the lowest valued unused color
        size_t lowest_unused_color = numeric_limits<size_t>::max();
        for (size_t j = 0; j < color_list.size(); ++j) {
            if (!color_list[j]) {
                lowest_unused_color = j;
                break;
            }
        }
        
        if (lowest_unused_color != numeric_limits<size_t>::max()) {
#ifdef debug_interchange
            cerr << "colorable with " << lowest_unused_color << endl;
#endif
            // we can use an already-existing color
            coloring[i] = lowest_unused_color;
        }
        else {
#ifdef debug_interchange
            cerr << "no available color, checking for interchange" << endl;
#endif
            
            // TODO: initalize this only one time for whole algorithm to remove O(n) here?
            // init some structures that we'll reuse for traveresals
            vector<bool> queued(graph.size(), false);
            vector<size_t> stack;
            vector<size_t> this_traversal;
            
            // attempt to find a pair of colors that can do an interchange
            size_t ix_color_1, ix_color_2;
            bool found_swap = false;
            for (ix_color_1 = 1; ix_color_1 < color_list.size() && !found_swap; ++ix_color_1) {
                // provisionally color this node
                coloring[i] = ix_color_1;
                for (ix_color_2 = 0; ix_color_2 < ix_color_1 && !found_swap; ++ix_color_2) {
                    // do an alternating-color DFS traversal with this color pair
                    
                    // init with the current node's color neighbors
                    for (auto edge = color_list[ix_color_2]; edge != nullptr; edge = edge->color_next) {
                        stack.push_back(edge->vertex);
                        queued[edge->vertex] = true;
                        this_traversal.emplace_back(edge->vertex);
                    }
                    while (!stack.empty()) {
                        size_t here = stack.back();
                        stack.pop_back();
                        size_t nbr_color = coloring[here] == ix_color_1 ? ix_color_2 : ix_color_1;
                        for (auto edge = color_lists[here][nbr_color]; edge != nullptr; edge = edge->color_next) {
                            // i should not have a color yet, so it shouldn't appear in color lists
                            assert(edge->vertex != i);
                            if (!queued[edge->vertex]) {
                                queued[edge->vertex] = true;
                                this_traversal.push_back(edge->vertex);
                                stack.emplace_back(edge->vertex);
                            }
                        }
                    }
                    
                    // did we encounter any of the current node's neighbors from the other color?
                    found_swap = true;
                    for (auto edge = color_list[ix_color_1]; edge != nullptr && found_swap; edge = edge->color_next) {
                        found_swap = !queued[edge->vertex];
                    }
                    
                    // reset the traversal marker
                    for (auto j : this_traversal) {
                        queued[j] = false;
                    }
                }
            }
            
            if (found_swap) {
#ifdef debug_interchange
                cerr << "found swappable pair" << endl;
#endif
                // we found a pair of colors that can do an interchange
                
                // undo the increment that happened at the end of the for loop
                --ix_color_1;
                --ix_color_2;
                
                // do a DFS traversal of the 2-color component to swap the colors
                InterchangeEdge* last_edge = nullptr;
                for (auto edge = color_list[ix_color_2]; edge != nullptr; edge = edge->color_next) {
                    stack.push_back(edge->vertex);
                    queued[edge->vertex] = true;
                    last_edge = edge;
                }
                while (!stack.empty()) {
                    size_t here = stack.back();
                    stack.pop_back();
                    size_t new_color = coloring[here] == ix_color_1 ? ix_color_2 : ix_color_1;
                    coloring[here] = new_color;
                    for (auto edge = color_lists[here][new_color]; edge != nullptr; edge = edge->color_next) {
                        if (!queued[edge->vertex]) {
                            queued[edge->vertex] = true;
                            stack.emplace_back(edge->vertex);
                        }
                    }
                    // the color adjacency lists can have their heads swapped
                    auto& color_list_here = color_lists[here];
                    swap(color_list_here[ix_color_1], color_list_here[ix_color_2]);
                }
                // recolor the edges from this node
                last_edge->color_next = color_list[ix_color_1];
                color_list[ix_color_1] = color_list[ix_color_2];
                color_list[ix_color_2] = nullptr;
                coloring[i] = ix_color_2;
            }
            else {
#ifdef debug_interchange
                cerr << "no swaps possible, adding a color: " << color_list.size() << endl;
#endif
                // no interchange is possible, choose the next color
                coloring[i] = color_list.size();
                // expand the color lists to accommodate the new color
                for (auto& clist : color_lists) {
                    clist.emplace_back(nullptr);
                }
            }
        }
        
        // update this node's coloring in the color lists of its neighbors
        auto color = coloring[i];
        for (auto edge = interx_adj_list[i]; edge != nullptr; edge = edge->next) {
            if (index[edge->vertex] > index[i]) {
                // we're only looking at nodes that are earlier in the ordering
                continue;
            }
            auto mate = edge->mate;
            auto& nbr_color_list = color_lists[edge->vertex];
            mate->color_next = nbr_color_list[color];
            nbr_color_list[color] = mate;
        }
    }
    
    // free the memory from the linked lists
    for (auto list_node : interx_adj_list) {
        while (list_node != nullptr) {
            auto next_node = list_node->next;
            delete list_node;
            list_node = next_node;
        }
    }
    
    return coloring;
    
}

size_t VertexColoring::lower_bound() const {
    // TODO: any clique finding algorithm or Hoffman's bound
    // TODO: Alon, Yuster, Zwick 1997 triangle detection O(m^(3/2))?
    for (const auto& adj : graph) {
        if (!adj.empty()) {
            return 1;
        }
    }

    return 0;
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
