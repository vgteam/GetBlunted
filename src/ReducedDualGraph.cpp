/**
 * \file ReducedDualGraph.cpp
 *
 * Implements an algorithm to compute biclique cover using a dual graph
 */
#include "ReducedDualGraph.hpp"

//#define debug_dual_graph

namespace bluntifier {

using std::numeric_limits;
using std::make_pair;
using std::max_element;
using std::remove_if;
using std::bitset;
using std::endl;
using std::cerr;
using std::stable_sort;

ReducedDualGraph::ReducedDualGraph(const BipartiteGraph& graph) : graph(&graph) {
    
#ifdef debug_dual_graph
    cerr << "constructing dual graph for bipartite graph with left size " << graph.left_size() << ", right size " << graph.right_size() << endl;
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
         cerr << "finding neighborhood of " << i << " (" << dual_node.first << " " << dual_node.second << ") using non-dual method" << endl;
#endif
    
        
        for (size_t j : left_edges[dual_node.first]) {
            auto& right_edge_idx = right_edge_index[j];
            for (size_t k : right_edges[dual_node.second]) {
                if (right_edge_idx.count(k)) {
                    // there is an edge between these two neighbors
                    size_t l = edge_to_dual_node[make_pair(k, j)];
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
    auto left_idx = left_adj_idx[dual_node.second];
    left_adj[left_idx] = left_adj.back();
    left_adj_idx[left_adj[left_idx]] = left_idx;
    left_adj.pop_back();
    left_adj_idx.erase(dual_node.second);
    
    // remove the corresponding edge in the right adjacency lists
    auto& right_adj = right_edges[dual_node.second];
    auto& right_adj_idx = right_edge_index[dual_node.second];
    auto right_idx = right_adj_idx[dual_node.first];
    right_adj[right_idx] = right_adj.back();
    right_adj_idx[right_adj[right_idx]] = right_idx;
    right_adj.pop_back();
    right_adj_idx.erase(dual_node.first);
    
    // move the dual node to the back
    size_t unreduced_back = reduced_size() - 1;
    swap(dual_node, dual_nodes[unreduced_back]);
    edge_to_dual_node[dual_nodes[i]] = i;
        
    // record the reduction
    reductions.emplace_back(i, dominatee);
    
#ifdef debug_dual_graph
    cerr << "reduced away " << i << " (" << dual_nodes[unreduced_back].first << " " << dual_nodes[unreduced_back].second << ") ";
    if (dominatee == numeric_limits<size_t>::max()) {
        cerr << "as isolated node";
    }
    else {
        cerr << "as dominator of " << dominatee;
    }
    cerr << ", dual node " << unreduced_back << " (" << dual_nodes[i].first << " " << dual_nodes[i].second << ") moved to position " << i << endl;
    cerr << "remaining graph:" << endl;
    print_dual_graph(cerr);
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
                    // a node can't dominate itself, and
                    if (i == j || j >= reduced_size()) {
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
                    
                    
                    
                    if (num_shared == neighborhood_size) {
                        // j dominates i, reduce by removing j
                        do_reduction(j, i);
                        fully_reduced = false;
                        return;
                    }
                    else if (num_shared == other_size) {
                        // i dominates j, reduce by removing i
                        do_reduction(i, j);
                        fully_reduced = false;
                        removed_i = true;
                        // we can't reduce with i anymore, it's gone
                        return;
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
    
#ifdef debug_dual_graph
    cerr << "making complement graph" << endl;
#endif
    
    // construct the complement graph
    vector<bool> is_left_neighbor(left_edges.size(), false);
    vector<bool> is_right_neighbor(right_edges.size(), false);
    vector<bool> is_dual_neighbor(reduced_size(), false);
    
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
    
    return VertexColoring(complement_graph).get(is_exact_out);
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
        cerr << "reversing reduction of dominator " << reduction.first << " by dominatee " << reduction.second << " swapping positions of " << endl;
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
            cerr << "assigned dominator " << reduction.first << " (" << dual_nodes[reduction.first].first << " " << dual_nodes[reduction.first].second << ") to dominatee " << reduction.second << " (" << dual_nodes[reduction.second].first << " " << dual_nodes[reduction.second].second << ")'s clique " << reduced_partition[reduction.second] << endl;
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

void ReducedDualGraph::print_dual_graph(ostream& out) {
    // prepare the banks of bools we need
    vector<bool> is_left_neighbor(left_edges.size(), false);
    vector<bool> is_right_neighbor(right_edges.size(), false);
    vector<bool> is_dual_neighbor(reduced_size(), false);
    
    for (size_t i = 0; i < reduced_size(); ++i) {
        out << i << " (" << dual_nodes[i].first << " " << dual_nodes[i].second << ")" << endl;
        dual_neighborhood_do(i, is_left_neighbor, is_right_neighbor, is_dual_neighbor,
                             [&](const vector<size_t>& dual_neighborhood) {
            for (auto j : dual_neighborhood) {
                if (i != j) {
                    out << "\t- " << j << "(" << dual_nodes[j].first << " " << dual_nodes[j].second << ")" << endl;
                }
            }
        });
    }
}

vector<bipartition> ReducedDualGraph::biclique_cover(bool& is_exact_out) {
    // get a clique partition of the reduced dual graph
    vector<size_t> partition = reduced_clique_partition(is_exact_out);
    // unwind the reduction
    reverse_reduction(partition);
    return convert_to_biclique_cover(partition);
}

}
