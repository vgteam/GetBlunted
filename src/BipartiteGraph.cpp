/**
 * \file BicliqueCover.cpp
 *
 * Implements algorithm for computing the biclique cover of a bipartite graph.
 */
#include "BipartiteGraph.hpp"

//#define debug_simplify

namespace bluntifier {

using std::cerr;
using std::endl;

// TODO: i don't know why thise needs the namespace here, but it does
BipartiteGraph::BipartiteGraph(const HandleGraph& graph,
                               const bluntifier::bipartition& partition) : graph(&graph)
{
    _partition.first.reserve(partition.first.size());
    _partition.second.reserve(partition.second.size());
    _partition.first.insert(_partition.first.end(), partition.first.begin(), partition.first.end());
    _partition.second.insert(_partition.second.end(), partition.second.begin(), partition.second.end());
    // sort to remove system dependent behavior
    sort(_partition.first.begin(), _partition.first.end());
    sort(_partition.second.begin(), _partition.second.end());
    // map the handles back to their index as well
    left_partition_index.reserve(_partition.first.size());
    right_partition_index.reserve(_partition.second.size());
    for (size_t i = 0; i < _partition.first.size(); ++i) {
        left_partition_index[_partition.first[i]] = i;
    }
    for (size_t i = 0; i < _partition.second.size(); ++i) {
        right_partition_index[_partition.second[i]] = i;
    }
    // make local adjacency lists
    left_edges.resize(_partition.first.size());
    right_edges.resize(_partition.second.size());
    for (size_t i = 0; i < _partition.first.size(); ++i) {
        graph.follow_edges(_partition.first[i], false, [&](const handle_t& right) {
            auto j = right_partition_index[graph.flip(right)];
            left_edges[i].push_back(j);
            right_edges[j].push_back(i);
        });
    }
}

BipartiteGraph::BipartiteGraph(const HandleGraph& graph,
                               const ordered_bipartition& partition) : graph(&graph), _partition(partition)
{
    // map the handles back to their index as well
    left_partition_index.reserve(_partition.first.size());
    right_partition_index.reserve(_partition.second.size());
    for (size_t i = 0; i < _partition.first.size(); ++i) {
        left_partition_index[_partition.first[i]] = i;
    }
    for (size_t i = 0; i < _partition.second.size(); ++i) {
        right_partition_index[_partition.second[i]] = i;
    }
    // make local adjacency lists
    left_edges.resize(_partition.first.size());
    right_edges.resize(_partition.second.size());
    for (size_t i = 0; i < _partition.first.size(); ++i) {
        graph.follow_edges(_partition.first[i], false, [&](const handle_t& right) {
            auto j = right_partition_index[graph.flip(right)];
            left_edges[i].push_back(j);
            right_edges[j].push_back(i);
        });
    }
}

BipartiteGraph::~BipartiteGraph() {
    
}

size_t BipartiteGraph::get_degree(handle_t node) const {
    auto it = left_partition_index.find(node);
    if (it != left_partition_index.end()) {
        return left_edges[it->second].size();
    }
    else {
        return right_edges[right_partition_index.at(node)].size();
    }
}

BipartiteGraph::const_iterator BipartiteGraph::left_begin() const {
    return _partition.first.begin();
}

BipartiteGraph::const_iterator BipartiteGraph::left_end() const {
    return _partition.first.end();
}

BipartiteGraph::const_iterator BipartiteGraph::left_iterator(const handle_t node) const {
    return _partition.first.begin() + left_partition_index.at(node);
}

size_t BipartiteGraph::left_size() const {
    return _partition.first.size();
}

BipartiteGraph::const_iterator BipartiteGraph::right_begin() const {
    return _partition.second.begin();
}

BipartiteGraph::const_iterator BipartiteGraph::right_end() const {
    return _partition.second.end();
}

BipartiteGraph::const_iterator BipartiteGraph::right_iterator(const handle_t node) const {
    return _partition.second.begin() + right_partition_index.at(node);
}

size_t BipartiteGraph::right_size() const {
    return _partition.second.size();
}

bool BipartiteGraph::is_left_side(const handle_t node) const {
    return left_partition_index.count(node);
}

void BipartiteGraph::for_each_adjacent_side(const handle_t& side,
                                            const function<void(handle_t)>& lambda) const {
    auto it = left_partition_index.find(side);
    if (it != left_partition_index.end()) {
        for (size_t j : left_edges[it->second]) {
            lambda(_partition.second[j]);
        }
    }
    else {
        for (size_t j : right_edges[right_partition_index.at(side)]) {
            lambda(_partition.first[j]);
        }
    }
}

BipartiteGraph BipartiteGraph::simplify(vector<pair<handle_t, vector<handle_t>>>& simplifications) const {
    BipartiteGraph simplifying(*graph, _partition);
    simplifying.simplify_side(simplifying._partition.first, simplifying._partition.second,
                              simplifying.left_edges, simplifying.right_edges, simplifications);
    simplifying.simplify_side(simplifying._partition.second, simplifying._partition.first,
                              simplifying.right_edges, simplifying.left_edges, simplifications);
    return simplifying;
}

void BipartiteGraph::simplify_side(const vector<handle_t>& simplifying_partition,
                                   const vector<handle_t>& opposite_partition,
                                   vector<vector<size_t>>& simplifying_edges,
                                   vector<vector<size_t>>& opposite_edges,
                                   vector<pair<handle_t, vector<handle_t>>>& simplifications) const {
#ifdef debug_simplify
    cerr << "simplifying a side" << endl;
    cerr << "left adj list:" << endl;
    for (size_t i = 0; i < _partition.first.size(); ++i) {
        cerr << i << ":";
        for (auto j : left_edges[i]) {
            cerr << " " << j;
        }
        cerr << endl;
    }
    cerr << "right adj list:" << endl;
    for (size_t i = 0; i < _partition.second.size(); ++i) {
        cerr << i << ":";
        for (auto j : right_edges[i]) {
            cerr << " " << j;
        }
        cerr << endl;
    }
#endif
    
    // matrix of successors (succ(u) in Amilhastre)
    vector<vector<bool>> successor;
    successor.reserve(simplifying_partition.size());
    for (size_t i = 0; i < simplifying_partition.size(); ++i) {
        successor.emplace_back(simplifying_partition.size(), false);
    }
    // keeps track how many successors a node has (and thereby also if it
    // has successors, i.e. LI in Amilhastre)
    vector<size_t> num_successors(simplifying_partition.size(), 0);
    
    // number of nodes in Nbd(i) \ Nbd(j)  (Delta(u,v) in Amilhastre)
    vector<vector<uint64_t>> neighbor_delta;
    
    // initialize the data structures above
    for (size_t i = 0; i < simplifying_partition.size(); ++i) {
        
        // get the neighborhood of i
        vector<bool> is_neighbor(opposite_edges.size(), false);
        for (size_t j : simplifying_edges[i]) {
            is_neighbor[j] = true;
        }
        
        // the size of this set difference starts at the degree
        neighbor_delta.emplace_back(simplifying_partition.size(), simplifying_edges[i].size());
        for (size_t j = 0; j < simplifying_partition.size(); ++j) {
            if (i == j) {
                // it's pointless to compare i to itself
                continue;
            }
            // subtract from the size of the set difference of the neighborhoods for
            // each node that these have in common
            for (size_t k : simplifying_edges[j]) {
                if (is_neighbor[k]) {
                    --neighbor_delta[i][j];
                }
            }
            
            if (neighbor_delta[i][j] == 0) {
                // the neighborhood of i is constained in the neighborhood of j, so the
                // containment preorder applies here
                successor[i][j] = true;
                ++num_successors[i];
            }
        }
    }
    
#ifdef debug_simplify
    cerr << "initial neighbor delta" << endl;
    for (auto& row : neighbor_delta) {
        for (auto val : row) {
            cerr << val << " ";
        }
        cerr << endl;
    }
    cerr << "num successors" << endl;
    for (auto val : num_successors) {
        cerr << val << " ";
    }
    cerr << endl;
    cerr << "successors" << endl;
    for (auto& row : successor) {
        for (auto val : row) {
            cerr << val << " ";
        }
        cerr << endl;
    }
#endif
    
    // TODO: there should be a more efficient way to do this in a system-independent
    // manner with a queue, but i think it doesn't affect the asymptotic run time
    
    // now we will start removing edges to simplify the graph
    bool fully_simplified = false;
    while (!fully_simplified)  {
        fully_simplified = true;
        
        // look for a simplification
        for (size_t i = 0; i < num_successors.size() && fully_simplified; ++i) {
            if (num_successors[i] == 0) {
                // we want to find a node with successors
                continue;
            }
            fully_simplified = false;
#ifdef debug_simplify
            cerr << "simplifying on " << i << endl;
#endif
            
            simplifications.emplace_back();
            simplifications.back().first = simplifying_partition[i];
            
            // find the next successor
            for (size_t j = 0; j < simplifying_partition.size(); ++j) {
                if (!successor[i][j]) {
                    continue;
                }
#ifdef debug_simplify
                cerr << "found successor " << j << endl;
#endif
                // we've found a successor of i, remove the edges in j that go to
                // neighbors of i
                simplifications.back().second.push_back(simplifying_partition[j]);
                
                for (size_t k : simplifying_edges[i]) {
                    
                    // TODO: i this could maybe blow up the run time?
                    // i might want to temporarily convert the adj lists to sets
                    
                    // remove the edge from the adjacency lists
                    auto& forward_edges = simplifying_edges[j];
                    auto& backward_edges = opposite_edges[k];
                    *find(forward_edges.begin(), forward_edges.end(), k) = forward_edges.back();
                    forward_edges.pop_back();
                    *find(backward_edges.begin(), backward_edges.end(), j) = backward_edges.back();
                    backward_edges.pop_back();
                    
#ifdef debug_simplify
                    cerr << "removed edge " << j << " -> " << k << endl;
#endif
                    
                    // update the state tracking variables
                    
                    // collect the neighbors of the other side of this edge
                    vector<bool> is_second_order_neighbor(simplifying_partition.size(), false);
                    for (size_t l : backward_edges) {
                        is_second_order_neighbor[l] = true;
                    }
                    
                    for (size_t k = 0; k < simplifying_partition.size(); ++k) {
                        if (is_second_order_neighbor[k]) {
                            // this is a neighbor of i and j's neighbor whose edge we just removed
                            
                            // there's now one more element in k's neighborhood being
                            // removed by the set difference with j's neighborhood
                            ++neighbor_delta[k][j];
#ifdef debug_simplify
                            cerr << "increment neighbor delta[" << k << "," << j << "] to " << neighbor_delta[k][j] << endl;
#endif
                            
                            if (successor[k][j]) {
                                // j can no longer be a successor of k because we took
                                // away a neighbor from j that they shared
                                successor[k][j] = false;
                                --num_successors[k];
#ifdef debug_simplify
                                cerr << j << " is no longer a successor of " << k << ", bringing num successors to " << num_successors[k] << endl;
#endif
                            }
                        }
                        else {
                            // there's now one less edge in j's neighborhood
                            --neighbor_delta[j][k];
#ifdef debug_simplify
                            cerr << "decrement neighbor delta[" << j << "," << k << "] to " << neighbor_delta[j][k] << endl;
#endif
                            
                            if (neighbor_delta[j][k] == 0 && !simplifying_edges[j].empty()) {
                                // j's neighbors are now a subset of k's neighbors
                                if (!successor[j][k]) {
                                    successor[j][k] = true;
                                    ++num_successors[j];
#ifdef debug_simplify
                                    cerr << k << " is now a successor of " << j << ", bringing num successors to " << num_successors[j] << endl;
#endif
                                }
                            }
                        }
                    }
                }
                if (successor[j][i]) {
                    // i and j were mutually successors, but no longer
                    successor[j][i] = false;
                    --num_successors[j];
#ifdef debug_simplify
                    cerr << "removing mutual successor relation on " << j << ", bringing num successors to " << num_successors[j] << endl;
#endif
                }
            }
        }
    }
}

const HandleGraph& BipartiteGraph::get_graph() const {
    return *graph;
}

}
