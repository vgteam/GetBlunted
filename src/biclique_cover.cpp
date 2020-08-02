/**
 * \file biclique_cover.cpp
 *
 * Implements algorithm for computing the biclique cover of a bipartite graph.
 */
#include "biclique_cover.hpp"

namespace bluntifier {

using std::unordered_set;

BicliqueCover::BicliqueCover(const HandleGraph& graph,
                             const bipartition& partition)
: graph(graph),
  left_partition(partition.first.begin() partition.first.end()),
  right_partition(partition.second.begin(), partition.second.end())
{
    // sort to remove system dependent behavior
    sort(left_partition.begin(), left_partition.end());
    sort(right_partition.begin(), right_partition.end());
    // map the handles back to their index as well
    left_partition_index.reserve(left_partition.size());
    right_partition_index.reserve(right_partition.size());
    for (size_t i = 0; i < left_partition.size(); ++i) {
        left_partition_index[left_partition[i]] = i;
    }
    for (size_t i = 0; i < left_partition.size(); ++i) {
        right_partition_index[right_partition[i]] = i;
    }
}

BicliqueCover::~BicliqueCover() {
    
}

vector<bipartition> BicliqueCover::get() const {
    vector<bipartition> return_val;
    size_t edge_count = 0;
    for (handle_t handle : left_partition) {
        edge_count += graph.get_degree(handle, false);
    }
    // TODO: magic number
    if (edge_count * (left_partition.size() + right_partition.size()) <= 65536
        && is_domino_free()) {
        
        // simplify the graph without affecting the biclique cover (Amilhastre, et al.
        // 1998 algorithm 2).
        SubtractiveHandleGraph simplified = simplify();
        
        //
        return_val = domino_free_cover(simplified);
    }
    else {
        return_val = heuristic_cover();
    }
    return return_val;
}


bool BicliqueCover::is_domino_free() const {
    
}

SubtractiveHandleGraph BicliqueCover::simplify() const {
    SubtractiveHandleGraph simplified(graph);
    simplify_side(left_partition, simplified);
    simplify_side(right_partition, simplified);
    return simplified;
}

vector<bipartition> BicliqueCover::domino_free_cover() const {
    
}
vector<bipartition> BicliqueCover::heuristic_cover() const {
    
}

bool BicliqueCover::for_each_adjacent_side(const handle_t& side,
                                           const function<bool(handle_t)>& lambda) const {
    return graph->follow_edges(side, false, [&](const handle_t& neighbor) {
        return lambda(graph->flip(neighbor));
    });
}

void BicliqueCover::simplify_side(const vector<handle_t>& simplifying_partition,
                                  SubtractiveHandleGraph& simplifying) const {
    
    // keeps track of which nodes have successors (LI in Amilhastre, et al. 1998)
    vector<bool> nonmaximal(simplifying_partition.size(), false);
    
    // matrix of successors (succ(u) in Amilhastre, et al. 1998)
    vector<vector<bool>> successor;
    successor.reserve(simplifying_partition.size());
    for (size_t i = 0; i < simplifying_partition.size(); ++i) {
        successor.emplace_back(simplifying_partition.size(), false);
    }
    vector<size_t> num_successors(simplifying_partition.size(), 0);
    
    // the degree of each node
    vector<size_t> degree(simplifying_partition.size());
    // number of nodes in Nbd(i) \ Nbd(j)  (Delta(u,v) in Amilhastre)
    vector<vector<uint64_t>> neighbor_delta;
    
    // initialize the data structures above
    for (size_t i = 0; i < simplifying_partition.size(); ++i) {
        
        // get the neighborhood of i
        unordered_set<handle_t> neighborhood;
        graph.follow_edges(simplifying_partition[i], false, [&](const handle_t& nbr) {
            neighborhood.insert(nbr);
        });
        degree[i] = neighborhood.size();
        
        // the size of this set difference starts at the degree
        neighbor_delta.emplace_back(simplifying_partition.size(), neighborhood.size());
        for (size_t j = 0; j < simplifying_partition.size(); ++j) {
            if (i == j) {
                // it's pointless to compare i to itself
                continue;
            }
            // subtract from the size of the set difference of the neighborhoods for
            // each node that these have in common
            graph.follow_edges(simplifying_partition[j], false, [&](const handle_t& nbr) {
                if (neighborhood.count(nbr)) {
                    --neighbor_delta[i][j];
                }
            });
            
            if (neighbor_delta[i][j] == 0) {
                // the neighborhood of i is constained in the neighborhood of j, so the
                // containment preorder applies here
                successor[i][j] = true;
                nonmaximal[i] = true;
                ++num_successors[i];
            }
        }
    }
    
    // TODO: there should be a more efficient way to do this in a system-independent
    // manner with a queue, but i think it doesn't affect the asymptotic run time
    
    // now we will start removing edges to simplify the graph
    bool fully_simplified = false;
    while (!fully_simplified)  {
        fully_simplified = true;
        
        // look for a simplification
        for (size_t i = 0; i < nonmaximal.size() && fully_simplified; ++i) {
            if (!nonmaximal[i]) {
                // we want to find a node with successors
                continue;
            }
            fully_simplified = false;
            // find the next successor
            for (size_t j = 0; j < simplifying_partition.size(); ++j) {
                if (!successors[i][j]) {
                    continue;
                }
                // we've found a successor of i, remove the edges in j that go to
                // neighbors of i
                simplifying.follow_edges(simplifying_partition[i], false, [&](const handle_t& nbr) {
                    
                    simplifying.subtract_edge(simplifying_partition[j], nbr);
                    --degree[j];
                    
                    // update the state tracking variables
                    
                    // collect the neighbors of the other side of this edge
                    unordered_set<handle_t> nbr_nbrs;
                    simplifying.follow_edges(nbr, true, [&](const handle_t& nbr_nbr) {
                        nbr_nbrs.insert(nbr_nbr);
                    }):
                    
                    for (size_t k = 0; k < simplifying_partition.size(); ++k) {
                        // there's now one less edge in j's neighborhood
                        --neighbor_delta[j][k];
                        
                        if (nbr_nbrs.count(simplifying_partition[k])) {
                            // this is a neighbor of i and j's neighbor whose edge we just removed
                            
                            // there's now one more element in k's neighborhood being
                            // removed by the set difference with j's neighborhood
                            ++neighbor_delta[k][j];
                            if (nonmaximal[k]) {
                                if (successors[k][j]) {
                                    // j can no longer be a successor of k because we took
                                    // away a neighbor from j that they shared
                                    successors[k][j] = false;
                                    --num_successors[k];
                                }
                                if (num_successors[k] == 0) {
                                    // j was the last successor of k, so k is now maximal
                                    nonmaximal[k] = false;
                                }
                            }
                        }
                        
                        if (neighbor_delta[j][k] == 0 && degree[j] > 0) {
                            // j's neighbors are now a subset of k's neighbors
                            nonmaximal[j] = true;
                            successors[k][j] = true;
                            ++num_successors[k];
                        }
                    }
                });
            }
        }
    }
}

}
