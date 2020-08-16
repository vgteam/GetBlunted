/**
 * \file biclique_cover.cpp
 *
 * Implements algorithm for computing the biclique cover of a bipartite graph.
 */
#include "biclique_cover.hpp"

namespace bluntifier {

using std::unordered_set;
using std::sort;

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
        // compute the biclique cover on the simplified graph using the algorithm
        // of Amilhastre, et al. 1998
        return_val = domino_free_cover();
    }
    else {
        // use the heuristic algorithm of Ene, et al. 2008
        return_val = heuristic_cover();
    }
    return return_val;
}

bool CenteredGaloisTree::CenteredGaloisTree(const BicliqueCover& cover,
                                            handle_t center) const {
    
    // get the two-hop subgraph starting at the center
    unordered_map<handle_t, size_t> left_idx;
    // we have to restrict rightward edges since some of them could
    // point outside the subgraph
    vector<vector<size_t>> left_edges;
    vector<handle_t> left_nodes, right_nodes;
    
    cover.for_each_adjacent_side(center, [&](handle_t right) {
        cover.for_each_adjacent_side(right, [&](handle_t left) {
            auto f = left_idx.find(left);
            if (f == left_idx.end()) {
                left_idx[left] = left_edges.size();
                left_edges.emplace_back(1, right_nodes.size());
                left_nodes.emplace_back(left);
            }
            else {
                left_edges[f->second].push_back(right_nodes.size());
            }
            return true;
        });
        right_nodes.emplace_back(right);
        return true;
    });
    
    // initialize every node on the left in the same equivalence class
    vector<size_t> equiv_class_assignment(left_nodes.size(),
                                          numeric_limits<size_t>::max());
    size_t next_equiv_class = 0;
    for (handle_t right : right_nodes) {
        // refine the classes using the edges of this node
        // TODO: this coud be done without an unordered_map by reseting
        // a vector after every iteration
        unordered_map<size_t, size_t> equiv_mapping;
        cover.for_each_adjacent_side(right, [&](handle_t left) {
            size_t eq_class = equiv_class_assignment[left_idx[left]];
            auto it = equiv_mapping.find(eq_class);
            if (it != equiv_mapping.end()) {
                // we've already refined this class, just look it up
                equiv_class_assignment[left_idx[left]] = it->second;
            }
            else {
                // refine with a new class
                equiv_mapping[eq_class] = next_equiv_class;
                equiv_class_assignment[left_idx[left]] = next_equiv_class;
                ++next_equiv_class;
            }
            return true;
        });
    }
    
    // quotient the nodes by the equivalence classes and compact the
    // equivalence class identifiers
    vector<vector<size_t>> equiv_classes_left_edges;
    equiv_classes_left_edges.reserve(next_equiv_class - right_nodes.size())
    
    vector<size_t> compacted_equiv_class(next_equiv_class, numeric_limits<size_t>::max());
    for (size_t i = 0; i < left_nodes.size(); ++i) {
        size_t eq_class = equiv_class_assignment[i];
        if (compacted_equiv_class[eq_class] == numeric_limits<size_t>::max()) {
            // we're coming across this class for the first time, assign it the next
            // compacted class identifier
            compacted_equiv_class[eq_class] = equiv_classes.size();
            eq_class = equiv_classes.size();
            equiv_classes.emplace_back();
            
            // bogart the edges from the original graph
            equiv_classes_left_edges.emplace_back(move(left_edges[i]));
            
            // collect and remember its neighborhood on the right side
            neighborhoods.emplace_back();
            auto& neighborhood = neighborhoods.back();
            neighborhood.reserve(equiv_classes_left_edges.back().size());
            for (auto j : equiv_classes_left_edges.back()) {
                neighborhood.push_back(right_nodes[j]);
            }
        }
        // add this left side node to the partition
        equiv_classes[eq_class].push_back(left_nodes[i]));
    }
    
    // partition left nodes by their degree (T_x(k) in Amilhastre)
    vector<vector<size_t>> degree_groups(right_nodes.size());
    for (size_t i = 0; i < neighborhoods.size(); ++i) {
        degree_groups[neighborhoods[i].size()].push_back(i);
    }
    
    // organize the neighborhoods of the right nodes in degree ordering
    // (V(y) in Amilhastre)
    vector<vector<size_t>> degree_ordered_nbds(right_nodes.size());
    for (const auto& degree_group : degree_groups) {
        for (auto left : degree_group) {
            for (auto right : equiv_classes_left_edges[left]) {
                degree_ordered_nbds[right].push_back(left);
            }
        }
    }
    
    // immediate successor in the equiv class ordering (Succ in Amilhastre)
    successors.resize(equiv_classes.size(), numeric_limits<size_t>::max());
    // all immediate predecessors in the equiv class ordering (Gamma_x^-
    // in Amilhastre)
    equiv_class_predecessors.resize(equiv_classes.size());
    
    // construct the predecessor lists and check for tree structure
    for (size_t i = 0; i < right_nodes.size(); ++i) {
        auto& degree_ordered_nbd = degree_ordered_nbds[i];
        size_t pred = degree_ordered_nbd.front();
        for (size_t j = 1; j < degree_ordered_nbd.size(); ++j) {
            size_t succ = degree_ordered_nbds[j];
            if (successors[pred] != numeric_limits<size_t>::max()) {
                successors[pred] = succ;
                equiv_class_predecessors[succ].push_back(pred);
            }
            else if (successors[pred] != succ) {
                // the successors don't form a tree, clear to mark as a failure
                clear();
                return;
            }
            pred = succ;
        }
    }
    
    // check for the proper containent relationships between the neighborhoods
    // of the nodes on the left
    for (size_t i = 0; i < equiv_classes_left_edges.size(); ++i) {
        auto& succ_nbd = equiv_classes_left_edges[i];
        for (auto j : equiv_class_predecessors[i]) {
            auto& pred_nbd = equiv_classes_left_edges[j];
            
            // note: all of the left edge lists are constructed in sorted order
            size_t p = 0;
            for (size_t s = 0; s < succ_nbd.size() && p < pred_nbd.size(); ++s) {
                if (succ_nbd[s] == pred_nbd[p]) {
                    ++p;
                }
            }
            if (p < pred_nbd.size()) {
                // the neighborhoods weren't contained, clear to mark as a failure
                clear();
                return;
            }
        }
    }
}
    
void CenteredGaloisTree::clear() {
    equiv_classes.clear();
    neighborhoods.clear();
    successors.clear();
    predecessors.clear();
}

bool CenteredGaloisTree::has_neighbor_ordering_property() {
    // did we clear everything and return in the constructor?
    return !equiv_classes.empty();
}

size_t CenteredGaloisTree::size() const {
    return equiv_classes.size();
}

const vector<size_t>& CenteredGaloisTree::predecessors(size_t i) const {
    return equiv_class_predecessors[i];
}

size_t CenteredGaloisTree::central_equivalence_class() const {
    size_t i = 0;
    while (successors[i] != numeric_limits<size_t>::max()) {
        i = successors[i];
    }
    return i;
}

size_t CenteredGaloisTree::successor(size_t i) const {
    return successors[i];
}

CenteredGaloisTree::edge_iterator CenteredGaloisTree::edge_begin(size_t i) const {
    return edge_iterator(0, 0, i, this);
}

CenteredGaloisTree::edge_iterator CenteredGaloisTree::edge_end(size_t i) const {
    return edge_iterator(equiv_classes[i].size(), 0, i, this);
}

CenteredGaloisTree::edge_iterator::edge_iterator(size_t left, size_t right, size_t eq_class,
                                                 const CenteredGaloisTree* iteratee)
    : left(left), right(right), eq_class(eq_class), iteratee(iteratee)
{
    
}


CenteredGaloisTree::edge_iterator& CenteredGaloisTree::edge_iterator::operator++() {
    if (right == iteratee->neighborhoods[eq_class].size()) {
        right = 0;
        ++left;
    }
    else {
        ++right;
    }
}

pair<handle_t, handle_t> operator*() const {
    return make_pair(iteratee->equiv_classes[eq_class][left],
                     iteratee->neighborhoods[equiv_classes][right])
}

bool operator==(const edge_iterator& other) const {
    return (left == other.left && right == other.right &&
            eq_class == other.eq_class && iteratee == other.iteratee)
}

bool operator!=(const edge_iterator& other) const {
    return !(*this == other);
}


pair<handle_t, handle_t> CenteredGaloisTree::arbitrary_edge_of(size_t i) const {
    return make_pair(equiv_classes[i].front(), neighborhood[i].front());
}

bipartition CenteredGaloisTree::biclique(size_t i) const {
    bipartition return_val;
    return_val.second.insert(neighborhood[i].begin(), neighborhood[i].end());
    while (i != numeric_limits<size_t>::max()) {
        return_val.first.insert(equiv_classes[i].begin(), equiv_classes[i].end());
        i = successors[i];
    }
    return return_val;
}

bool BicliqueCover::is_domino_free() const {
    bool domino_free = true;
    // TODO: embarassingly parallel
    for (size_t i = 0; i < left_partition.size() && domino_free; ++i) {
        CenteredGaloisTree tree(*this, left_partition[i]);
        domino_free = tree.has_neighbor_ordering_property();
    }
    return domino_free;
}

SubtractiveHandleGraph BicliqueCover::simplify() const {
    SubtractiveHandleGraph simplified(graph);
    simplify_side(left_partition, simplified);
    simplify_side(right_partition, simplified);
    return simplified;
}

GaloisLattice BicliqueCover::get_galois_lattice(const HandleGraph& simple_graph) const {
    
}

vector<bipartition> BicliqueCover::domino_free_cover() const {
    // simplify the graph without affecting the biclique cover (Amilhastre, et al.
    // 1998 algorithm 2).
    SubtractiveHandleGraph simplified = simplify();
    GaloisLattice galois_lattice = get_galois_lattice(simplified);
    // TODO: get separator of galois lattice
    // TODO: convert separator to biclique cover
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
    
    // matrix of successors (succ(u) in Amilhastre)
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
                            if (!successors[k][j]) {
                                successors[k][j] = true;
                                ++num_successors[k];
                            }
                        }
                    }
                });
            }
            nonmaximal[i] = false;
        }
    }
}

}
