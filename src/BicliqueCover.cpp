/**
 * \file BicliqueCover.cpp
 *
 * Implements algorithm for computing the biclique cover of a bipartite graph.
 */
#include "BicliqueCover.hpp"

#define debug_apx_cover
 
namespace bluntifier {

using std::unordered_set;
using std::sort;
using std::deque;
using std::make_heap;
using std::pop_heap;
using std::push_heap;
using std::priority_queue;
using std::greater;
using std::tuple;
using std::cerr;
using std::endl;

BicliqueCover::BicliqueCover(const BipartiteGraph& graph) : graph(graph) {

}

BicliqueCover::~BicliqueCover() {
    
}

vector<bipartition> BicliqueCover::get() const {
    
    vector<bipartition> return_val;

    // attempt the exact solution for domino-free graphs
    vector<pair<handle_t, vector<handle_t>>> simplifications;
    BipartiteGraph simplified = graph.simplify(simplifications);
    GaloisLattice galois_lattice(simplified);
    if (galois_lattice.is_domino_free()) {
        // the graph is domino-free, we can compute the exact cover
        return_val = galois_lattice.biclique_separator();
        unsimplify(return_val, simplifications);
    }
    else {
        // the graph is not domino-free, try the reduction algorithm
        ReducedDualGraph dual_graph(graph);
        bool is_exact;
        return_val = dual_graph.biclique_cover(is_exact);
        if (!is_exact || return_val.empty()) {
            // the reduced graph was too large to get an exact result, let's
            // also try another heuristic
            vector<bipartition> alt_cover = biclique_cover_apx();
            // locally improve both heuristic solutions
            lattice_polish(alt_cover);
            lattice_polish(return_val);
            // choose the smallest heuristic solution as the final result
            if (return_val.empty() || alt_cover.size() < return_val.size()) {
                return_val = move(alt_cover);
            }
        }
    }
    return return_val;
}

void BicliqueCover::unsimplify(vector<bipartition>& simplified_cover,
                               const vector<pair<handle_t, vector<handle_t>>>& simplifications) const {
    
    unordered_map<handle_t, unordered_set<size_t>> left_clique_membership, right_clique_membership;
    
    for (size_t i = 0; i < simplified_cover.size(); ++i) {
        for (auto side : simplified_cover[i].first) {
            left_clique_membership[side].insert(i);
        }
        for (auto side : simplified_cover[i].second) {
            right_clique_membership[side].insert(i);
        }
    }
    
    // undo the simplifications in reverse order
    for (auto rit = simplifications.rbegin(); rit != simplifications.rend(); ++rit) {
        auto it = left_clique_membership.find(rit->first);
        if (it != left_clique_membership.end()) {
            // simplification is on the left side
            for (size_t biclique_idx : it->second) {
                for (handle_t side : rit->second) {
                    simplified_cover[biclique_idx].first.insert(side);
                    left_clique_membership[side].insert(biclique_idx);
                }
            }
        }
        else {
            // simplification is on the right side
            for (size_t biclique_idx : right_clique_membership[rit->first]) {
                for (handle_t side : rit->second) {
                    simplified_cover[biclique_idx].second.insert(side);
                    right_clique_membership[side].insert(biclique_idx);
                }
            }
        }
    }
}

void BicliqueCover::lattice_polish(vector<bipartition>& cover) const {
    
    // TODO: a limit on the number of local search steps?
    
    bool completely_polished = false;
    while (!completely_polished) {
        completely_polished = true;
        for (bool on_right : {true, false}) {
            for (size_t i = 0; i < cover.size(); ) {
                auto& biclique = cover[i];
                // copy the right side to keep track of uncovered nodes
                unordered_set<handle_t> remaining = on_right ? biclique.second : biclique.first;
                vector<size_t> contained;
                for (size_t j = 0; j < cover.size() && !remaining.empty(); ++j) {
                    if (i == j) {
                        continue;
                    }
                    auto& other = cover[j];
                    // check if the right side of this biclique is contained in the other
                    bool is_contained = true;
                    auto it = on_right ? other.second.begin() : other.first.begin();
                    auto end = on_right ? other.second.end() : other.first.end();
                    for (; it != end && is_contained; ++it) {
                        is_contained = on_right ? biclique.second.count(*it) : biclique.first.count(*it);
                    }
                    if (is_contained) {
                        // mark the right side nodes in the parent as having been found
                        bool any_new = false;
                        for (auto node : on_right ? other.second : other.first) {
                            if (remaining.count(node)) {
                                remaining.erase(node);
                                any_new = true;
                            }
                        }
                        if (any_new) {
                            contained.push_back(j);
                        }
                    }
                }
                
                if (remaining.empty()) {
                    // we found a lattice cover of this biclique, expand its successors
                    auto& merging = on_right ? biclique.first : biclique.second;
                    for (size_t j : contained) {
                        auto& merge_into = on_right ? cover[j].first : cover[j].second;
                        merge_into.insert(merging.begin(), merging.end());
                    }
                    // and remove this biclique
                    cover[i] = move(cover.back());
                    cover.pop_back();
                    completely_polished = false;
                }
                else {
                    ++i;
                }
            }
        }
    }
}

vector<bipartition> BicliqueCover::biclique_cover_apx() const {
    
    vector<bipartition> return_val;
    
    // gather the edges and queue up the left side nodes based on their degree
    vector<unordered_set<handle_t>> left_uncovered_edges(graph.left_size());
    vector<unordered_set<handle_t>> right_uncovered_edges(graph.right_size());
    priority_queue<tuple<size_t, size_t, bool>, vector<tuple<size_t, size_t, bool>>, greater<tuple<size_t, size_t, bool>>> queue;
    for (auto it = graph.left_begin(); it != graph.left_end(); ++it) {
        auto& left_edges = left_uncovered_edges[it - graph.left_begin()];
        graph.for_each_adjacent_side(*it, [&](handle_t right_side) {
            right_uncovered_edges[graph.right_iterator(right_side) - graph.right_begin()].insert(*it);
            left_edges.insert(right_side);
        });
        queue.emplace(left_edges.size(), it - graph.left_begin(), true);
    }
    for (size_t i = 0; i < right_uncovered_edges.size(); ++i) {
        queue.emplace(right_uncovered_edges[i].size(), i, false);
    }
    
    // we'll keep track of how many nodes don't have all of their edges covered
    vector<bool> covered_left(graph.left_size());
    vector<bool> covered_right(graph.right_size());
    size_t num_covered_left = 0;
        
    while (num_covered_left < graph.left_size()) {
        
        // dequeue the record with the smallest degree
        auto top = queue.top();
        queue.pop();
                
        if (std::get<2>(top) ? covered_left[std::get<0>(top)] : covered_right[std::get<0>(top)]) {
            continue;
        }
        
        // make a new biclique to fill out
        return_val.emplace_back();
        auto& biclique = return_val.back();
        auto& same_side = std::get<2>(top) ? biclique.first : biclique.second;
        auto& other_side = std::get<2>(top) ? biclique.second : biclique.first;
        
        handle_t pivot = std::get<2>(top) ? *(graph.left_begin() + std::get<1>(top)) : *(graph.right_begin() + std::get<1>(top));
        
        // the right side consists of the neighbors of this node
        graph.for_each_adjacent_side(pivot, [&](handle_t other_side_node) {
            other_side.insert(other_side_node);
        });
        
        // set the left side to be the intersection of the right side's neighborhoods
        auto it = other_side.begin();
        graph.for_each_adjacent_side(*it, [&](handle_t same_side_node) {
            same_side.insert(same_side_node);
        });
        ++it;
        for (; it != other_side.end(); ++it) {
            unordered_set<handle_t> neighborhood;
            graph.for_each_adjacent_side(*it, [&](handle_t same_side_node) {
                neighborhood.insert(same_side_node);
            });
            vector<handle_t> to_erase;
            for (auto same_side_node : same_side) {
                if (!neighborhood.count(same_side_node)) {
                    to_erase.push_back(same_side_node);
                }
            }
            for (auto same_side_node : to_erase) {
                same_side.erase(same_side_node);
            }
        }
        
        // mark any newly covered edges and enqueue newly reduced nodes
        for (auto left_side : biclique.first) {
            auto& left_edges = left_uncovered_edges[graph.left_iterator(left_side) - graph.left_begin()];
            for (auto right_side : biclique.second) {
                auto& right_edges = right_uncovered_edges[graph.right_iterator(right_side) - graph.right_begin()];
                if (left_edges.count(right_side)) {
                    left_edges.erase(right_side);
                }
                if (right_edges.count(left_side)) {
                    right_edges.erase(left_side);
                }
            }
            if (left_edges.size() == 0) {
                if (!covered_left[graph.left_iterator(left_side) - graph.left_begin()]) {
                    covered_left[graph.left_iterator(left_side) - graph.left_begin()] = true;
                    ++num_covered_left;
                }
            }
            else {
                queue.emplace(left_edges.size(), graph.left_iterator(left_side) - graph.left_begin());
            }
        }
        for (auto right_side : biclique.second) {
            auto& right_edges = right_uncovered_edges[graph.right_iterator(right_side) - graph.right_begin()];
            if (right_edges.size() == 0) {
                if (!covered_right[graph.right_iterator(right_side) - graph.right_begin()]) {
                    covered_right[graph.right_iterator(right_side) - graph.right_begin()] = true;
                }
            }
            else {
                queue.emplace(right_edges.size(), graph.right_iterator(right_side) - graph.right_begin());
            }
        }
    }
    
    return return_val;
}

}
