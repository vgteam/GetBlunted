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
using std::greater;
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
    // TODO
}

vector<bipartition> BicliqueCover::biclique_cover_apx() const {
    
    vector<bipartition> return_val;
    
    static greater<pair<size_t, size_t>> cmp;
    
    // records of (remaining edge, node idx)
    vector<pair<size_t, size_t>> queue;
    queue.reserve(graph.left_size());
    
    // collect the uncovered edges (i.e. all of them, at this point) and compute
    // the set of nodes with a given degree
    vector<unordered_set<handle_t>> uncovered_edges_left(graph.left_size());
    vector<unordered_set<handle_t>> uncovered_edges_right(graph.right_size());
    for (auto it = graph.left_begin(); it != graph.left_end(); ++it) {
        auto& edges = uncovered_edges_left[it - graph.left_begin()];
        graph.for_each_adjacent_side(*it, [&](handle_t side) {
            edges.insert(side);
            uncovered_edges_right[graph.right_iterator(side) - graph.right_begin()].insert(*it);
        });
        queue.emplace_back(edges.size(), it - graph.left_begin());
    }
    
    vector<bool> covered(graph.left_size(), false);
    size_t num_covered = 0;
    
    make_heap(queue.begin(), queue.end(), cmp);
    while ((num_covered < graph.left_size()) && !queue.empty()) {
        // TODO: use a heap with optimal asymptotics (e.g. fibonacci)?
        cerr << num_covered << " of " << graph.left_size() << " " << (num_covered < graph.left_size()) << " " << queue.size() << " " << queue.empty() << " " << ((num_covered < graph.left_size()) && !queue.empty()) << endl;
        // get the node with the fewest uncovered edges
        auto top = queue.front();
        pop_heap(queue.begin(), queue.end(), cmp);
        queue.pop_back();
        
#ifdef debug_apx_cover
        cerr << "next node is " << graph.get_graph().get_id(*(graph.left_begin() + top.second)) << " with " << top.first << " uncovered edges" << endl;
#endif
        if (covered[top.second]) {
#ifdef debug_apx_cover
            cerr << "this node has been covered since being queued" << endl;
#endif
            continue;
        }
        
        // the uncovered neighbors of this pivot element will be the right side of the biclique
        return_val.emplace_back();
        auto& biclique = return_val.back();
        graph.for_each_adjacent_side(*(graph.left_begin() + top.second), [&](handle_t side) {
            biclique.second.insert(side);
        });
        
        // start with the neighbors of an arbitrary right node as the left side of the clique
        auto it = biclique.second.begin();
        graph.for_each_adjacent_side(*it, [&](handle_t side) {
            biclique.first.insert(side);
        });
        ++it;
        // reduce the left side to the intersection of the neighborhoods
        for (; it != biclique.second.end(); ++it) {
            unordered_set<handle_t> neighborhood;
            graph.for_each_adjacent_side(*it, [&](handle_t side) {
                neighborhood.insert(side);
            });
            vector<handle_t> to_remove;
            for (auto side : biclique.first) {
                if (!neighborhood.count(side)) {
                    to_remove.push_back(side);
                }
            }
            for (auto side : to_remove) {
                biclique.first.erase(side);
            }
        }
        
        // mark edges as covered and add new heap records for nodes that
        // still have uncovered edges
        for (auto left_side : biclique.first) {
            auto& left_edges = uncovered_edges_left[graph.left_iterator(left_side) - graph.left_begin()];
            for (auto right_side : biclique.second) {
                auto& right_edges = uncovered_edges_right[graph.right_iterator(right_side) - graph.right_begin()];
                auto lit = left_edges.find(right_side);
                if (lit != left_edges.end()) {
                    left_edges.erase(lit);
                }
                auto rit = right_edges.find(left_side);
                if (rit != right_edges.end()) {
                    right_edges.erase(rit);
                }
            }
            if (left_edges.size() == 0) {
                if (!covered[graph.left_iterator(left_side) - graph.left_begin()]) {
#ifdef debug_apx_cover
                    cerr << "marking " << graph.get_graph().get_id(left_side) << " covered" << endl;
#endif
                    covered[graph.left_iterator(left_side) - graph.left_begin()] = true;
                    ++num_covered;
                }
            }
            else {
#ifdef debug_apx_cover
                cerr << "enqueuing " << graph.get_graph().get_id(left_side) << " with " << left_edges.size() << " uncovered edges" << endl;
#endif
                queue.emplace_back(left_edges.size(), graph.left_iterator(left_side) - graph.left_begin());
                push_heap(queue.begin(), queue.end(), cmp);
            }
        }
    }
}

}
