/**
 * \file BicliqueCover.cpp
 *
 * Implements algorithm for computing the biclique cover of a bipartite graph.
 */
#include "BicliqueCover.hpp"

//#define debug_galois_tree
//#define debug_galois_lattice
//#define debug_max_flow
 
namespace bluntifier {

using std::unordered_set;
using std::sort;
using std::deque;

BicliqueCover::BicliqueCover(const BipartiteGraph& graph) : graph(graph) {

}

BicliqueCover::~BicliqueCover() {
    
}

vector<bipartition> BicliqueCover::get() const {
    
    vector<bipartition> return_val;

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
            // choose the smallest solution as the final result
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
    // TODO
}

}
