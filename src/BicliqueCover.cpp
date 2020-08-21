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
    
    size_t edge_count = 0;
    for (auto it = graph.left_begin(), end = graph.left_end(); it != end; ++it) {
        edge_count += graph.get_degree(*it);
    }
    // TODO: magic number (2^16)
    if (edge_count * (graph.left_size() + graph.right_size()) <= 65536) {
        // compute the biclique cover on the simplified graph using the algorithm
        // of Amilhastre, et al. 1998
        vector<pair<handle_t, vector<handle_t>>> simplifications;
        BipartiteGraph simplified = graph.simplify(simplifications);
        GaloisLattice galois_lattice(simplified);
        if (galois_lattice.is_domino_free()) {
            return_val = galois_lattice.biclique_separator();
            unsimplify(return_val, simplifications);
        }
    }
    
    if (return_val.empty()){
        // use the heuristic algorithm of Ene, et al. 2008
        return_val = heuristic_cover();
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

vector<bipartition> BicliqueCover::heuristic_cover() const {
    // TODO
    return vector<bipartition>();
}

}
