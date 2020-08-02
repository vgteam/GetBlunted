/**
 * \file biclique_cover.cpp
 *
 * Implements algorithm for computing the biclique cover of a bipartite graph.
 */
#include "biclique_cover.hpp"

namespace bluntifier {

BicliqueCover::BicliqueCover(const HandleGraph& graph,
                             const bipartition& partition) {
    
}

BicliqueCover::~BicliqueCover() {
    
}

vector<bipartition> BicliqueCover::get() const {
    
}


bool BicliqueCover::is_domino_free() const {
    
}

SubtractiveHandleGraph BicliqueCover::simplify() const {
    
    SubtractiveGraph sim
    
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

}
