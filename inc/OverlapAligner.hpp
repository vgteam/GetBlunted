#ifndef BLUNTIFIER_OVERLAPALIGNER_HPP
#define BLUNTIFIER_OVERLAPALIGNER_HPP

#include "handlegraph/handle_graph.hpp"
#include "Subgraph.hpp"
#include "Biclique.hpp"
#include "utility.hpp"
#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"

#include <unordered_map>

using handlegraph::HandleGraph;
using spoa::AlignmentEngine;
using spoa::AlignmentType;
using spoa::Graph;

using std::unordered_map;
using std::unique_ptr;
using std::to_string;
using std::string;


namespace bluntifier {

class OverlapAligner {
public:

    OverlapAligner();
};


void harmonize_biclique_orientations(HandleGraph& gfa_graph, Bicliques& bicliques);

}

#endif //BLUNTIFIER_OVERLAPALIGNER_HPP
