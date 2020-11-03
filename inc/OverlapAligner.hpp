#ifndef BLUNTIFIER_OVERLAPALIGNER_HPP
#define BLUNTIFIER_OVERLAPALIGNER_HPP

#include "handlegraph/handle_graph.hpp"
#include "Biclique.hpp"

using handlegraph::HandleGraph;

namespace bluntifier {

class OverlapAligner {
public:

    OverlapAligner();
};


void harmonize_biclique_orientations(HandleGraph& gfa_graph, Bicliques& bicliques);

}

#endif //BLUNTIFIER_OVERLAPALIGNER_HPP
