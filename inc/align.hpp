#ifndef BLUNTIFIER_ALIGN_HPP
#define BLUNTIFIER_ALIGN_HPP

#include "handlegraph/handle_graph.hpp"
#include "OverlappingOverlap.hpp"
#include "Subgraph.hpp"
#include "Biclique.hpp"
#include "utility.hpp"
#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"

#include <unordered_map>

using handlegraph::HandleGraph;
using handlegraph::nid_t;
using spoa::AlignmentEngine;
using spoa::AlignmentType;
using spoa::Graph;

using std::unordered_map;
using std::unique_ptr;
using std::to_string;
using std::string;


namespace bluntifier {

void add_alignments_to_poa(
        const HandleGraph& gfa_graph,
        Subgraph& subgraph,
        Graph& spoa_graph,
        unique_ptr<AlignmentEngine>& alignment_engine,
        const vector<edge_t>& biclique);


void convert_spoa_to_bdsg(
        const HandleGraph& gfa_graph,
        Subgraph& subgraph,
        Graph& spoa_graph);


void align_biclique_overlaps(
        size_t i,
        const HandleGraph& gfa_graph,
        const Bicliques& bicliques,
        vector <Subgraph>& subgraphs);


void harmonize_biclique_orientations(HandleGraph& gfa_graph, Bicliques& bicliques);

}

#endif //BLUNTIFIER_ALIGN_HPP
