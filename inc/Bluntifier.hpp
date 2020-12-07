#ifndef BLUNTIFIER_BLUNTIFIER_HPP
#define BLUNTIFIER_BLUNTIFIER_HPP

#include "AdjacencyComponent.hpp"
#include "BicliqueCover.hpp"
#include "Biclique.hpp"
#include "OverlapMap.hpp"
#include "Duplicator.hpp"
#include "Subgraph.hpp"
#include "OverlappingOverlap.hpp"
#include "OverlappingOverlapSplicer.hpp"
#include "gfa_to_handle.hpp"
#include "utility.hpp"
#include "unchop.hpp"

#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"

#include <unordered_map>

#include "bdsg/hash_graph.hpp"


using handlegraph::MutablePathDeletableHandleGraph;
using handlegraph::as_integer;
using handlegraph::handle_t;
using bdsg::HashGraph;
using spoa::AlignmentEngine;
using spoa::AlignmentType;
using spoa::Graph;


namespace bluntifier {

class ProvenanceInfo{
public:
    size_t start;
    size_t stop;
    bool reversal;

    ProvenanceInfo(size_t start, size_t stop, bool reversal);
};


class Bluntifier {
private:
    /// Attributes ///
    string gfa_path;

    HashGraph gfa_graph;
    IncrementalIdMap<string> id_map;
    OverlapMap overlaps;

    // Where all the ACs go
    vector<AdjacencyComponent> adjacency_components;

    Bicliques bicliques;
    mutex biclique_mutex;
    vector <vector <BicliqueEdgeIndex> > node_to_biclique_edge;

    // Store the mapping from children to parent, and a boolean to tell whether that child is a suffix/prefix or the
    // original node material
    map <nid_t, pair<nid_t, bool> > child_to_parent;
    map <nid_t, set<nid_t> > parent_to_children;

    vector <Subgraph> subgraphs;
    map <nid_t, OverlappingNodeInfo> overlapping_overlap_nodes;

    // Child node -> start_index -> (parent_node, stop_index)
    unordered_map<nid_t, multimap <nid_t, ProvenanceInfo> > provenance_map;

    unordered_set <nid_t> to_be_destroyed;

public:
    /// Methods ///
    Bluntifier(string gfa_path);

    void bluntify();

    void write_provenance(string& output_path);

private:
    void deduplicate_and_canonicalize_biclique_cover(
            vector<bipartition>& biclique_cover,
            vector<vector<edge_t> >& deduplicated_biclique_cover);

    void compute_biclique_cover(size_t i);

    void print_adjacency_components_stats(size_t i);

    void map_splice_sites_by_node();

    void harmonize_biclique_orientations();

    void align_biclique_overlaps(size_t i);

    void add_alignments_to_poa(
            Graph& spoa_graph,
            unique_ptr<AlignmentEngine>& alignment_engine,
            size_t i);

    void convert_spoa_to_bdsg(Graph& spoa_graph, size_t i);

    void splice_subgraphs();

    tuple <bool, bool> is_oo_node(nid_t node_id);

    void compute_provenance();

    // Lol
    void update_path_provenances(
            nid_t parent_node_id,
            size_t parent_index,
            bool parent_side,
            bool reversal,
            size_t parent_length,
            OverlapInfo& overlap_info,
            edge_t& canonical_edge,
            edge_t& edge,
            nid_t child_id);
};


}

#endif //BLUNTIFIER_BLUNTIFIER_HPP
