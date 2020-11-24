#ifndef BLUNTIFIER_BLUNTIFIER_HPP
#define BLUNTIFIER_BLUNTIFIER_HPP

#include "AdjacencyComponent.hpp"
#include "BicliqueCover.hpp"
#include "Biclique.hpp"
#include "OverlapMap.hpp"
#include "Duplicator.hpp"
#include "Subgraph.hpp"
#include "align.hpp"
#include "OverlappingOverlap.hpp"
#include "OverlappingOverlapSplicer.hpp"
#include "gfa_to_handle.hpp"
#include "utility.hpp"
#include "unchop.hpp"

#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"

#include <unordered_map>

#include "bdsg/hash_graph.hpp"


using bluntifier::gfa_to_path_handle_graph_in_memory;
using bluntifier::gfa_to_path_handle_graph;
using bluntifier::gfa_to_handle_graph;
using bluntifier::handle_graph_to_gfa;
using bluntifier::parent_path;
using bluntifier::join_paths;
using bluntifier::IncrementalIdMap;
using bluntifier::OverlapMap;
using bluntifier::Duplicator;
using bluntifier::Alignment;
using bluntifier::for_each_adjacency_component;
using bluntifier::AdjacencyComponent;
using bluntifier::BipartiteGraph;
using bluntifier::BicliqueCover;
using bluntifier::BicliqueEdgeIndex;
using bluntifier::Bicliques;
using bluntifier::NodeInfo;
using bluntifier::Subgraph;
using bluntifier::OverlappingOverlapSplicer;
using bluntifier::OverlappingNodeInfo;
using bluntifier::OverlappingChild;
using bluntifier::bipartition;
using bluntifier::copy_path_handle_graph;
using bluntifier::duplicate_prefix;
using bluntifier::duplicate_suffix;
using bluntifier::run_command;
using bluntifier::unchop;

using handlegraph::MutablePathDeletableHandleGraph;
using handlegraph::as_integer;
using handlegraph::handle_t;
using bdsg::HashGraph;
using spoa::AlignmentEngine;
using spoa::AlignmentType;
using spoa::Graph;


namespace bluntifier {

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

    map <nid_t, nid_t> child_to_parent;
    map <nid_t, set<nid_t> > parent_to_children;

    vector <Subgraph> subgraphs;
    map <nid_t, OverlappingNodeInfo> overlapping_overlap_nodes;

    // Child node -> start_index -> (parent_node, stop_index)
    unordered_map<nid_t, map <nid_t, pair <size_t, size_t> > > provenance_map;

    unordered_set <handle_t> to_be_destroyed;

public:
    /// Methods ///
    Bluntifier(string gfa_path);

    void bluntify();

    void write_provenance(string& output_path);

    void write_to_gfa();

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

    bool is_oo_node_child(nid_t node_id);

    bool is_oo_node_parent(nid_t node_id);

    void find_path_info(Subgraph& subgraph, handle_t handle, PathInfo& path_info, string& path_name);

    void compute_provenance();

    void find_child_provenance(nid_t child_node, nid_t parent_node_id, size_t parent_index, bool side);

};


}

#endif //BLUNTIFIER_BLUNTIFIER_HPP
