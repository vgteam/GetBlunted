#ifndef BLUNTIFIER_PILEUPGENERATOR_HPP
#define BLUNTIFIER_PILEUPGENERATOR_HPP

#include "IncrementalIdMap.hpp"
#include "BicliqueCover.hpp"
#include "OverlapMap.hpp"
#include "Pileup.hpp"
#include <vector>
#include <deque>
#include <stack>
#include <mutex>

#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"

using spoa::AlignmentEngine;
using spoa::AlignmentType;
using spoa::Graph;

using bluntifier::bipartition;
using handlegraph::handle_t;
using handlegraph::edge_t;

using std::vector;
using std::deque;
using std::queue;
using std::stack;
using std::mutex;


namespace bluntifier {


class BicliqueIterator {
public:
    /// Node attributes ///
    handle_t node;
    unordered_set<handle_t> visited;
    stack<handle_t> node_stack;
    bool first_step;
    bool is_left;

    /// Edge attributes ///
    edge_t edge;
    handle_t prev_left_node;
    handle_t prev_right_node;
    bool prev_left_is_valid = false;
    bool prev_right_is_valid = false;

    /// Methods ///
    BicliqueIterator();
};


class PoaPileupGenerator {
public:
    /// Attributes ///
    const BipartiteGraph& bipartite_graph;
    const IncrementalIdMap<string>& id_map;
    OverlapMap& overlaps;
    HandleGraph& gfa_graph;
    PoaPileup& pileup;
    vector <vector <SpliceData> >& splice_sites;
    vector <mutex>& splice_site_mutexes;

    /// Methods ///
    PoaPileupGenerator(
            const BipartiteGraph& bipartite_graph,
            const IncrementalIdMap<string>& id_map,
            OverlapMap& overlaps,
            HandleGraph& gfa_graph,
            PoaPileup& pileup,
            vector <vector <SpliceData> >& splice_sites,
            vector <mutex>& splice_site_mutexes
    );

    // Basically just throw all the sequences into a POA alignment and see what happens
    void generate_from_bipartition(const bipartition& bipartition, size_t component_index);

    // Do POA with spoa for an arbitrary collection of edges
    void generate_from_edges(const vector<edge_t>& edges, size_t component_index);

private:
    void sort_alignment_data_by_length();

    SpliceData& get_splice_data(bool is_left, size_t id, size_t i);

    void add_alignments_to_poa(Graph& spoa_graph, unique_ptr<AlignmentEngine>& alignment_engine);

    void convert_spoa_to_bdsg(Graph& spoa_graph);
};


class PileupGenerator {
public:
    /// Attributes ///
    const BipartiteGraph& bipartite_graph;
    const IncrementalIdMap<string>& id_map;
    OverlapMap& overlaps;
    HandleGraph& gfa_graph;
    Pileup& pileup;
    vector <vector <SpliceData> >& splice_sites;
    vector <mutex>& splice_site_mutexes;

    /// Methods ///
    PileupGenerator(
            const BipartiteGraph& bipartite_graph,
            const IncrementalIdMap<string>& id_map,
            OverlapMap& overlaps,
            HandleGraph& gfa_graph,
            Pileup& pileup,
            vector <vector <SpliceData> >& splice_sites,
            vector <mutex>& splice_site_mutexes
    );

    // Given a bipartition, build a multiple sequence alignment by projecting a set of pairwise alignments
    // through each other, using an arbitrary subset of pairs
    void generate_from_bipartition();

    // A generator-style DFS walk of the nodes in the partition
    bool traverse_bipartition_nodes(BicliqueIterator& iterator);

    // A generator-style walk along the bipartition such that each step fills in the "edge" object with an edge
    // containing one sequence that was in a previous alignment allowing a projected MSA to be built.
    bool traverse_bipartition_edges(BicliqueIterator& iterator);

    void update_pseudoref(
            handle_t pseudo_reference,
            size_t pseudo_ref_length,
            size_t prev_pseudo_ref_length,
            int64_t pseudo_ref_id,
            bool is_left);

private:
    void debug_print(BicliqueIterator& biclique_iterator);

    bool pseudoref_is_reversed(
            const edge_t& canonical_edge,
            const BicliqueIterator& biclique_iterator);
};


}

#endif //BLUNTIFIER_PILEUPGENERATOR_HPP
