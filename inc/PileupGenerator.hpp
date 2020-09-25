#ifndef BLUNTIFIER_PILEUPGENERATOR_HPP
#define BLUNTIFIER_PILEUPGENERATOR_HPP

#include "IncrementalIdMap.hpp"
#include "BicliqueCover.hpp"
#include "OverlapMap.hpp"
#include "Pileup.hpp"
#include <vector>
#include <deque>
#include <stack>

using bluntifier::bipartition;
using handlegraph::handle_t;
using handlegraph::edge_t;

using std::vector;
using std::deque;
using std::queue;
using std::stack;


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



class AlignmentData{
public:
    uint64_t start;
    uint64_t stop;

    AlignmentData(uint64_t start, uint64_t stop);
};


class PileupGenerator {
public:
    /// Attributes ///

    /// Methods ///
    PileupGenerator();

    // Given a bipartition, build a multiple sequence alignment by projecting a set of pairwise alignments
    // through each other, using an arbitrary subset of pairs
    static void generate_from_bipartition(
            const BipartiteGraph& bipartite_graph,
            const IncrementalIdMap<string>& id_map,
            OverlapMap& overlaps,
            HandleGraph& graph,
            Pileup& pileup);

    // Basically just throw all the sequences into a POA alignment and see what happens
    static void generate_spoa_graph_from_bipartition(
            const BipartiteGraph& bipartite_graph,
            const IncrementalIdMap<string>& id_map,
            OverlapMap& overlaps,
            HandleGraph& graph,
            Pileup& pileup);

    // A generator-style DFS walk of the nodes in the partition
    static bool traverse_bipartition_nodes(
            const HandleGraph& graph,
            const OverlapMap& overlaps,
            const BipartiteGraph& bipartite_graph,
            BicliqueIterator& iterator);

    // A generator-style walk along the bipartition such that each step fills in the "edge" object with an edge
    // containing one sequence that was in a previous alignment allowing a projected MSA to be built.
    static bool traverse_bipartition_edges(
            const HandleGraph& graph,
            const OverlapMap& overlaps,
            const BipartiteGraph& bipartite_graph,
            BicliqueIterator& iterator);

    static void update_pseudoref(
            Pileup& pileup,
            HandleGraph& graph,
            handle_t pseudo_reference,
            size_t pseudo_ref_length,
            size_t prev_pseudo_ref_length,
            int64_t pseudo_ref_id,
            bool is_left);

private:
    static void debug_print(
            const IncrementalIdMap<string>& id_map,
            HandleGraph& graph,
            OverlapMap& overlaps,
            BicliqueIterator& biclique_iterator);

    static bool pseudoref_is_reversed(
            const Pileup& pileup,
            const edge_t& canonical_edge,
            const BicliqueIterator& biclique_iterator);
};


}

#endif //BLUNTIFIER_PILEUPGENERATOR_HPP
