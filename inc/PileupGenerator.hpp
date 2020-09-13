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

private:
    static void debug_print(
            const IncrementalIdMap<string>& id_map,
            HandleGraph& graph,
            OverlapMap& overlaps,
            BicliqueIterator& biclique_iterator);
};


}

#endif //BLUNTIFIER_PILEUPGENERATOR_HPP
