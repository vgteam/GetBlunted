#ifndef BLUNTIFIER_PILEUPGENERATOR_HPP
#define BLUNTIFIER_PILEUPGENERATOR_HPP

#include "BicliqueCover.hpp"
#include "OverlapMap.hpp"
#include "Pileup.hpp"
#include <vector>
#include <deque>

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
    /// Attributes ///
    edge_t& edge;
    unordered_set<handle_t> visited;
    stack<handle_t> nodes;
    bool first_step;

    /// Methods ///
    void update(edge_t edge);
    void next();
};


class PileupGenerator {
public:
    /// Attributes ///

    /// Methods ///
    PileupGenerator();

    // Given a bipartition, build a multiple sequence alignment by projecting a set of pairwise alignments
    // through each other, using an arbitrary subset of pairs
    void generate_from_bipartition(bipartition& partition, HandleGraph& graph, Pileup& pileup);

    // A generator-style function which returns false when all nodes are visited via some edge.
    // The goal is to walk along the bipartition such that each step fills in the "edge" object
    // with an edge containing one sequence that was in a previous alignment.
    // This allows a projected MSA to be built.
    bool traverse_bipartition(HandleGraph& graph, bipartition& partition, BicliqueIterator& iterator);


};


}

#endif //BLUNTIFIER_PILEUPGENERATOR_HPP
