#ifndef BLUNTIFIER_PILEUP_HPP
#define BLUNTIFIER_PILEUP_HPP

#include "handlegraph/handle_graph.hpp"
#include "bdsg/packed_graph.hpp"
#include "IncrementalIdMap.hpp"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <deque>
#include <string>

using bluntifier::IncrementalIdMap;
using handlegraph::handle_t;
using handlegraph::edge_t;
using bdsg::PackedGraph;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::string;
using std::deque;


namespace bluntifier{


class Pileup {
public:
    /// Attributes ///

    // The data structure containing the multiple sequence alignment. Columns (arrays) correspond to coordinates in
    // alignment. Rows correspond to a single sequence's path through the alignment when adjusted for inserts/deletes
    deque <vector <char> > matrix;

    // A new graph to represent the multiple sequence alignment
    PackedGraph graph;

    // In terms of nodes in the graph, what path does each alignment take?
    vector <deque <handle_t> > paths;

    // List of edges that have been traversed in order of traversal
    vector <edge_t> edges_traversed;

    // Convert nodes to integer IDs
    IncrementalIdMap<handle_t> id_map;

    // For each node in the original gfa graph (stored here in order of their appearance) what are the nodes in the
    // pileup graph (from left to right) that will be spliced back into the gfa
    vector <vector <handle_t> > splice_nodes;

    // Filler character
    static const char space;

    /// Methods ///
    Pileup();
    void to_string(string& s);

};

}

#endif //BLUNTIFIER_PILEUP_HPP
