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

    // Coverage in terms of sequences that contributed to each node in the alignment graph
    unordered_map <handle_t, unordered_set <handle_t> > coverage;

    // Convert nodes to integer IDs
    IncrementalIdMap<handle_t> id_map;

    // Filler character
    static const char space;

    /// Methods ///
    void to_string(string& s);

};

}

#endif //BLUNTIFIER_PILEUP_HPP
