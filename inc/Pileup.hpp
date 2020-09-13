#ifndef BLUNTIFIER_PILEUP_HPP
#define BLUNTIFIER_PILEUP_HPP

#include "handlegraph/handle_graph.hpp"
#include "bdsg/packed_graph.hpp"
#include <vector>
#include <deque>
#include <string>

using handlegraph::handle_t;
using bdsg::PackedGraph;
using std::vector;
using std::deque;
using std::string;


namespace bluntifier{


class Pileup {
public:
    /// Attributes ///

    // The data structure containing the multiple sequence alignment. Columns (arrays) correspond to coordinates in
    // alignment. Rows correspond to a single sequence's path through the alignment when adjusted for inserts/deletes
    deque <vector <char> > matrix;

    //
    PackedGraph graph;

    // In terms of column coordinates in the matrix, what path each alignment take?
    vector <deque <handle_t> > paths;

    vector <handle_t> nodes;

    // Filler character
    static const char space;

    /// Methods ///
    void to_string(string& s);

};

}

#endif //BLUNTIFIER_PILEUP_HPP
