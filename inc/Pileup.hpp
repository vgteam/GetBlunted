#ifndef BLUNTIFIER_PILEUP_HPP
#define BLUNTIFIER_PILEUP_HPP

#include <handlegraph/handle_graph.hpp>
#include <vector>
#include <deque>

using handlegraph::handle_t;
using std::vector;
using std::deque;


namespace bluntifier{


class Pileup {
    /// Attributes ///

    // The data structure containing the multiple sequence alignment. Columns (arrays) correspond to coordinates in
    // alignment. Rows correspond to a single sequence's path through the alignment when adjusted for inserts/deletes
    deque <vector <char> > matrix;

    // In terms of column coordinates in the matrix, what path did the last/current alignment take?
    vector<size_t> path;

    // What was the last/current sequence that was added to the alignment pileup?
    handle_t sequence_id;

    // Has this path been invalidated since it was created (inserts added)? If so it needs to be recomputed
    bool path_valid;

    /// Methods ///


};

}

#endif //BLUNTIFIER_PILEUP_HPP
