#ifndef BLUNTIFIER_APPLY_BULK_MODIFICATIONS_HPP_INCLUDED
#define BLUNTIFIER_APPLY_BULK_MODIFICATIONS_HPP_INCLUDED

/**
 * \file apply_bulk_modifications.hpp
 *
 * Defines utility algorithms for applying mutable graph operations in bulk.
 */

#include "handlegraph/handle_graph.hpp"

#include <vector>
#include <unordered_set>

#include "handlegraph/mutable_handle_graph.hpp"
#include "handlegraph/handle_graph.hpp"

namespace bluntifier{

using handlegraph::MutableHandleGraph;
using handlegraph::handle_t;
using std::unordered_set;
using std::vector;

/// Modifies underlying graph so that any node whose handle is given in the reverse orientation
/// is flipped so that all locally forward orientations match the orientation of the provided handles.
/// Returns a set of IDs for nodes that were flipped. Invalid if vector contains multiple handles to
/// the same node. May change the ordering of the underlying graph.
unordered_set<id_t> apply_orientations(MutableHandleGraph& graph, const vector<handle_t>& orientations);

}


#endif
