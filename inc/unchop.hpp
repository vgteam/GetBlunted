#ifndef VG_ALGORITHMS_UNCHOP_HPP_INCLUDED
#define VG_ALGORITHMS_UNCHOP_HPP_INCLUDED


#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <handlegraph/util.hpp>

#include <vector>

namespace bluntifier {

/**
 * Unchop by gluing abutting handles with just a single edge between them and
 * compatible path steps together.
 */
void unchop(handlegraph::MutablePathDeletableHandleGraph* graph);

/**
 * Chop the graph so nodes are at most max_node_length
 */
void chop(handlegraph::MutableHandleGraph* graph, size_t max_node_length);


}

#endif
