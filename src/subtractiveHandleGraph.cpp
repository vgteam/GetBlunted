/**
 * \file subtractiveHandleGraph.cpp: contains the implementation of SubtractiveHandleGraph
 */


#include "subtractiveHandleGraph.hpp"

#include <atomic>

namespace bluntifier {

using namespace std;

SubtractiveHandleGraph::SubtractiveHandleGraph(const HandleGraph& super) : super(&super) {
        // nothing to do
    }

void SubtractiveHandleGraph::subtract_edge(const handle_t& left, const handle_t& right) {
    subtracted_edges.insert(super->edge_handle(left, right));
}

bool SubtractiveHandleGraph::has_node(nid_t node_id) const {
    return super->has_node(node_id);
}

handle_t SubtractiveHandleGraph::get_handle(const nid_t& node_id, bool is_reverse) const {
    return super->get_handle(node_id, is_reverse);
}

nid_t SubtractiveHandleGraph::get_id(const handle_t& handle) const {
    return super->get_id(handle);
}

bool SubtractiveHandleGraph::get_is_reverse(const handle_t& handle) const {
    return super->get_is_reverse(handle);
}

handle_t SubtractiveHandleGraph::flip(const handle_t& handle) const {
    return super->flip(handle);
}

size_t SubtractiveHandleGraph::get_length(const handle_t& handle) const {
    return super->get_length(handle);
}

string SubtractiveHandleGraph::get_sequence(const handle_t& handle) const {
    return super->get_sequence(handle);
}

bool SubtractiveHandleGraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                               const function<bool(const handle_t&)>& iteratee) const {
    // only let it travel along edges whose endpoints are in the subgraph
    return super->follow_edges(handle, go_left, [&](const handle_t& next) {
        edge_t edge = go_left ? super->edge_handle(next, handle) : super->edge_handle(handle, next);
        if (subtracted_edges.count(edge)) {
            // this edge is removed, iterate past it
            return true;
        }
        else {
            // let the lambda pass through to the underlying graph
            return iteratee(next);
        }
    });
}

bool SubtractiveHandleGraph::for_each_handle_impl(const function<bool(const handle_t&)>& iteratee,
                                                  bool parallel) const {
    super->for_each_handle(iteratee, parallel);
}

size_t SubtractiveHandleGraph::get_node_count() const {
    return super->get_node_count();
}

nid_t SubtractiveHandleGraph::min_node_id() const {
    return super->min_node_id();
}

nid_t SubtractiveHandleGraph::max_node_id() const {
    return super->max_node_id();
}

handle_t SubtractiveHandleGraph::get_underlying_handle(const handle_t& handle) const {
    return handle;
}
}

