#ifndef BLUNTIFIER_DUPLICATE_TERMINUS_HPP
#define BLUNTIFIER_DUPLICATE_TERMINUS_HPP

#include "handlegraph/mutable_path_mutable_handle_graph.hpp"
#include "handlegraph/handle_graph.hpp"
#include <deque>
#include <queue>

using handlegraph::MutablePathMutableHandleGraph;
using handlegraph::handle_t;
using std::deque;
using std::queue;


namespace bluntifier{


void duplicate_terminus(
        MutablePathMutableHandleGraph& graph,
        queue<size_t>& sizes,
        deque<handle_t>& children,
        handle_t parent_node);


}

#endif //BLUNTIFIER_DUPLICATE_TERMINUS_HPP
