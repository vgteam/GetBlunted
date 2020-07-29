#include "canonicalize.hpp"


pair<handle_t, handle_t> canonicalize(handle_t left, handle_t right, const HandleGraph& graph) {
    if (graph.get_id(left) < graph.get_id(right) or (graph.get_id(left) == graph.get_id(right) and graph.get_is_reverse(left))) {
        return make_pair(left, right);
    }
    else {
        return make_pair(graph.flip(right), graph.flip(left));
    }
}
