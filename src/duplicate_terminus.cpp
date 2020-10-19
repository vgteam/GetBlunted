#include "duplicate_terminus.hpp"

using std::runtime_error;


namespace bluntifier {


// Recursively duplicate the terminus of a node at the specified loci. Assume the node is in the correct orientation
// so that the desired terminus side is left
void duplicate_terminus(
        MutablePathMutableHandleGraph& graph,
        deque<size_t>& sizes,
        deque<handle_t>& children,
        handle_t parent_node) {

    auto size = sizes.front();
    sizes.pop_front();

    if (size < sizes.front()){
        throw runtime_error("ERROR: recursive duplicator only operates on sizes sorted in descending order");
    }

    auto fragments = graph.divide_handle(parent_node, size);

    // Store the right side of the node, which will not be further modified by successive recursions
    if (children.empty()) {
        children.emplace_front(fragments.second);
    }

    while (true) {
        // If this is the terminal recursion, then store the final node (and don't duplicate it)
        if (sizes.empty()) {
            children.emplace_front(fragments.first);
            return;
        }

        // Copy the left fragment and connect it to the right side
        auto dupe = graph.create_handle(graph.get_sequence(fragments.first));
        graph.create_edge(dupe, fragments.second);

        children.emplace_front(dupe);

        if (size == sizes.front()){
            size = sizes.front();
            sizes.pop_front();
        }
        else{
            break;
        }
    }

    // Pass on the leftovers (LEFT-overs! ha!) to be split and duplicated again
    duplicate_terminus(graph, sizes, children, fragments.first);
}


}
