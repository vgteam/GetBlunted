#include "duplicate_terminus.hpp"

using std::runtime_error;
using std::to_string;
using std::pair;


namespace bluntifier {


// Recursively duplicate the terminus of a node at the specified loci. Assume the node is in the correct orientation
// so that the desired terminus side is left
void duplicate_prefix(
        MutablePathMutableHandleGraph& graph,
        deque<size_t>& sizes,
        deque<handle_t>& children,
        handle_t parent_handle) {

    size_t size = sizes.front();
    sizes.pop_front();

    if (not sizes.empty() and size < sizes.front()){
        throw runtime_error("ERROR: recursive duplicator only operates on sizes sorted in descending order");
    }

    if (size == 0){
        std::cerr << "WARNING: 0 length overlap is ignored in duplicator\n";

        // Store the right side of the node, which will not be further modified by successive recursions
        if (children.empty()) {
            children.emplace_back(parent_handle);
        }

        // Make a placeholder so that the number of children is still sizes.size() + 1
        children.emplace_back(parent_handle);

        // Overlaps are sorted so this could only be the last overlap. So just end here.
        return;
    }

    pair<handle_t,handle_t> fragments;

    if (size < graph.get_length(parent_handle)) {
        fragments = graph.divide_handle(parent_handle, size);
    }
    else if (size == graph.get_length(parent_handle)) {
        fragments = {parent_handle, parent_handle};
    }
    else {
        throw runtime_error("ERROR: cannot duplicate overlap of length " + to_string(size) +
                            " which is longer than parent node (id): " + to_string(graph.get_id(parent_handle)));
    }

    // Store the right side of the node, which will not be further modified by successive recursions
    if (children.empty()) {
        children.emplace_back(fragments.second);
    }

    while (true) {
        // If this is the terminal recursion, then store the final node (and don't duplicate it)
        if (sizes.empty()) {
            children.emplace_back(fragments.first);
            return;
        }

        // Copy the left fragment and connect it to the right side
        auto dupe = graph.create_handle(graph.get_sequence(fragments.first));
        graph.create_edge(dupe, fragments.second);

        children.emplace_back(dupe);

        if (size == sizes.front()){
            size = sizes.front();
            sizes.pop_front();
        }
        else{
            break;
        }
    }

    // Pass on the leftovers (LEFT-overs! ha!) to be split and duplicated again
    duplicate_prefix(graph, sizes, children, fragments.first);
}


// Recursively duplicate the terminus of a node at the specified loci. Assume the node is in the correct orientation
// so that the desired terminus side is left
void duplicate_suffix(
        MutablePathMutableHandleGraph& graph,
        deque<size_t>& sizes,
        deque<handle_t>& children,
        handle_t parent_handle) {

    size_t size = sizes.front();
    sizes.pop_front();

    if (not sizes.empty() and size < sizes.front()){
        throw runtime_error("ERROR: recursive duplicator only operates on sizes sorted in descending order");
    }

    if (size == 0){
        std::cerr << "WARNING: 0 length overlap is ignored in duplicator\n";

        // Store the right side of the node, which will not be further modified by successive recursions
        if (children.empty()) {
            children.emplace_back(parent_handle);
        }

        // Make a placeholder so that the number of children is still sizes.size() + 1
        children.emplace_back(parent_handle);

        // Overlaps are sorted so this could only be the last overlap. So just end here.
        return;
    }

    pair<handle_t,handle_t> fragments;

    if (size < graph.get_length(parent_handle)) {
        size_t coordinate = graph.get_length(parent_handle) - size;
        fragments = graph.divide_handle(parent_handle, {coordinate});
    }
    else if (size == graph.get_length(parent_handle)) {
        fragments = {parent_handle, parent_handle};
    }
    else {
        throw runtime_error("ERROR: cannot duplicate overlap of length " + to_string(size) +
                            " which is longer than parent node (id): " + to_string(graph.get_id(parent_handle)));
    }

    // Store the right side of the node, which will not be further modified by successive recursions
    if (children.empty()) {
        children.emplace_back(fragments.first);
    }

    while (true) {
        // If this is the terminal recursion, then store the final node (and don't duplicate it)
        if (sizes.empty()) {
            children.emplace_back(fragments.second);
            return;
        }

        // Copy the Right fragment and connect it to the Left side
        auto dupe = graph.create_handle(graph.get_sequence(fragments.second));
        graph.create_edge(fragments.first, dupe);

        children.emplace_back(dupe);

        if (size == sizes.front()){
            size = sizes.front();
            sizes.pop_front();
        }
        else{
            break;
        }
    }

    // Pass on the leftovers (LEFT-overs! ha!) to be split and duplicated again
    duplicate_suffix(graph, sizes, children, fragments.second);
}


}
