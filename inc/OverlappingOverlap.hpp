#ifndef BLUNTIFIER_OVERLAPPINGOVERLAP_HPP
#define BLUNTIFIER_OVERLAPPINGOVERLAP_HPP

#include "handlegraph/handle_graph.hpp"

#include <map>

using handlegraph::HandleGraph;
using handlegraph::path_handle_t;
using handlegraph::handle_t;
using handlegraph::nid_t;

using std::map;


namespace bluntifier{

class OverlappingChild{
public:
    const handle_t handle;
    const size_t biclique_index;
    const bool side;

    OverlappingChild(handle_t handle, size_t biclique_index, bool side);
    void print(HandleGraph& gfa_graph);
};

class OverlappingNodeInfo{
public:
    /// Attributes ///

    // Assuming a forward oriented parent node, map each child of this node by its termination index in the parent node
    map <size_t, OverlappingChild> overlapping_children;
    array<map<size_t, handle_t>, 2> biclique_side_to_child;

    // What was the original GFA node id?
    nid_t parent_node;

    // What length was the original node?
    size_t length;

    /// Methods ///

    OverlappingNodeInfo(nid_t parent_node);

private:

};
}

#endif //BLUNTIFIER_OVERLAPPINGOVERLAP_HPP
