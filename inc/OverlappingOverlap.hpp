#ifndef BLUNTIFIER_OVERLAPPINGOVERLAP_HPP
#define BLUNTIFIER_OVERLAPPINGOVERLAP_HPP

#include "handlegraph/handle_graph.hpp"

#include <vector>
#include <deque>
#include <array>
#include <map>

using handlegraph::HandleGraph;
using handlegraph::path_handle_t;
using handlegraph::handle_t;
using handlegraph::nid_t;

using std::vector;
using std::deque;
using std::array;
using std::map;


namespace bluntifier{

class OverlappingChild{
public:
    const handle_t handle;
    const size_t biclique_index;

    OverlappingChild(handle_t handle, size_t biclique_index);
    void print(HandleGraph& gfa_graph) const;
};


class OverlappingNodeInfo{
public:
    /// Attributes ///

    // Assuming a forward oriented parent node, map each child of this node by its termination index in the parent node
    array <map <size_t, OverlappingChild>, 2> overlapping_children;
    array <map <size_t, OverlappingChild>, 2> normal_children;

    // If there is anything leftover of the original node, it goes here
    vector <handle_t> leftover_parent;

    // What was the original GFA node id?
    nid_t parent_node;

    // What length was the original node?
    size_t length;

    /// Methods ///

    OverlappingNodeInfo(nid_t parent_node);
    void print(HandleGraph& graph);

    // Return a vector containing all the children whose extents exceed or match the given index
    void find_overlapping_children_by_base_index(size_t index, bool side, deque<OverlappingChild>& children);
    void find_normal_children_by_base_index(size_t index, bool side, deque<OverlappingChild>& children);

private:

};
}

#endif //BLUNTIFIER_OVERLAPPINGOVERLAP_HPP
