#ifndef BLUNTIFIER_DUPLICATOR_HPP
#define BLUNTIFIER_DUPLICATOR_HPP

#include "BicliqueCover.hpp"
#include "Biclique.hpp"
#include "NodeInfo.hpp"
#include "OverlapMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "duplicate_terminus.hpp"
#include "utility.hpp"

#include "bdsg/hash_graph.hpp"

using handlegraph::MutablePathDeletableHandleGraph;
using handlegraph::as_integer;
using handlegraph::handle_t;
using handlegraph::nid_t;
using bdsg::HashGraph;

namespace bluntifier{


class OverlappingNodeInfo{
public:
    /// Attributes ///

    // Assuming a forward oriented parent node, map each biclique (index) that this node participates in to each of its
    // duplicated segments, for each side
    array <map <size_t, handle_t> , 2> biclique_side_to_child;

    // What was the original GFA node id?
    nid_t parent_node;

    // What length was the original node?
    size_t length;

    /// Methods ///

    OverlappingNodeInfo(nid_t parent_node);

private:

};


class Duplicator{
public:
    /// Attributes ///
    const vector <vector <BicliqueEdgeIndex> >& node_to_biclique_edge;
    map<nid_t, OverlappingNodeInfo> overlapping_overlap_nodes;
    Bicliques& bicliques;
    OverlapMap& overlaps;
    map <nid_t, set<nid_t> > parent_to_children;
    map <nid_t, nid_t> child_to_parent;


    /// Methods ///
    Duplicator(
            const vector <vector <BicliqueEdgeIndex> >& node_to_biclique_edge,
            Bicliques& bicliques,
            OverlapMap& overlaps);

    void duplicate_all_node_termini(MutablePathDeletableHandleGraph& gfa_graph);

private:
    void repair_edges(
            MutablePathDeletableHandleGraph& gfa_graph,
            const array <map <size_t, handle_t>, 2>& biclique_side_to_child,
            handle_t old_handle,
            handle_t old_handle_flipped);

    void remove_participating_edges(
            const array <deque <size_t>, 2>& sorted_bicliques_per_side,
            MutablePathDeletableHandleGraph& gfa_graph,
            nid_t parent_node);

    void duplicate_termini(
            MutablePathDeletableHandleGraph& gfa_graph,
            array <deque <size_t>, 2>& sorted_sizes_per_side,
            const array <deque <size_t>, 2>& sorted_bicliques_per_side,
            array<map<size_t, handle_t>, 2>& biclique_side_to_child,
            const NodeInfo& node_info);

    void duplicate_overlapping_termini(
            MutablePathDeletableHandleGraph& gfa_graph,
            const array <deque <size_t>, 2>& sorted_sizes_per_side,
            const array <deque <size_t>, 2>& sorted_bicliques_per_side,
            array<map<size_t, handle_t>, 2>& biclique_side_to_child,
            const NodeInfo& node_info);

};


}

#endif //BLUNTIFIER_DUPLICATOR_HPP
