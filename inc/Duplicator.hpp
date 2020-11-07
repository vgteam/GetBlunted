#ifndef BLUNTIFIER_DUPLICATOR_HPP
#define BLUNTIFIER_DUPLICATOR_HPP

#include "BicliqueCover.hpp"
#include "Biclique.hpp"
#include "NodeInfo.hpp"
#include "OverlapMap.hpp"
#include "OverlappingOverlap.hpp"
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

    map<nid_t, OverlappingNodeInfo>::iterator preprocess_overlapping_overlaps(
            MutablePathDeletableHandleGraph& gfa_graph,
            array <deque <size_t>, 2>& sorted_sizes_per_side,
            array <deque <size_t>, 2>& sorted_bicliques_per_side,
            array<map<size_t, handle_t>, 2>& biclique_side_to_child,
            const NodeInfo& node_info);

    void postprocess_overlapping_overlap(
            const HandleGraph& gfa_graph,
            map<nid_t, OverlappingNodeInfo>::iterator iter,
            array<map<size_t, handle_t>, 2> biclique_side_to_child);

    bool contains_overlapping_overlaps(
            const HandleGraph& gfa_graph,
            handle_t parent_handle,
            const array <deque <size_t>, 2>& sorted_sizes_per_side);
};


}

#endif //BLUNTIFIER_DUPLICATOR_HPP
