#ifndef BLUNTIFIER_NODEINFO_HPP
#define BLUNTIFIER_NODEINFO_HPP

#include "BicliqueCover.hpp"
#include "Biclique.hpp"
#include "OverlapMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "duplicate_terminus.hpp"
#include "copy_graph.hpp"
#include "utility.hpp"

#include "bdsg/hash_graph.hpp"

using handlegraph::MutablePathDeletableHandleGraph;
using handlegraph::as_integer;
using handlegraph::handle_t;
using bdsg::HashGraph;

namespace bluntifier {

class OverlapInfo {
public:
    size_t edge_index;
    size_t length;

    OverlapInfo(size_t edge_index, size_t length);
};


class NodeInfo {
public:
    array<map<size_t, vector<OverlapInfo> >, 2> factored_overlaps;
    const vector<vector<BicliqueEdgeIndex> >& node_to_biclique_edge;
    const Bicliques& bicliques;
    const HandleGraph& gfa_graph;
    const OverlapMap& overlaps;
    const nid_t node_id;

    NodeInfo(
            const vector<vector<BicliqueEdgeIndex> >& node_to_biclique_edge,
            const Bicliques& bicliques,
            const HandleGraph& gfa_graph,
            const OverlapMap& overlaps,
            nid_t node_id);

    NodeInfo(
            const vector<vector<BicliqueEdgeIndex> >& node_to_biclique_edge,
            const map <nid_t, nid_t>& child_to_parent,
            const Bicliques& bicliques,
            const HandleGraph& gfa_graph,
            const OverlapMap& overlaps,
            nid_t node_id);

    void factor_overlaps_by_biclique_and_side();
    void factor_overlaps_by_biclique_and_side(const map <nid_t, nid_t>& child_to_parent);

    void sort_factored_overlaps();

    size_t get_overlap_length(edge_t edge, bool side);

    void get_sorted_biclique_extents(
            array<deque<size_t>, 2>& sorted_extents_per_side,
            array<deque<size_t>, 2>& sorted_bicliques_per_side) const;

    void print_stats() const;
};


}

#endif //BLUNTIFIER_NODEINFO_HPP
