
#ifndef BLUNTIFIER_OVERLAPPINGOVERLAPSPLICER_HPP
#define BLUNTIFIER_OVERLAPPINGOVERLAPSPLICER_HPP

#include "handlegraph/mutable_path_deletable_handle_graph.hpp"
#include "handlegraph/handle_graph.hpp"
#include "OverlappingOverlap.hpp"
#include "Subgraph.hpp"
#include "utility.hpp"
#include <utility>
#include <set>

using std::runtime_error;
using std::to_string;
using std::string;
using std::vector;
using std::pair;
using std::set;

using handlegraph::MutablePathDeletableHandleGraph;
using handlegraph::HandleGraph;
using handlegraph::handle_t;
using handlegraph::nid_t;

namespace bluntifier {


class OverlappingSplicePair{
public:
    size_t left_parent_index;
    size_t left_child_index;
    size_t left_length;
    string left_child_path_name;
    size_t right_parent_index;
    size_t right_child_index;
    size_t right_length;
    string right_child_path_name;

    OverlappingSplicePair()=default;
};


class OverlappingOverlapSplicer {
public:
    map <nid_t, OverlappingNodeInfo>& overlapping_overlap_nodes;
    map <nid_t, set<nid_t> >& parent_to_children;
    const vector<Subgraph>& subgraphs;

    OverlappingOverlapSplicer(
            map<nid_t, OverlappingNodeInfo>& overlapping_overlap_nodes,
            map <nid_t, set<nid_t> >& parent_to_children,
            const vector<Subgraph>& subgraphs);

    void splice_overlapping_overlaps(
            MutablePathDeletableHandleGraph& gfa_graph);

private:
    void find_path_info(
            const HandleGraph& gfa_graph,
            size_t biclique_index,
            handle_t handle,
            PathInfo& path_info,
            string& path_name);

    void find_parent_path_bounds(
            MutablePathDeletableHandleGraph& gfa_graph,
            nid_t parent_id,
            pair<size_t, size_t>& bounds);

    void find_splice_pairs(
            HandleGraph& gfa_graph,
            OverlappingNodeInfo& overlap_info,
            vector <OverlappingSplicePair>& oo_splice_pairs);

    // Iterate along a path until the cumulative number of bases iterated would equal "target_base_index" and then
    // return the handle and intra-handle index
    pair<handle_t, size_t> seek_to_path_base(
            MutablePathDeletableHandleGraph& gfa_graph,
            string& path_name,
            size_t target_base_index);

    // Convenience wrapper for seek_to_path_base
    pair<handle_t, size_t> seek_to_child_path_base(
            MutablePathDeletableHandleGraph& gfa_graph,
            OverlappingChild& overlapping_child,
            size_t target_base_index);

};

}

#endif //BLUNTIFIER_OVERLAPPINGOVERLAPSPLICER_HPP
