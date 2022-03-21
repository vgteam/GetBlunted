
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
using std::tuple;
using std::pair;
using std::set;

using handlegraph::MutablePathDeletableHandleGraph;
using handlegraph::HandleGraph;
using handlegraph::handle_t;
using handlegraph::nid_t;

namespace bluntifier {


class OverlappingSplicePair{
public:
    int64_t left_parent_index;
    int64_t left_child_index;
    string left_child_path_name;
    int64_t right_parent_index;
    int64_t right_child_index;
    string right_child_path_name;

    // From which side of the OO node is this pair? 0 = left, 1 = right. This is needed for full-node overlaps
    bool side;

    bool left_reversal;
    bool right_reversal;

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
    bool find_path_info(
            const HandleGraph& gfa_graph,
            size_t biclique_index,
            handle_t handle,
            PathInfo& path_info,
            string& path_name);

    void find_splice_pairs(
            HandleGraph& gfa_graph,
            OverlappingNodeInfo& overlap_info,
            vector <OverlappingSplicePair>& oo_splice_pairs);

    // Iterate along a path until the cumulative number of bases iterated would equal "target_base_index" and then
    // return the handle and intra-handle index
    tuple<handle_t, int64_t, int64_t, bool> seek_to_path_base(
            MutablePathDeletableHandleGraph& gfa_graph,
            string& path_name,
            int64_t target_base_index);

    tuple<handle_t, int64_t, int64_t, bool> seek_to_reverse_path_base(
            MutablePathDeletableHandleGraph& gfa_graph,
            string& path_name,
            int64_t target_base_index);

};

}

#endif //BLUNTIFIER_OVERLAPPINGOVERLAPSPLICER_HPP
