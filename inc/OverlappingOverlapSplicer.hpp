
#ifndef BLUNTIFIER_OVERLAPPINGOVERLAPSPLICER_HPP
#define BLUNTIFIER_OVERLAPPINGOVERLAPSPLICER_HPP

#include "handlegraph/mutable_path_deletable_handle_graph.hpp"
#include "handlegraph/handle_graph.hpp"
#include "OverlappingOverlap.hpp"
#include "Subgraph.hpp"
#include "utility.hpp"

using std::vector;
using std::string;
using std::runtime_error;
using std::to_string;

using handlegraph::MutablePathDeletableHandleGraph;
using handlegraph::HandleGraph;
using handlegraph::handle_t;
using handlegraph::nid_t;

namespace bluntifier {

class OverlappingOverlapSplicer {
public:
    map<nid_t, OverlappingNodeInfo>& overlapping_overlap_nodes;
    const vector<Subgraph>& subgraphs;

    OverlappingOverlapSplicer(
            map<nid_t, OverlappingNodeInfo>& overlapping_overlap_nodes,
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
};


void find_path_info(
        const vector<Subgraph>& subgraphs,
        const HandleGraph& gfa_graph,
        size_t biclique_index,
        handle_t handle,
        PathInfo& path_info,
        string& path_name);

void splice_overlapping_overlaps(
        MutablePathDeletableHandleGraph& gfa_graph,
        vector <Subgraph>& subgraphs,
        map<nid_t, OverlappingNodeInfo>& overlapping_overlap_nodes);


}

#endif //BLUNTIFIER_OVERLAPPINGOVERLAPSPLICER_HPP
