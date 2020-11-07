#include "OverlappingOverlapSplicer.hpp"
#include "handlegraph/util.hpp"


using std::cout;

namespace bluntifier {


OverlappingOverlapSplicer::OverlappingOverlapSplicer(
        map<nid_t, OverlappingNodeInfo>& overlapping_overlap_nodes,
        const vector<Subgraph>& subgraphs):
    overlapping_overlap_nodes(overlapping_overlap_nodes),
    subgraphs(subgraphs)
{}


void OverlappingOverlapSplicer::find_path_info(
        const HandleGraph& gfa_graph,
        size_t biclique_index,
        handle_t handle,
        PathInfo& path_info,
        string& path_name){

    auto& subgraph = subgraphs[biclique_index];

    // Don't know which side of the biclique this overlap was on until we search for it in the subgraph
    auto result = subgraph.paths_per_handle[0].find(handle);
    if (result == subgraph.paths_per_handle[0].end()) {

        result = subgraph.paths_per_handle[1].find(handle);

        // Sanity check
        if (result == subgraph.paths_per_handle[1].end()) {
            throw runtime_error("ERROR: node not found in biclique subgraph: " +
                                to_string(gfa_graph.get_id(handle)));
        }
    }

    path_name = subgraph.graph.get_path_name(result->second.path_handle);
    path_info = result->second;
}


void OverlappingOverlapSplicer::splice_overlapping_overlaps(MutablePathDeletableHandleGraph& gfa_graph) {

    for (auto& item: overlapping_overlap_nodes) {
        auto node_id = item.first;
        auto parent_handle = gfa_graph.get_handle(node_id, false);
        auto& overlap_info = item.second;
        auto parent_size = overlap_info.length;

        cout << "Splicing overlapping overlap: " << node_id << '\n';


            // Splice into the remaining body of the node if any of the other overlaps don't overlap this one
    }
}



}
