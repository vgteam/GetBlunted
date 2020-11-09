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


// TODO: if this ever becomes rate limiting, this whole loop could be replaced by a stepwise iterator with a "next"
// method. Then the base indexes would be iterated, and for all overlaps, the iterators are advanced in sync until they
// expire
pair<handle_t, size_t> OverlappingOverlapSplicer::seek_to_path_base(
        MutablePathDeletableHandleGraph& gfa_graph,
        OverlappingChild& overlapping_child,
        size_t target_base_index){

    PathInfo path_info;
    string path_name;
    find_path_info(gfa_graph, overlapping_child.biclique_index, overlapping_child.handle, path_info, path_name);

    auto path_handle = gfa_graph.get_path_handle(path_name);
    auto step = gfa_graph.path_begin(path_handle);

    size_t cumulative_index = 0;
    size_t intra_handle_index = 0;
    handle_t step_handle;

    bool fail = true;
    while (step != gfa_graph.path_end(path_handle)){
        step_handle = gfa_graph.get_handle_of_step(step);
        auto step_length = gfa_graph.get_length(step_handle);

        if (cumulative_index + step_length > target_base_index){
            intra_handle_index = target_base_index - cumulative_index;
            fail = false;
            break;
        }
        else{
            gfa_graph.get_next_step(step);
        }

        cumulative_index += step_length;
    }

    if (fail){
        throw runtime_error("ERROR: path base index " + to_string(target_base_index) + " exceeds sum of handle lengths");
    }

    return {step_handle, intra_handle_index};
}


void OverlappingOverlapSplicer::splice_overlapping_overlaps(MutablePathDeletableHandleGraph& gfa_graph) {

    for (auto& item: overlapping_overlap_nodes) {
        auto node_id = item.first;
        auto parent_handle = gfa_graph.get_handle(node_id, false);
        auto& overlap_info = item.second;
        auto parent_size = overlap_info.length;

        cout << "Splicing overlapping overlap: " << node_id << '\n';

        overlap_info.print(gfa_graph);

        for (size_t i=0; i<overlap_info.length; i++){

//            seek_to_path_base(gfa_graph, )
        }
    }
}



}
