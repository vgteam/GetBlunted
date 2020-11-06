#include "OverlappingOverlapSplicer.hpp"


using std::cout;

namespace bluntifier {


void find_path_info(
        const vector<Subgraph>& subgraphs,
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


void splice_overlapping_overlaps(
        MutablePathDeletableHandleGraph& gfa_graph,
        vector <Subgraph>& subgraphs,
        map<nid_t, OverlappingNodeInfo>& overlapping_overlap_nodes) {

    for (auto& a: overlapping_overlap_nodes) {
        auto node_id = a.first;
        auto parent_handle = gfa_graph.get_handle(node_id, false);
        auto& overlap_info = a.second;
        auto parent_size = overlap_info.length;

        cout << "Splicing overlapping overlap: " << node_id << '\n';

        for (auto& b: overlap_info.overlapping_children){
            auto& overlapping_child_info = b.second;
            auto overlapping_child_size = gfa_graph.get_length(overlapping_child_info.handle);

            PathInfo overlapping_path_info;
            string overlapping_path_name;

            find_path_info(
                    subgraphs,
                    gfa_graph,
                    overlapping_child_info.biclique_index,
                    overlapping_child_info.handle,
                    overlapping_path_info,
                    overlapping_path_name);

            overlapping_child_info.print(gfa_graph);

            cout << "This side\n";
            for (auto& item: overlap_info.biclique_side_to_child[overlapping_child_info.side]) {
                auto other_biclique = item.first;
                auto other_handle = item.second;
                auto other_size = gfa_graph.get_length(other_handle);

                cout << '\t' << gfa_graph.get_id(other_handle) << " Biclique=" << other_biclique << " Length="
                     << other_size << '\n';
            }

            cout << "Other side\n";
            // Start looking for overlaps on the other side of the NODE (not biclique).
            // Perform splicing for each one that overlaps this overlap
            bool other_side = 1 - overlapping_child_info.side;
            for (auto& item: overlap_info.biclique_side_to_child[other_side]){
                auto other_biclique = item.first;
                auto other_handle = item.second;
                auto other_size = gfa_graph.get_length(other_handle);

                cout << '\t' << gfa_graph.get_id(other_handle) << " Biclique=" << other_biclique << " Length=" << other_size << '\n';

                // Dont splice this node if it's actually not overlapping
                if (other_size + overlapping_child_size < parent_size){
                    continue;
                }

                PathInfo other_path_info;
                string other_path_name;
                find_path_info(
                        subgraphs,
                        gfa_graph,
                        other_biclique,
                        other_handle,
                        other_path_info,
                        other_path_name);

                auto other_path_handle = gfa_graph.get_path_handle(other_path_name);
                auto step = gfa_graph.path_back(other_path_handle);
                size_t traversed_bases = 0;

                while (step != gfa_graph.path_front_end(other_path_handle)){
                    auto h = gfa_graph.get_handle_of_step(step);
                    auto s = gfa_graph.get_length(h);

                    cout << "Traversing node in subgraph: " << gfa_graph.get_id(h) << '\n';

                    size_t i = s - 1;
                    while (i >= 0 and other_size - traversed_bases - i + overlapping_child_size > parent_size){
//                        gfa_graph.divide_handle(h);
                        cout << gfa_graph.get_base(h,i) << '\n';
                        i++;
                    }

                    step = gfa_graph.get_previous_step(step);
                    traversed_bases += i;
                }
            }

            // Splice into the remaining body of the node if any of the other overlaps don't overlap this one
        }
    }
}



}
