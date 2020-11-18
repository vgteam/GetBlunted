#include "OverlappingOverlapSplicer.hpp"
#include "handlegraph/util.hpp"


using std::cout;

namespace bluntifier {

typedef map <size_t, OverlappingChild>::iterator overlapping_child_iter;


OverlappingOverlapSplicer::OverlappingOverlapSplicer(
        map<nid_t, OverlappingNodeInfo>& overlapping_overlap_nodes,
        map <nid_t, set<nid_t> >& parent_to_children,
        const vector<Subgraph>& subgraphs):
    overlapping_overlap_nodes(overlapping_overlap_nodes),
    parent_to_children(parent_to_children),
    subgraphs(subgraphs)
{}


void OverlappingOverlapSplicer::find_parent_path_bounds(
        MutablePathDeletableHandleGraph& gfa_graph,
        nid_t parent_id,
        pair<size_t, size_t>& bounds) {

    size_t i_step = 0;

    auto parent_path_handle = gfa_graph.get_path_handle(to_string(parent_id));
    for (auto h: gfa_graph.scan_path(parent_path_handle)) {
        auto id = gfa_graph.get_id(h);
        std::cout << id << " " << gfa_graph.get_sequence(h);

        // If this node is not one of the suffixes/prefixes
        if (parent_to_children.at(parent_id).count(id) == 0) {

        }
        // Otherwise remember the amount of trimmed sequence
        else {
            std::cout << " <- terminus";
            if (i_step == 0) {
                bounds.first = gfa_graph.get_length(h);
            } else {
            // If it's not the first step and this handle is a prefix/suffix then the path ends
                bounds.second = gfa_graph.get_length(h);
                break;
            }
        }
        std::cout << '\n';

        i_step++;
    }
    std::cout << '\n';
    std::cout << '\n';

    gfa_graph.destroy_path(parent_path_handle);

//    if (sorted_extents_per_side[0].empty()) {
//        if (left_trimmed != 0) {
//            throw runtime_error("ERROR: parent path trim size does not equal smallest overlap size: " +
//                                to_string(node_info.node_id) + " on left side\n" +
//                                "\tTrimmed: " + to_string(left_trimmed) + "\n"
//                                                                          "\tNo overlap: " +
//                                to_string(sorted_extents_per_side[0].empty()));
//        }
//    }
//    else if (left_trimmed != sorted_extents_per_side[0].back()){
//        throw runtime_error("ERROR: parent path trim size does not equal smallest overlap size: " +
//                            to_string(node_info.node_id) + " on left side\n" +
//                            "\tTrimmed: " + to_string(left_trimmed) + "\n"
//                            "\tMin overlap: " + to_string(sorted_extents_per_side[0].back()));
//    }
//
//    if (sorted_extents_per_side[1].empty()) {
//        if (right_trimmed != 0) {
//            throw runtime_error("ERROR: parent path trim size does not equal smallest overlap size: " +
//                                to_string(node_info.node_id) + " on right side\n" +
//                                "\tTrimmed: " + to_string(right_trimmed) + "\n"
//                                                                           "\tNo overlap: " +
//                                to_string(sorted_extents_per_side[1].empty()));
//        }
//    }
//    else if (right_trimmed != sorted_extents_per_side[1].back()){
//        throw runtime_error("ERROR: parent path trim size does not equal smallest overlap size: " +
//                            to_string(node_info.node_id) + " on right side\n" +
//                            "\tTrimmed: " + to_string(right_trimmed) + "\n"
//                            "\tMin overlap: " + to_string(sorted_extents_per_side[1].back()));
//    }
}

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

    cout << "Seeking path base: " << target_base_index << '\n';

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
            cout << gfa_graph.get_id(step_handle) << " " << cumulative_index << " " << target_base_index << " " << gfa_graph.get_sequence(step_handle) << '\n';
            step = gfa_graph.get_next_step(step);
        }

        cumulative_index += step_length;
    }
    cout << gfa_graph.get_id(step_handle) << " " << gfa_graph.get_sequence(step_handle) << '\n';
    cout << '\n';

    if (fail){
        throw runtime_error("ERROR: path base index " + to_string(target_base_index) + " exceeds sum of handle lengths");
    }

    return {step_handle, intra_handle_index};
}


void OverlappingOverlapSplicer::find_splice_pairs(
        HandleGraph& gfa_graph,
        OverlappingNodeInfo& overlap_info,
        vector <OverlappingSplicePair>& oo_splice_pairs){

    array <vector <overlapping_child_iter>, 2> other_children;
    set <pair <handle_t,handle_t> > visited;

    for (auto side: {0,1}) {
        for (auto oo: overlap_info.overlapping_children[side]){

            if (oo.first == 0 or oo.first == overlap_info.length-1){
                throw runtime_error("ERROR: Overlapping overlap meets or exceeds node length");
            }

            if (1-side == 0) {
                greater_than_or_equal(overlap_info.overlapping_children[1-side], oo.first - 1, other_children[1-side]);
                greater_than_or_equal(overlap_info.normal_children[1-side], oo.first - 1, other_children[1-side]);
            }
            else{
                less_than_or_equal(overlap_info.overlapping_children[1-side], oo.first + 1, other_children[1-side]);
                less_than_or_equal(overlap_info.normal_children[1-side], oo.first + 1, other_children[1-side]);
            }

            cout << "Splice site: " << oo.first << '\n';

            if (not other_children[1-side].empty()) {
                for (auto& normie: other_children[1 - side]) {
                    cout << "\tNormal  - side: " << 1 - side << " - " << "Index=" << normie->first << " ";
                    normie->second.print(gfa_graph);
                }
                cout << "\tOverlap - side: " << side << " - " << "Index=" << oo.first << " ";
                oo.second.print(gfa_graph);

                for (auto& other: other_children[1 - side]) {
                    auto& oo_child = oo.second;
                    auto& other_child = other->second;

                    if (oo_child.handle == other_child.handle){
                        continue;
                    }

                    OverlappingSplicePair splice_pair_a;
                    OverlappingSplicePair splice_pair_b;

                    if (side == 0) {
                        //
                        //  oo     [0]-[1]-[2]-[3]-[4]
                        //                \ a         \ b
                        //  other          [0]-[1]-[2]-[3]-[4]
                        //
                        //
                        //  parent [0] [1] [2] [3] [4] [5] [6]
                        //

                        splice_pair_a.left_child = oo_child;
                        splice_pair_b.left_child = oo_child;

                        splice_pair_a.right_child = other_child;
                        splice_pair_b.right_child = other_child;

                        splice_pair_a.left_parent_index = other->first - 1;
                        splice_pair_b.left_parent_index = oo.first;

                        splice_pair_a.right_parent_index = splice_pair_a.left_parent_index + 1;
                        splice_pair_b.right_parent_index = splice_pair_b.left_parent_index + 1;

                        splice_pair_a.left_child_index = splice_pair_a.left_parent_index;
                        splice_pair_b.left_child_index = splice_pair_b.left_parent_index;

                        splice_pair_a.right_child_index = 0;
                        splice_pair_b.right_child_index = oo.first - splice_pair_a.left_child_index;

                        splice_pair_a.left_length = gfa_graph.get_length(splice_pair_a.left_child.handle);
                        splice_pair_b.left_length = gfa_graph.get_length(splice_pair_b.left_child.handle);

                        splice_pair_a.right_length = gfa_graph.get_length(splice_pair_a.right_child.handle);
                        splice_pair_b.right_length = gfa_graph.get_length(splice_pair_b.right_child.handle);
                    }
                    else{
                        //
                        //  other  [0]-[1]-[2]-[3]-[4]
                        //                \ a         \ b
                        //  oo             [0]-[1]-[2]-[3]-[4]
                        //
                        //  parent [0] [1] [2] [3] [4] [5] [6]
                        //

                        splice_pair_a.left_child = other_child;
                        splice_pair_b.left_child = other_child;

                        splice_pair_a.right_child = oo_child;
                        splice_pair_b.right_child = oo_child;

                        splice_pair_a.left_parent_index = oo.first - 1;
                        splice_pair_b.left_parent_index = other->first;

                        splice_pair_a.right_parent_index = splice_pair_a.left_parent_index + 1;
                        splice_pair_b.right_parent_index = splice_pair_b.left_parent_index + 1;

                        splice_pair_a.left_child_index = splice_pair_a.left_parent_index;
                        splice_pair_b.left_child_index = splice_pair_a.left_parent_index;

                        splice_pair_a.right_child_index = 0;
                        splice_pair_b.right_child_index = other->first - splice_pair_a.left_child_index;

                        splice_pair_a.left_length = gfa_graph.get_length(splice_pair_a.left_child.handle);
                        splice_pair_b.left_length = gfa_graph.get_length(splice_pair_b.left_child.handle);

                        splice_pair_a.right_length = gfa_graph.get_length(splice_pair_a.right_child.handle);
                        splice_pair_b.right_length = gfa_graph.get_length(splice_pair_b.right_child.handle);
                    }

                    oo_splice_pairs.emplace_back(splice_pair_a);
                    oo_splice_pairs.emplace_back(splice_pair_b);
                    cout << "";
                }
            }
        }
    }
}


void OverlappingOverlapSplicer::splice_overlapping_overlaps(MutablePathDeletableHandleGraph& gfa_graph) {

    for (auto& oo_item: overlapping_overlap_nodes) {
        auto node_id = oo_item.first;
        auto parent_handle = gfa_graph.get_handle(node_id, false);
        auto& overlap_info = oo_item.second;
        auto parent_size = overlap_info.length;

        cout << "Splicing overlapping overlap: " << node_id << '\n';

        overlap_info.print(gfa_graph);

        map<handle_t, set<size_t> > division_sites;
        vector<OverlappingSplicePair> oo_splice_pairs;

        // Find all division points
        find_splice_pairs(gfa_graph, overlap_info, oo_splice_pairs);

        // Break the nodes if there is a splice into their midsection
        for (auto& splice_pair: oo_splice_pairs) {
            if (splice_pair.left_child_index < splice_pair.left_length and splice_pair.left_child_index > 1) {
                auto[left_handle, left_index] = seek_to_path_base(gfa_graph, splice_pair.left_child, splice_pair.left_child_index + 1);
                division_sites[left_handle].emplace(left_index);
            }
            if (splice_pair.right_child_index < splice_pair.right_length and splice_pair.right_child_index > 1) {
                auto[right_handle, right_index] = seek_to_path_base(gfa_graph, splice_pair.right_child, splice_pair.right_child_index);
                division_sites[right_handle].emplace(right_index);
            }
        }


        for (auto& item: division_sites) {
            vector<size_t> sites;
            for (auto& i: item.second){
                if (i > 0) {
                    sites.emplace_back(i);
                }
            }

            gfa_graph.divide_handle(item.first, sites);
        }

        for (auto& splice_pair: oo_splice_pairs) {
            auto[left_handle, left_index] = seek_to_path_base(gfa_graph, splice_pair.left_child, splice_pair.left_child_index);
            auto[right_handle, right_index] = seek_to_path_base(gfa_graph, splice_pair.right_child, splice_pair.right_child_index);

            cout << "Creating (" << gfa_graph.get_id(left_handle);
            cout << (gfa_graph.get_is_reverse(left_handle) ? "-" : "+");
            cout << ") -> (" << gfa_graph.get_id(right_handle);
            cout << (gfa_graph.get_is_reverse(right_handle) ? "-" : "+") << ")" << '\n';
            gfa_graph.create_edge(left_handle, right_handle);
        }
    }
}



}
