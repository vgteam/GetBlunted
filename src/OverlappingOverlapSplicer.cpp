#include "OverlappingOverlapSplicer.hpp"
#include "handle_to_gfa.hpp"
#include "handlegraph/util.hpp"


using std::cerr;
using std::tie;

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


bool OverlappingOverlapSplicer::find_path_info(
        const HandleGraph& gfa_graph,
        size_t biclique_index,
        handle_t handle,
        PathInfo& path_info,
        string& path_name){

    auto& subgraph = subgraphs[biclique_index];

    bool reversal = false;

    // Don't know which side of the biclique this overlap was on until we search for it in the subgraph
    auto result = subgraph.paths_per_handle[0].find(handle);

    if (result == subgraph.paths_per_handle[0].end()) {
        result = subgraph.paths_per_handle[1].find(handle);

        if (result == subgraph.paths_per_handle[1].end()) {
            result = subgraph.paths_per_handle[0].find(gfa_graph.flip(handle));

            if (result == subgraph.paths_per_handle[0].end()) {
                result = subgraph.paths_per_handle[1].find(gfa_graph.flip(handle));

                // Sanity check
                if (result == subgraph.paths_per_handle[1].end()) {
                    throw runtime_error("ERROR: node not found in biclique subgraph. Node id: " +
                                        to_string(gfa_graph.get_id(handle)));
                }
                else{
                    reversal = true;
                }
            }
            else{
                reversal = true;
            }
        }
    }

    path_name = subgraph.graph.get_path_name(result->second.path_handle);
    path_info = result->second;

    return reversal;
}


tuple<handle_t, size_t, size_t, bool> OverlappingOverlapSplicer::seek_to_path_base(
        MutablePathDeletableHandleGraph& gfa_graph,
        string& path_name,
        size_t target_base_index){

    cerr << "seeking path: " << path_name << '\n';

    auto path_handle = gfa_graph.get_path_handle(path_name);
    auto step = gfa_graph.path_begin(path_handle);

    size_t cumulative_index = 0;
    size_t intra_handle_index = 0;
    size_t remainder = 0;
    handle_t step_handle;

    bool fail = true;
    while (step != gfa_graph.path_end(path_handle)){
        step_handle = gfa_graph.get_handle_of_step(step);
        auto step_length = gfa_graph.get_length(step_handle);

        cerr << target_base_index << " " << cumulative_index << " " << intra_handle_index << " " << gfa_graph.get_sequence(step_handle) << '\n';

        if (cumulative_index + step_length > target_base_index){
            intra_handle_index = target_base_index - cumulative_index;
            fail = false;
            break;
        }
        else{
            step = gfa_graph.get_next_step(step);
        }

        cumulative_index += step_length;
    }

    remainder = target_base_index - cumulative_index;

    cerr << (fail ? "FAIL with remainder " + to_string(remainder) : "PASS with remainder " + to_string(remainder)) << '\n';

    return {step_handle, intra_handle_index, remainder, fail};
}




tuple<handle_t, size_t, size_t, bool> OverlappingOverlapSplicer::seek_to_reverse_path_base(
        MutablePathDeletableHandleGraph& gfa_graph,
        string& path_name,
        size_t target_base_index){

    cerr << "seeking REVERSE path: " << path_name << '\n';

    auto path_handle = gfa_graph.get_path_handle(path_name);
    auto step = gfa_graph.path_back(path_handle);

    size_t cumulative_index = 0;
    size_t intra_handle_index = 0;
    size_t remainder = 0;
    handle_t step_handle;

    bool fail = true;
    while (step != gfa_graph.path_front_end(path_handle)){
        step_handle = gfa_graph.flip(gfa_graph.get_handle_of_step(step));
        auto step_length = gfa_graph.get_length(step_handle);

        cerr << target_base_index << " " << cumulative_index << " " << intra_handle_index << " " << gfa_graph.get_sequence(step_handle) << '\n';

        if (cumulative_index + step_length > target_base_index){
            intra_handle_index = target_base_index - cumulative_index;
            fail = false;
            break;
        }
        else{
            step = gfa_graph.get_previous_step(step);
        }

        cumulative_index += step_length;
    }

    remainder = target_base_index - cumulative_index;

    cerr << (fail ? "FAIL with remainder " + to_string(remainder) : "PASS with remainder " + to_string(remainder)) << '\n';

    return {step_handle, intra_handle_index, remainder, fail};
}



void OverlappingOverlapSplicer::find_splice_pairs(
        HandleGraph& gfa_graph,
        OverlappingNodeInfo& overlap_info,
        vector <OverlappingSplicePair>& oo_splice_pairs){

    array <vector <overlapping_child_iter>, 2> other_children;

    overlap_info.print(gfa_graph);

    for (auto side: {0,1}) {
        for (auto oo: overlap_info.overlapping_children[side]){

//            if (oo.first == 0 or oo.first == overlap_info.length-1){
//                throw runtime_error("ERROR: Overlapping overlap meets or exceeds node length: "
//                                     + to_string(overlap_info.parent_node));
//            }

            size_t n_normal_children_exceeding_overlap;
            bool splice_to_parent = false;

            if (1-side == 0) {
                greater_than_or_equal(overlap_info.normal_children[1-side], oo.first - 1, other_children[1-side]);
                n_normal_children_exceeding_overlap = other_children[1-side].size();

                greater_than_or_equal(overlap_info.overlapping_children[1-side], oo.first - 1, other_children[1-side]);
            }
            else{
                less_than_or_equal(overlap_info.normal_children[1-side], oo.first + 1, other_children[1-side]);
                n_normal_children_exceeding_overlap = other_children[1-side].size();

                less_than_or_equal(overlap_info.overlapping_children[1-side], oo.first + 1, other_children[1-side]);
            }

            // Check if there are some overlaps that would be unreachable by this OO
            if (n_normal_children_exceeding_overlap < overlap_info.normal_children[1-side].size()){
                splice_to_parent = true;
            }

            if (not other_children[1-side].empty()) {
                for (auto& other: other_children[1 - side]) {
                    auto& oo_child = oo.second;
                    auto& other_child = other->second;

                    OverlappingChild left_child;
                    OverlappingChild right_child;

                    OverlappingSplicePair splice_pair_a;
                    OverlappingSplicePair splice_pair_b;
                    splice_pair_a.side = 0;
                    splice_pair_b.side = 1;

                    bool left_reversal;
                    bool right_reversal;

                    if (side == 0) {
                        //
                        //  oo     [0]-[1]-[2]-[3]-[4]
                        //                \ a         \ b
                        //  other          [0]-[1]-[2]-[3]-[4]
                        //
                        //  parent [0] [1] [2] [3] [4] [5] [6]
                        //

                        PathInfo path_info;
                        string left_child_path_name;
                        string right_child_path_name;

                        left_child = oo_child;
                        right_child = other_child;

                        left_reversal = find_path_info(
                                gfa_graph,
                                left_child.biclique_index,
                                left_child.handle,
                                path_info,
                                left_child_path_name);

                        right_reversal = find_path_info(
                                gfa_graph,
                                right_child.biclique_index,
                                right_child.handle,
                                path_info,
                                right_child_path_name);

                        splice_pair_a.left_reversal = left_reversal;
                        splice_pair_b.left_reversal = left_reversal;

                        splice_pair_a.right_reversal = right_reversal;
                        splice_pair_b.right_reversal = right_reversal;

                        splice_pair_a.left_child_path_name = left_child_path_name;
                        splice_pair_b.left_child_path_name = left_child_path_name;

                        splice_pair_a.right_child_path_name = right_child_path_name;
                        splice_pair_b.right_child_path_name = right_child_path_name;

                        splice_pair_a.left_parent_index = other->first - 1;                             // 1
                        splice_pair_b.left_parent_index = oo.first;                                     // 4

                        splice_pair_a.right_parent_index = splice_pair_a.left_parent_index + 1;         // 2
                        splice_pair_b.right_parent_index = splice_pair_b.left_parent_index + 1;         // 5

                        splice_pair_a.left_child_index = splice_pair_a.left_parent_index;               // 1
                        splice_pair_b.left_child_index = splice_pair_b.left_parent_index;               // 4

                        splice_pair_a.right_child_index = 0;                                            // 0
                        splice_pair_b.right_child_index = oo.first - splice_pair_a.left_parent_index;   // 4 - 1 = 3

                    }
                    else{
                        //
                        //  other  [0]-[1]-[2]-[3]-[4]
                        //                \ a         \ b
                        //  oo             [0]-[1]-[2]-[3]-[4]
                        //
                        //  parent [0] [1] [2] [3] [4] [5] [6]
                        //

                        //
                        //  other  [0]-[1]-[2]-[3]
                        //            \ a         \ b
                        //  oo         [0]-[1]-[2]-[3]
                        //
                        //  parent [0] [1] [2] [3] [4]
                        //


                        PathInfo path_info;
                        string left_child_path_name;
                        string right_child_path_name;

                        left_child = other_child;
                        right_child = oo_child;

                        left_reversal = find_path_info(
                                gfa_graph,
                                left_child.biclique_index,
                                left_child.handle,
                                path_info,
                                left_child_path_name);

                        right_reversal = find_path_info(
                                gfa_graph,
                                right_child.biclique_index,
                                right_child.handle,
                                path_info,
                                right_child_path_name);

                        splice_pair_a.left_reversal = left_reversal;
                        splice_pair_b.left_reversal = left_reversal;

                        splice_pair_a.right_reversal = right_reversal;
                        splice_pair_b.right_reversal = right_reversal;

                        splice_pair_a.left_child_path_name = left_child_path_name;
                        splice_pair_b.left_child_path_name = left_child_path_name;

                        splice_pair_a.right_child_path_name = right_child_path_name;
                        splice_pair_b.right_child_path_name = right_child_path_name;

                        splice_pair_a.left_parent_index = oo.first - 1;                                    // 1
                        splice_pair_b.left_parent_index = other->first;                                    // 4

                        splice_pair_a.right_parent_index = splice_pair_a.left_parent_index + 1;            // 2  1
                        splice_pair_b.right_parent_index = splice_pair_b.left_parent_index + 1;            // 5  4

                        splice_pair_a.left_child_index = splice_pair_a.left_parent_index;                  // 1  0
                        splice_pair_b.left_child_index = splice_pair_b.left_parent_index;                  // 4  3

                        splice_pair_a.right_child_index = 0;
                        splice_pair_b.right_child_index = other->first - splice_pair_a.left_parent_index;  // 4 - 1 = 3  3 - 0 = 3
                    }

                    oo_splice_pairs.emplace_back(splice_pair_a);
                    oo_splice_pairs.emplace_back(splice_pair_b);
                }
            }

            // If there are some overlaps on the other side of the OO that are not reached by the OO,
            // the OO must be spliced into the parent node so they are reachable
            if (splice_to_parent){
                auto& oo_child = oo.second;
                OverlappingSplicePair splice_pair;

                bool left_reversal;
                bool right_reversal;

                if (side == 0) {
                    //
                    //  oo     [0]-[1]
                    //                \
                    //  parent [0]-[1]-[2]-[3]
                    //

                    PathInfo path_info;
                    string left_child_path_name;

                    left_reversal = find_path_info(
                            gfa_graph,
                            oo_child.biclique_index,
                            oo_child.handle,
                            path_info,
                            left_child_path_name);

                    right_reversal = false;

                    splice_pair.left_reversal = left_reversal;
                    splice_pair.right_reversal = right_reversal;

                    splice_pair.left_child_path_name = left_child_path_name;
                    splice_pair.right_child_path_name = overlap_info.parent_path_name;

                    splice_pair.left_parent_index = oo.first;
                    splice_pair.right_parent_index = splice_pair.left_parent_index + 1;

                    splice_pair.left_child_index = splice_pair.left_parent_index;
                    splice_pair.right_child_index = splice_pair.right_parent_index;

                }
                else {
                    //
                    //  oo             [0]-[1]
                    //                /
                    //  parent [0]-[1]-[2]-[3]
                    //

                    PathInfo path_info;
                    string right_child_path_name;

                    right_reversal = find_path_info(
                            gfa_graph,
                            oo_child.biclique_index,
                            oo_child.handle,
                            path_info,
                            right_child_path_name);

                    left_reversal = false;

                    splice_pair.left_reversal = left_reversal;
                    splice_pair.right_reversal = right_reversal;

                    splice_pair.left_child_path_name = overlap_info.parent_path_name;
                    splice_pair.right_child_path_name = right_child_path_name;

                    splice_pair.left_parent_index = oo.first - 1;
                    splice_pair.right_parent_index = splice_pair.left_parent_index + 1;

                    splice_pair.left_child_index = splice_pair.left_parent_index;
                    splice_pair.right_child_index = 0;
                }

                oo_splice_pairs.emplace_back(splice_pair);
            }
        }
    }
}


void OverlappingOverlapSplicer::splice_overlapping_overlaps(MutablePathDeletableHandleGraph& gfa_graph) {

    for (auto& oo_item: overlapping_overlap_nodes) {
        // Output an image of the graph, can be uncommented for debugging
        {
            string test_path_prefix = "test_bluntify_splice_" + to_string(oo_item.first);
            std::ofstream test_output(test_path_prefix + ".gfa");
            handle_graph_to_gfa(gfa_graph, test_output);
            test_output.close();

            if (gfa_graph.get_node_count() < 200) {
                string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                                 + test_path_prefix + ".png";

                cerr << "Running: " << command << '\n';

                run_command(command);
            }
        }

        auto node_id = oo_item.first;
        auto& overlap_info = oo_item.second;

        map<handle_t, set<size_t> > division_sites;
        vector<OverlappingSplicePair> oo_splice_pairs;

        // Find all division points
        find_splice_pairs(gfa_graph, overlap_info, oo_splice_pairs);

        // Aggregate the splice pairs by the handles that they splice
        for (auto& splice_pair: oo_splice_pairs) {
            handle_t left_handle;
            handle_t right_handle;
            size_t left_index;
            size_t right_index;
            size_t left_remainder;
            size_t right_remainder;
            bool left_fail;
            bool right_fail;

            if (splice_pair.left_reversal) {
                tie(left_handle, left_index, left_remainder, left_fail) = seek_to_reverse_path_base(
                        gfa_graph,
                        splice_pair.left_child_path_name,
                        splice_pair.left_child_index);
            }
            else{
                tie(left_handle, left_index, left_remainder, left_fail) = seek_to_path_base(
                        gfa_graph,
                        splice_pair.left_child_path_name,
                        splice_pair.left_child_index);
            }

            if (splice_pair.right_reversal) {
                tie(right_handle, right_index, right_remainder, right_fail) = seek_to_reverse_path_base(
                        gfa_graph,
                        splice_pair.right_child_path_name,
                        splice_pair.right_child_index);
            }
            else{
                tie(right_handle, right_index, right_remainder, right_fail) = seek_to_path_base(
                        gfa_graph,
                        splice_pair.right_child_path_name,
                        splice_pair.right_child_index);
            }

            if (left_index + 1 < gfa_graph.get_length(left_handle) and not left_fail) {
                division_sites[left_handle].emplace(left_index + 1);
                cerr << "Splice site for parent node " << node_id << " with children " << gfa_graph.get_id(left_handle) << " (" << gfa_graph.get_id(right_handle) << ") = " << left_index + 1 << '\n';
            }

            if (right_index < gfa_graph.get_length(right_handle) and right_index > 0 and not right_fail) {
                division_sites[right_handle].emplace(right_index);
                cerr << "Splice site for parent node " << node_id << " with children " << gfa_graph.get_id(right_handle) << " (" << gfa_graph.get_id(left_handle) << ") = " << right_index << '\n';
            }
        }

        cerr << "Dividing nodes at oo splice sites... " << '\n';

        // Do the divisions in bulk
        for (auto& item: division_sites) {
            vector<size_t> sites;
            for (auto& i: item.second){
                if (i > 0) {
                    sites.emplace_back(i);
                }
            }

            gfa_graph.divide_handle(item.first, sites);
        }

        // Splice each division point
        // (doesn't really need to use seek_path but also this entire project doesn't really need to exist so...)
        for (auto& splice_pair: oo_splice_pairs) {
            handle_t left_handle;
            handle_t right_handle;
            size_t left_index;
            size_t right_index;
            size_t left_remainder;
            size_t right_remainder;
            size_t left_fail;
            size_t right_fail;

            if (splice_pair.left_reversal) {
                tie(left_handle, left_index, left_remainder, left_fail) = seek_to_reverse_path_base(
                        gfa_graph,
                        splice_pair.left_child_path_name,
                        splice_pair.left_child_index);
            } else {
                tie(left_handle, left_index, left_remainder, left_fail) = seek_to_path_base(
                        gfa_graph,
                        splice_pair.left_child_path_name,
                        splice_pair.left_child_index);
            }

            if (splice_pair.right_reversal) {
                tie(right_handle, right_index, right_remainder, right_fail) = seek_to_reverse_path_base(
                        gfa_graph,
                        splice_pair.right_child_path_name,
                        splice_pair.right_child_index);
            } else {
                tie(right_handle, right_index, right_remainder, right_fail) = seek_to_path_base(
                        gfa_graph,
                        splice_pair.right_child_path_name,
                        splice_pair.right_child_index);
            }


            bool left_side_is_end = (left_fail and left_remainder == 0) and (not right_fail and right_remainder > 0);
            bool right_side_is_end = (not left_fail and left_remainder > 0) and (right_fail and right_remainder == 0);
            bool both_sides_are_end = (left_remainder == 0 and right_remainder == 0);

            //
            if (not left_fail and not right_fail) {
                gfa_graph.create_edge(left_handle, right_handle);
            }
            // In the following case, the overlapping overlap meets the end of the node, which means its edges
            // are going to need to be spliced to the neighboring nodes
            else if (both_sides_are_end or left_side_is_end or right_side_is_end){
                vector <handle_t> left_splice_handles;
                vector <handle_t> right_splice_handles;

                if (splice_pair.side == 0) {
                    gfa_graph.follow_edges(left_handle, true, [&](handle_t h) {
                        left_splice_handles.emplace_back(h);
                    });

                    right_splice_handles.emplace_back(right_handle);
                }
                else {
                    gfa_graph.follow_edges(right_handle, false, [&](handle_t h) {
                        right_splice_handles.emplace_back(h);
                    });

                    left_splice_handles.emplace_back(left_handle);
                }

                for (auto& l: left_splice_handles){
                    for (auto& r: right_splice_handles){
                        gfa_graph.create_edge(l, r);
                    }
                }
            }
            else{
                throw runtime_error("ERROR: overlap length is > parent node/path length by "
                                    + to_string(std::max(right_remainder, left_remainder)));
            }
        }
    }
}



}