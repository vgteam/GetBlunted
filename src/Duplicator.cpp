#include "Duplicator.hpp"

namespace bluntifier{


OverlappingNodeInfo::OverlappingNodeInfo(nid_t parent_node):
    parent_node(parent_node)
{}


Duplicator::Duplicator(
        const vector <vector <BicliqueEdgeIndex> >& node_to_biclique_edge,
        Bicliques& bicliques,
        OverlapMap& overlaps):
        node_to_biclique_edge(node_to_biclique_edge),
        bicliques(bicliques),
        overlaps(overlaps)
{}


void Duplicator::duplicate_overlapping_termini(
        MutablePathDeletableHandleGraph& gfa_graph,
        const array <deque <size_t>, 2>& sorted_sizes_per_side,
        const array <deque <size_t>, 2>& sorted_bicliques_per_side,
        array<map<size_t, handle_t>, 2>& biclique_side_to_child,
        const NodeInfo& node_info){

    auto result = overlapping_overlap_nodes.emplace(node_info.node_id, node_info.node_id);
    auto& overlap_node_info = result.first->second;

    auto parent_handle = gfa_graph.get_handle(node_info.node_id, false);

    // Duplicate all the LEFT side biclique overlaps (using only the maximum length overlap)
    for (size_t i=0; i<sorted_sizes_per_side[0].size(); i++){
        size_t s = sorted_sizes_per_side[0][i];
        string sequence = gfa_graph.get_subsequence(parent_handle, 0, s);

        auto child = gfa_graph.create_handle(sequence);

        size_t biclique_index = sorted_bicliques_per_side[0][i];
        biclique_side_to_child[0].emplace(biclique_index, child);

        // Update provenance map
        auto child_node = gfa_graph.get_id(child);
        child_to_parent[child_node] = node_info.node_id;
        parent_to_children[node_info.node_id].emplace(child_node);
    }

    // Duplicate all the RIGHT side biclique overlaps (using only the maximum length overlap)
    for (size_t i=0; i<sorted_sizes_per_side[1].size(); i++){
        size_t start = gfa_graph.get_length(parent_handle) - sorted_sizes_per_side[1][i];
        string sequence = gfa_graph.get_subsequence(parent_handle, start, sorted_sizes_per_side[1][i]);

        auto child = gfa_graph.create_handle(sequence);

        size_t biclique_index = sorted_bicliques_per_side[1][i];
        biclique_side_to_child[1].emplace(biclique_index, child);

        // Update provenance map
        auto child_node = gfa_graph.get_id(child);
        child_to_parent[child_node] = node_info.node_id;
        parent_to_children[node_info.node_id].emplace(child_node);
    }

    // Copy the mapping from side -> biclique -> child node into the overlapping_node_info so that it can be used later
    // because each child is an orphaned from the others until after POA and re-splicing
    // This could be done without copying but I just really dgaf anymore
    overlap_node_info.biclique_side_to_child = biclique_side_to_child;
    overlap_node_info.length = gfa_graph.get_length(parent_handle);

    // Die a horrible death you vile monster
    gfa_graph.destroy_handle(parent_handle);
}


void Duplicator::repair_edges(
        MutablePathDeletableHandleGraph& gfa_graph,
        const array <map <size_t, handle_t>, 2>& biclique_side_to_child,
        handle_t old_handle,
        handle_t old_handle_flipped) {

//    cout << "Children generated per biclique:" << '\n';
//    for (const size_t side: {0, 1}) {
//        for (const auto&[biclique_index, child_handle]: biclique_side_to_child[side]) {
//            cout << side << " " << biclique_index << " " << gfa_graph.get_id(child_handle) << " ";
//
//            if (gfa_graph.get_length(child_handle) < 60){
//                cout << gfa_graph.get_sequence(child_handle);
//            }
//
//            cout << '\n';
//        }
//    }


    for (const size_t side: {0, 1}) {
        for (const auto&[biclique_index, child_handle]: biclique_side_to_child[side]) {
            for (auto& edge: bicliques[biclique_index]) {
                auto old_edge = edge;

                bool edge_found = false;
                if (side == 0) {
                    if (edge.second == old_handle) {
                        // If this is a loop, look for the corresponding biclique child on the other side of the node
                        if (edge.first == edge.second){
//                            cout << "Loop ";
                            edge.first = biclique_side_to_child[1-side].at(biclique_index);
                        }

                        edge.second = child_handle;

                        edge_found = true;

//                        cout << "Case A\n";
                    } else if (edge.first == old_handle_flipped) {
                        // If this is a loop, look for the corresponding biclique child on the other side of the node
                        if (edge.first == edge.second){
//                            cout << "Loop ";
                            edge.second = biclique_side_to_child[1-side].at(biclique_index);
                        }

                        edge.first = gfa_graph.flip(child_handle);

                        edge_found = true;

//                        cout << "Case B\n";
                    }
                }
                else {
                    if (edge.first == old_handle) {
                        // If this is a loop, look for the corresponding biclique child on the other side of the node
                        if (edge.first == edge.second){
//                            cout << "Loop ";
                            edge.second = biclique_side_to_child[1-side].at(biclique_index);
                        }

                        edge.first = child_handle;

                        edge_found = true;

//                        cout << "Case C\n";
                    } else if (edge.second == old_handle_flipped) {
                        // If this is a loop, look for the corresponding biclique child on the other side of the node
                        if (edge.first == edge.second){
//                            cout << "Loop ";
                            edge.first = biclique_side_to_child[1-side].at(biclique_index);
                        }

                        edge.second = gfa_graph.flip(child_handle);

                        edge_found = true;

//                        cout << "Case D\n";
                    }
                }

                if (edge_found) {
//                    cout << "Creating (" << gfa_graph.get_id(edge.first);
//                    cout << (gfa_graph.get_is_reverse(edge.first) ? "-" : "+");
//                    cout << ") -> (" << gfa_graph.get_id(edge.second);
//                    cout << (gfa_graph.get_is_reverse(edge.second) ? "-" : "+") << ")" << '\n';
                    overlaps.update_edge(old_edge, edge);
                    gfa_graph.create_edge(edge);
                }
            }
        }
    }
}


void Duplicator::remove_participating_edges(
        const array <deque <size_t>, 2>& sorted_bicliques_per_side,
        MutablePathDeletableHandleGraph& gfa_graph,
        nid_t parent_node
){
    for (auto side: {0,1}) {
        for (auto& biclique_index: sorted_bicliques_per_side[side]) {
            for (auto& edge: bicliques[biclique_index]) {
                if (gfa_graph.get_id(edge.first) == parent_node or gfa_graph.get_id(edge.second) == parent_node) {
//                    cout << "Deleting (" << gfa_graph.get_id(edge.first);
//                    cout << (gfa_graph.get_is_reverse(edge.first) ? "-" : "+");
//                    cout << ") -> (" << gfa_graph.get_id(edge.second);
//                    cout << (gfa_graph.get_is_reverse(edge.second) ? "-" : "+") << ")" << '\n';
                    gfa_graph.destroy_edge(edge);
                }
            }
        }
    }
}


void Duplicator::duplicate_termini(
        MutablePathDeletableHandleGraph& gfa_graph,
        array <deque <size_t>, 2>& sorted_sizes_per_side,
        const array <deque <size_t>, 2>& sorted_bicliques_per_side,
        array<map<size_t, handle_t>, 2>& biclique_side_to_child,
        const NodeInfo& node_info
){
    handle_t parent_handle = gfa_graph.get_handle(node_info.node_id, 0);

    deque<handle_t> left_children;
    deque<handle_t> right_children;

    remove_participating_edges(sorted_bicliques_per_side, gfa_graph, node_info.node_id);

    // Do left duplication
    if (not sorted_sizes_per_side[0].empty()) {
        duplicate_prefix(gfa_graph, sorted_sizes_per_side[0], left_children, parent_handle);
    } else {
        // If no duplication was performed, then the children are the parent
        left_children = {parent_handle, parent_handle};
    }

    // Update edge repair info
    for (size_t i = 0; i < sorted_bicliques_per_side[0].size(); i++) {
        auto& biclique = sorted_bicliques_per_side[0][i];
        auto& child = left_children.at(i + 1);

        biclique_side_to_child[0][biclique] = child;
    }

    // Check if right dupe is necessary
    bool right_dupe_is_trivial = false;
    if (not sorted_sizes_per_side[1].empty()) {
        if (sorted_sizes_per_side[1].size() == 1 and
            sorted_sizes_per_side[1][0] == gfa_graph.get_length(left_children.front())) {
            right_dupe_is_trivial = true;
        }
    } else {
        right_dupe_is_trivial = true;
    }

    // Do right duplication if necessary
    if (not right_dupe_is_trivial) {
        duplicate_suffix(gfa_graph, sorted_sizes_per_side[1], right_children, left_children.front());
    } else {
        // If no duplication was performed, then the children are the parent (spooky)
        right_children = {left_children.front(), left_children.front()};
    }

    // Update edge repair info
    for (size_t i = 0; i < sorted_bicliques_per_side[1].size(); i++) {
        auto& biclique = sorted_bicliques_per_side[1][i];
        auto& child = right_children[i + 1];

        biclique_side_to_child[1][biclique] = child;
    }

    // Update provenance maps
    for (size_t i=1; i<left_children.size(); i++){
        auto child_node = gfa_graph.get_id(left_children[i]);
        child_to_parent[child_node] = node_info.node_id;
        parent_to_children[node_info.node_id].emplace(child_node);
    }

    for (size_t i=1; i<right_children.size(); i++){
        auto child_node = gfa_graph.get_id(right_children[i]);
        child_to_parent[child_node] = node_info.node_id;
        parent_to_children[node_info.node_id].emplace(child_node);
    }

}


void Duplicator::duplicate_all_node_termini(MutablePathDeletableHandleGraph& gfa_graph){
    for (size_t node_id=1; node_id<node_to_biclique_edge.size(); node_id++){
        // Factor the overlaps into hierarchy: side -> biclique -> (overlap, length)
        const NodeInfo node_info(node_to_biclique_edge, bicliques, gfa_graph, overlaps, node_id);

        // Keep track of which biclique is in which position once sorted
        array <deque <size_t>, 2> sorted_sizes_per_side;
        array <deque <size_t>, 2> sorted_bicliques_per_side;

        node_info.get_sorted_biclique_extents(sorted_sizes_per_side, sorted_bicliques_per_side);

        // When repairing edges, which biclique corresponds to which child node?
        array<map<size_t, handle_t>, 2> biclique_side_to_child;

        handle_t parent_handle = gfa_graph.get_handle(node_id, 0);
        handle_t parent_handle_flipped = gfa_graph.flip(parent_handle);

        // If this is an overlapping overlap node, duping is different (and simpler)
        if (not (sorted_sizes_per_side[0].empty() and not sorted_sizes_per_side[1].empty())
                and (sorted_sizes_per_side[0][0] > gfa_graph.get_length(parent_handle) - sorted_sizes_per_side[1][0])){

            duplicate_overlapping_termini(
                    gfa_graph,
                    sorted_sizes_per_side,
                    sorted_bicliques_per_side,
                    biclique_side_to_child,
                    node_info);
        }
        else {
            duplicate_termini(
                    gfa_graph,
                    sorted_sizes_per_side,
                    sorted_bicliques_per_side,
                    biclique_side_to_child,
                    node_info);
        }

        repair_edges(
            gfa_graph,
            biclique_side_to_child,
            parent_handle,
            parent_handle_flipped);

        if (gfa_graph.get_node_count() < 30){
            string test_path_prefix = "test_bluntify_duplication_" + std::to_string(node_info.node_id);
            handle_graph_to_gfa(gfa_graph, test_path_prefix + ".gfa");
            string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                             + test_path_prefix + ".png";
            run_command(command);
        }

    }
}


}