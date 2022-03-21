#include "Duplicator.hpp"

namespace bluntifier{


Duplicator::Duplicator(
        const vector <vector <BicliqueEdgeIndex> >& node_to_biclique_edge,
        OverlapMap& overlaps,
        Bicliques& bicliques,
        map <nid_t, set<nid_t> >& parent_to_children,
        map <nid_t, pair<nid_t, bool> >& child_to_parent,
        map<nid_t, OverlappingNodeInfo>& overlapping_overlap_nodes
        ):
        node_to_biclique_edge(node_to_biclique_edge),
        overlaps(overlaps),
        bicliques(bicliques),
        parent_to_children(parent_to_children),
        child_to_parent(child_to_parent),
        overlapping_overlap_nodes(overlapping_overlap_nodes)
{}


map<nid_t, OverlappingNodeInfo>::iterator Duplicator::preprocess_overlapping_overlaps(
        MutablePathDeletableHandleGraph& gfa_graph,
        array <deque <size_t>, 2>& sorted_sizes_per_side,
        array <deque <size_t>, 2>& sorted_bicliques_per_side,
        array<map<size_t, handle_t>, 2>& biclique_side_to_child,
        const NodeInfo& node_info){

    auto result = overlapping_overlap_nodes.emplace(node_info.node_id, node_info.node_id);
    auto& overlap_node_info = result.first->second;

    auto parent_handle = gfa_graph.get_handle(node_info.node_id, false);
    overlap_node_info.length = gfa_graph.get_length(parent_handle);

    // Iteratively remove and document the longest overlaps from alternating sides, until the node is effectively a
    // normal node. Additionally, duplicate orphan segments from the parent node for the overlapping overlaps
    while (contains_overlapping_overlaps(gfa_graph, parent_handle, sorted_sizes_per_side)){
        if (sorted_sizes_per_side[0][0] > sorted_sizes_per_side[1][0]){
            int64_t s = sorted_sizes_per_side[0][0];
            string sequence = gfa_graph.get_subsequence(parent_handle, 0, s);

            auto child = gfa_graph.create_handle(sequence);

            // Storing this child so the other participant in the overlap can be spliced
            size_t biclique_index = sorted_bicliques_per_side[0][0];
            biclique_side_to_child[0].emplace(biclique_index, child);

            // More info needs to be stored until AFTER all POA subgraphs are done, to enable splicing within this node
            OverlappingChild child_info(child, biclique_index, 0);
            overlap_node_info.overlapping_children[0].emplace(s-1, child_info);

            if (s-1 < 0){
                child_info.print(gfa_graph);
                throw runtime_error("ERROR: negative start index for overlapping overlap" );
            }

            // Update provenance map
            auto child_node = gfa_graph.get_id(child);
            child_to_parent[child_node] = {node_info.node_id, true};
            parent_to_children[node_info.node_id].emplace(child_node);

            sorted_sizes_per_side[0].pop_front();
            sorted_bicliques_per_side[0].pop_front();
        }
        else {
            int64_t start = gfa_graph.get_length(parent_handle) - sorted_sizes_per_side[1][0];
            string sequence = gfa_graph.get_subsequence(parent_handle, start, sorted_sizes_per_side[1][0]);

            auto child = gfa_graph.create_handle(sequence);

            // Storing this child so the other participant in the overlap can be spliced
            size_t biclique_index = sorted_bicliques_per_side[1][0];
            biclique_side_to_child[1].emplace(biclique_index, child);

            // More info needs to be stored until AFTER all POA subgraphs are done, to enable splicing within this node
            OverlappingChild child_info(child, biclique_index, 1);
            overlap_node_info.overlapping_children[1].emplace(start, child_info);

            if (start < 0){
                child_info.print(gfa_graph);
                throw runtime_error("ERROR: negative start index for overlapping overlap" );
            }

            // Update provenance map
            auto child_node = gfa_graph.get_id(child);
            child_to_parent[child_node] = {node_info.node_id, true};
            parent_to_children[node_info.node_id].emplace(child_node);

            sorted_sizes_per_side[1].pop_front();
            sorted_bicliques_per_side[1].pop_front();
        }
    }

    return result.first;
}


void Duplicator::postprocess_overlapping_overlap(
        const HandleGraph& gfa_graph,
        map<nid_t, OverlappingNodeInfo>::iterator iter,
        array<map<size_t, handle_t>, 2> biclique_side_to_child){

    auto& overlapping_node_info = iter->second;

    // Deduplicate the non-OO biclique info
    array<set<size_t>, 2> oo_biclique_indexes;
    for (auto side: {0,1}) {
        for (auto& item:overlapping_node_info.overlapping_children[side]) {
            oo_biclique_indexes[side].emplace(item.second.biclique_index);
        }
    }

    // Refactor the "biclique_side_to_child" object to match the structure of the overlapping children (mapped by pos)
    for (auto side: {0,1}) {
        for (auto& item: biclique_side_to_child[side]){
            auto biclique_index = item.first;
            if (oo_biclique_indexes[side].count(biclique_index) == 0){
                auto handle = item.second;

                size_t position;
                if (side == 0){
                    position = gfa_graph.get_length(handle) - 1;
                }
                else{
                    position = overlapping_node_info.length - gfa_graph.get_length(handle);
                }

                OverlappingChild o(handle, biclique_index, side);
                overlapping_node_info.normal_children[side].emplace(position, o);
            }
        }
    }

    overlapping_node_info.parent_path_name = to_string(overlapping_node_info.parent_node);
//    find_leftover_parent(gfa_graph, overlapping_node_info);
}


void Duplicator::repair_edges(
        MutablePathDeletableHandleGraph& gfa_graph,
        const array <map <size_t, handle_t>, 2>& biclique_side_to_child,
        handle_t old_handle,
        handle_t old_handle_flipped) {

    for (const size_t side: {0, 1}) {
        for (const auto&[biclique_index, child_handle]: biclique_side_to_child[side]) {
            for (auto& edge: bicliques[biclique_index]) {
                auto old_edge = edge;

                bool edge_found = false;
                if (side == 0) {
                    if (edge.second == old_handle) {
                        // If this is a loop, look for the corresponding biclique child on the other side of the node
                        if (edge.first == edge.second){
                            edge.first = biclique_side_to_child[1-side].at(biclique_index);
                        }

                        edge.second = child_handle;

                        edge_found = true;

                    } else if (edge.first == old_handle_flipped) {
                        // If this is a loop, look for the corresponding biclique child on the other side of the node
                        if (edge.first == edge.second){
                            edge.second = biclique_side_to_child[1-side].at(biclique_index);
                        }

                        edge.first = gfa_graph.flip(child_handle);

                        edge_found = true;
                    }
                }
                else {
                    if (edge.first == old_handle) {
                        // If this is a loop, look for the corresponding biclique child on the other side of the node
                        if (edge.first == edge.second){
                            edge.second = biclique_side_to_child[1-side].at(biclique_index);
                        }

                        edge.first = child_handle;

                        edge_found = true;

                    } else if (edge.second == old_handle_flipped) {
                        // If this is a loop, look for the corresponding biclique child on the other side of the node
                        if (edge.first == edge.second){
                            edge.first = biclique_side_to_child[1-side].at(biclique_index);
                        }

                        edge.second = gfa_graph.flip(child_handle);

                        edge_found = true;
                    }
                }

                if (edge_found) {
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
                    gfa_graph.destroy_edge(edge);
                }
            }
        }
    }
}


void Duplicator::duplicate_termini(
        MutablePathDeletableHandleGraph& gfa_graph,
        array <deque <size_t>, 2> sorted_sizes_per_side,
        const array <deque <size_t>, 2>& sorted_bicliques_per_side,
        array<map<size_t, handle_t>, 2>& biclique_side_to_child,
        const NodeInfo& node_info
){
    handle_t parent_handle = gfa_graph.get_handle(node_info.node_id, 0);

    deque<handle_t> left_children;
    deque<handle_t> right_children;

    // Do left duplication
    bool left_dupes_exist = true;
    if (not sorted_sizes_per_side[0].empty()) {
        duplicate_prefix(gfa_graph, sorted_sizes_per_side[0], left_children, parent_handle);
    } else {
        // If no duplication was performed, then the children are the parent
        left_children = {parent_handle, parent_handle};
        left_dupes_exist = false;
    }

    // Update edge repair info
    for (size_t i = 0; i < sorted_bicliques_per_side[0].size(); i++) {
        auto& biclique = sorted_bicliques_per_side[0][i];
        auto& child = left_children.at(i + 1);

        // The Overlapping Overlap preprocessing step may have already added something here
        biclique_side_to_child[0].try_emplace(biclique, child);
    }

    // Check if right dupe is necessary
    bool right_dupe_is_trivial = false;
    bool right_dupes_exist = true;
    if (not sorted_sizes_per_side[1].empty()) {
        if (sorted_sizes_per_side[1].size() == 1 and
            sorted_sizes_per_side[1][0] == gfa_graph.get_length(left_children.front())) {
            right_dupe_is_trivial = true;
        }
    } else {
        right_dupes_exist = false;
    }

    // Do right duplication if necessary
    if (not right_dupe_is_trivial and right_dupes_exist) {
        duplicate_suffix(gfa_graph, sorted_sizes_per_side[1], right_children, left_children.front());
    } else {
        // If no duplication was performed, then the children are the parent (spooky)
        right_children = {left_children.front(), left_children.front()};
    }

    // Update edge repair info
    for (size_t i = 0; i < sorted_bicliques_per_side[1].size(); i++) {
        auto& biclique = sorted_bicliques_per_side[1][i];
        auto& child = right_children[i + 1];

        // The Overlapping Overlap preprocessing step may have already added something here
        biclique_side_to_child[1].try_emplace(biclique, child);
    }

    // Update provenance maps
    if (left_dupes_exist) {
        auto parent_node = gfa_graph.get_id(left_children[0]);

        for (size_t i = 1; i < left_children.size(); i++) {
            auto child_node = gfa_graph.get_id(left_children[i]);

            child_to_parent[child_node] = {node_info.node_id, child_node != parent_node};
            parent_to_children[node_info.node_id].emplace(child_node);
        }
    }

    if (right_dupes_exist) {
        auto parent_node = gfa_graph.get_id(right_children[0]);

        for (size_t i = 1; i < right_children.size(); i++) {
            auto child_node = gfa_graph.get_id(right_children[i]);

            child_to_parent[child_node] = {node_info.node_id, child_node != parent_node};
            parent_to_children[node_info.node_id].emplace(child_node);
        }
    }

}


bool Duplicator::contains_overlapping_overlaps(
        const HandleGraph& gfa_graph,
        handle_t parent_handle,
        const array <deque <size_t>, 2>& sorted_sizes_per_side){

    bool contains_oo = false;

    // If either side is empty, it's impossible to have an Overlapping Overlap
    if ((not sorted_sizes_per_side[0].empty()) and (not sorted_sizes_per_side[1].empty())) {
        // If there are overlaps on both sides, simply test their longest overlaps for conflict
        if (sorted_sizes_per_side[0][0] > gfa_graph.get_length(parent_handle) - sorted_sizes_per_side[1][0]) {
            contains_oo = true;
        }
    }

    return contains_oo;
}


void Duplicator::duplicate_all_node_termini(MutablePathDeletableHandleGraph& gfa_graph){
    for (size_t node_id=1; node_id<node_to_biclique_edge.size(); node_id++){
        // Factor the overlaps into hierarchy: side -> biclique -> (overlap, length)
        const NodeInfo node_info(node_to_biclique_edge, bicliques, gfa_graph, overlaps, node_id);

        // Keep track of which biclique is in which position once sorted
        array <deque <size_t>, 2> sorted_sizes_per_side;
        array <deque <size_t>, 2> sorted_bicliques_per_side;

        node_info.get_sorted_biclique_extents(sorted_sizes_per_side, sorted_bicliques_per_side);

        remove_participating_edges(sorted_bicliques_per_side, gfa_graph, node_info.node_id);

        // When repairing edges, which biclique corresponds to which child node?
        array<map<size_t, handle_t>, 2> biclique_side_to_child;

        handle_t parent_handle = gfa_graph.get_handle(node_id, 0);
        handle_t parent_handle_flipped = gfa_graph.flip(parent_handle);

        // Set a path that only describes the parent node
        string parent_path_name = to_string(node_info.node_id);
        auto parent_path_handle = gfa_graph.create_path_handle(parent_path_name);
        gfa_graph.append_step(parent_path_handle, parent_handle);

        set <size_t> overlapping_bicliques;

        auto overlapping_overlap_iter = overlapping_overlap_nodes.end();
        // If this is an overlapping overlap node, the offending overlaps need to be processed separately
        if (contains_overlapping_overlaps(gfa_graph, parent_handle, sorted_sizes_per_side)){
            overlapping_overlap_iter = preprocess_overlapping_overlaps(
                    gfa_graph,
                    sorted_sizes_per_side,
                    sorted_bicliques_per_side,
                    biclique_side_to_child,
                    node_info);
        }

        duplicate_termini(
                gfa_graph,
                sorted_sizes_per_side,
                sorted_bicliques_per_side,
                biclique_side_to_child,
                node_info);

        repair_edges(
                gfa_graph,
                biclique_side_to_child,
                parent_handle,
                parent_handle_flipped);


        if (overlapping_overlap_iter != overlapping_overlap_nodes.end()){
            postprocess_overlapping_overlap(gfa_graph, overlapping_overlap_iter, biclique_side_to_child);
        }
    }
}


}