#include "Duplicator.hpp"

namespace bluntifier{


Duplicator::Duplicator(
        const vector <vector <BicliqueEdgeIndex> >& node_to_biclique_edge,
        Bicliques& bicliques,
        OverlapMap& overlaps):
        node_to_biclique_edge(node_to_biclique_edge),
        bicliques(bicliques),
        overlaps(overlaps)
{}


void Duplicator::update_biclique_edges(
        MutablePathDeletableHandleGraph& gfa_graph,
        nid_t old_node_id,
        handle_t old_handle,
        handle_t old_handle_flipped,
        const array <deque <size_t>, 2>& sorted_bicliques_per_side,
        const deque <handle_t>& children,
        bool duped_side){

    for (auto& item: children) {
        std::cout << gfa_graph.get_id(item) << " " << as_integer(item) << "F " << as_integer(gfa_graph.flip(item)) << "R " << gfa_graph.get_sequence(item) << '\n';
    }

    for (bool side: {0,1}) {
        for (size_t i = 0; i < sorted_bicliques_per_side[side].size(); i++) {
            size_t biclique_index = sorted_bicliques_per_side[side][i];

            for (auto& edge: bicliques[biclique_index]) {
                auto old_edge = edge;

                cout << "Replacing " << old_node_id << '\n';
                cout << "Replacing " << as_integer(old_handle) <<"F or " << as_integer(old_handle_flipped) << "R" << '\n';

                cout << "[" << biclique_index << "] " << duped_side << " " << side << as_integer(old_edge.first) << "h->" << as_integer(old_edge.second) << "h" << '\n';

                if (duped_side == 0) {
                    if (side == 0) {
                        if (old_edge.second == old_handle){
                            // Account for self loops (non-reversing)
                            if (old_edge.first == old_edge.second){
                                edge.first = children[0];
                            }

                            edge.second = children[i+1];
                            cout << "Creating " << gfa_graph.get_id(edge.first) << "->" << gfa_graph.get_id(edge.second) << '\n';
                            gfa_graph.create_edge(edge);
                            overlaps.update_edge(old_edge, edge);

                        }
                        else if (old_edge.second == old_handle_flipped){
                            // Account for self loops (non-reversing)
                            if (old_edge.first == old_edge.second){
                                edge.first = gfa_graph.flip(children[i+1]);
                            }

                            edge.second = gfa_graph.flip(children[0]);
                            cout << "Creating " << gfa_graph.get_id(edge.first) << "->" << gfa_graph.get_id(edge.second) << '\n';
                            gfa_graph.create_edge(edge);
                            overlaps.update_edge(old_edge, edge);

                        }
                    }
                    else{
                        if (old_edge.first == old_handle){
                            // Account for self loops (non-reversing)
                            if (old_edge.first == old_edge.second){
                                edge.second = children[i+1];
                            }

                            edge.first = children[0];
                            cout << "Creating " << gfa_graph.get_id(edge.first) << "->" << gfa_graph.get_id(edge.second) << '\n';
                            gfa_graph.create_edge(edge);
                            overlaps.update_edge(old_edge, edge);
                        }
                        else if (old_edge.first == old_handle_flipped){
                            // Account for self loops (non-reversing)
                            if (old_edge.first == old_edge.second){
                                edge.second = gfa_graph.flip(children[0]);
                            }

                            edge.first = gfa_graph.flip(children[i+1]);
                            cout << "Creating " << gfa_graph.get_id(edge.first) << "->" << gfa_graph.get_id(edge.second) << '\n';
                            gfa_graph.create_edge(edge);
                            overlaps.update_edge(old_edge, edge);
                        }
                    }
                }
                else{
                    if (side == 0) {
                        if (old_edge.second == old_handle){
                            // Account for self loops (non-reversing)
                            if (old_edge.first == old_edge.second){
                                edge.first = children[i+1];
                            }

                            edge.second = children[0];
                            cout << "Creating " << gfa_graph.get_id(edge.first) << "->" << gfa_graph.get_id(edge.second) << '\n';
                            gfa_graph.create_edge(edge);
                            overlaps.update_edge(old_edge, edge);

                        }
                        else if (old_edge.second == old_handle_flipped){
                            // Account for self loops (non-reversing)
                            if (old_edge.first == old_edge.second){
                                edge.first = gfa_graph.flip(children[0]);
                            }

                            edge.second = gfa_graph.flip(children[i+1]);
                            cout << "Creating " << gfa_graph.get_id(edge.first) << "->" << gfa_graph.get_id(edge.second) << '\n';
                            gfa_graph.create_edge(edge);
                            overlaps.update_edge(old_edge, edge);

                        }
                    }
                    else{
                        if (old_edge.first == old_handle){
                            // Account for self loops (non-reversing)
                            if (old_edge.first == old_edge.second){
                                edge.second = children[0];
                            }

                            edge.first = children[i+1];
                            cout << "Creating " << gfa_graph.get_id(edge.first) << "->" << gfa_graph.get_id(edge.second) << '\n';
                            gfa_graph.create_edge(edge);
                            overlaps.update_edge(old_edge, edge);

                        }
                        else if (old_edge.first == old_handle_flipped){
                            // Account for self loops (non-reversing)
                            if (old_edge.first == old_edge.second){
                                edge.second = gfa_graph.flip(children[i+1]);
                            }

                            edge.first = gfa_graph.flip(children[0]);
                            cout << "Creating " << gfa_graph.get_id(edge.first) << "->" << gfa_graph.get_id(edge.second) << '\n';
                            gfa_graph.create_edge(edge);
                            overlaps.update_edge(old_edge, edge);

                        }
                    }
                }
            }
        }
    }
    cout << '\n';
}


void Duplicator::remove_participating_edges(
        const array <deque <size_t>, 2>& sorted_bicliques_per_side,
        MutablePathDeletableHandleGraph& gfa_graph,
        nid_t parent_node
){

    for (bool side: {0,1}) {
        for (auto& biclique_index: sorted_bicliques_per_side[side]) {
            for (auto& edge: bicliques[biclique_index]) {
                if (gfa_graph.get_id(edge.first) == parent_node or gfa_graph.get_id(edge.second) == parent_node) {
//                    cout << "Deleting " << gfa_graph.get_id(edge.first) << "->" << gfa_graph.get_id(edge.second)
//                         << '\n';
                    gfa_graph.destroy_edge(edge);
                }
            }
        }
    }

}


void Duplicator::duplicate_termini(MutablePathDeletableHandleGraph& gfa_graph){

    for (size_t node_id=1; node_id<node_to_biclique_edge.size(); node_id++){
        // Factor the overlaps into hierarchy: side -> biclique -> (overlap, length)
        NodeInfo node_info(node_to_biclique_edge, bicliques, gfa_graph, overlaps, node_id);

        node_info.print_stats();

        // Keep track of which biclique is in which position once sorted
        array <deque <size_t>, 2> sorted_sizes_per_side;
        array <deque <size_t>, 2> sorted_bicliques_per_side;

        node_info.get_sorted_biclique_extents(sorted_sizes_per_side, sorted_bicliques_per_side);

        {
            string test_path_prefix = "test_bluntify_" + std::to_string(node_id) + "_";
            handle_graph_to_gfa(gfa_graph, test_path_prefix + ".gfa");
            string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                             + test_path_prefix + ".png";
            run_command(command);
        }

        handle_t parent_handle = gfa_graph.get_handle(node_id, 0);
        handle_t parent_handle_flipped = gfa_graph.flip(parent_handle);
        nid_t parent_node = node_id;

        remove_participating_edges(sorted_bicliques_per_side, gfa_graph, parent_node);

        deque<handle_t> left_children;
        deque<handle_t> right_children;

        if (not sorted_sizes_per_side[0].empty()) {
            duplicate_prefix(gfa_graph, sorted_sizes_per_side[0], left_children, parent_handle);

            update_biclique_edges(
                    gfa_graph,
                    parent_node,
                    parent_handle,
                    parent_handle_flipped,
                    sorted_bicliques_per_side,
                    left_children,
                    0);

            parent_handle = left_children.front();
            parent_handle_flipped = gfa_graph.flip(parent_handle);
            parent_node = gfa_graph.get_id(parent_handle);

            {
                string test_path_prefix = "test_bluntify_" + std::to_string(node_id) + "_" + std::to_string(0);
                handle_graph_to_gfa(gfa_graph, test_path_prefix + ".gfa");
                string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                                 + test_path_prefix + ".png";
                run_command(command);
            }
        }


        if (not sorted_sizes_per_side[1].empty()) {
            vector <edge_t> right_edges;
            gfa_graph.follow_edges(parent_handle, false, [&](const handle_t& h){
                edge_t e = {parent_handle, h};
                right_edges.emplace_back(e);
            });

            for (auto& e: right_edges){
                gfa_graph.destroy_edge(e);
            }

            // Skip trivial duplication
            if (sorted_sizes_per_side[1].size() == 1 and sorted_sizes_per_side[1][0] == gfa_graph.get_length(parent_handle)){
//                cout << "Skipping trivial duplication\n";
                continue;
            }

            duplicate_suffix(gfa_graph, sorted_sizes_per_side[1], right_children, parent_handle);

            update_biclique_edges(
                    gfa_graph,
                    parent_node,
                    parent_handle,
                    parent_handle_flipped,
                    sorted_bicliques_per_side,
                    right_children,
                    1);

            {
                string test_path_prefix = "test_bluntify_" + std::to_string(node_id) + "_" + std::to_string(1);
                handle_graph_to_gfa(gfa_graph, test_path_prefix + ".gfa");
                string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                                 + test_path_prefix + ".png";
                run_command(command);
            }
        }
    }
}

}