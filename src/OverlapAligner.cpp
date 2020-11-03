#include "OverlapAligner.hpp"
#include <array>
#include <map>

using handlegraph::nid_t;
using std::array;
using std::map;


namespace bluntifier{


handle_t& get_side(edge_t& e, bool side){
    if (side == 0){
        return e.first;
    }
    else{
        return e.second;
    }
}


/// For all the edges in a biclique, reorient them by matching the majority orientation of the node with the most
/// edges. In the case where there is no such orientation, pick arbitrarily
void harmonize_biclique_orientations(HandleGraph& gfa_graph, Bicliques& bicliques){
    for (auto& biclique: bicliques.bicliques){
        map <nid_t, array<uint64_t, 2> > n_edges_per_node;
        map <nid_t, vector<size_t> > primary_edge_indexes;

        uint64_t max_edges = 0;
        nid_t max_node;

        for (size_t edge_index=0; edge_index<biclique.size(); edge_index++){
            auto& edge = biclique[edge_index];

            nid_t left_node = gfa_graph.get_id(edge.first);
            nid_t right_node = gfa_graph.get_id(edge.second);

            auto iter_left = n_edges_per_node.find(left_node);
            auto iter_right = n_edges_per_node.find(right_node);

            if (iter_left == n_edges_per_node.end()){
                array<uint64_t,2> a = {0,0};
                iter_left = n_edges_per_node.emplace(left_node, a).first;
            }

            if (iter_right == n_edges_per_node.end()){
                array<uint64_t,2> a = {0,0};
                iter_right = n_edges_per_node.emplace(right_node, a).first;
            }

            iter_left->second[gfa_graph.get_is_reverse(edge.first)]++;
            iter_right->second[gfa_graph.get_is_reverse(edge.second)]++;

            uint64_t total_left = iter_left->second[0] + iter_left->second[1];
            uint64_t total_right = iter_right->second[0] + iter_right->second[1];

            if (total_left > max_edges){
                max_edges = total_left;
                max_node = left_node;
            }
            if (total_right > max_edges){
                max_edges = total_right;
                max_node = right_node;
            }

            primary_edge_indexes[left_node].emplace_back(edge_index);
            primary_edge_indexes[right_node].emplace_back(edge_index);
        }

        // Decide what orientation to use for the seed node
        auto max_edge_info = n_edges_per_node.at(max_node);
        auto n_reverse = max_edge_info[1];
        auto n_forward = max_edge_info[0];

        bool seed_majority_reversal = false;
        if (n_reverse > n_forward){
            seed_majority_reversal = true;
        }

        // Iterate the edges directly linked with the seed node and flip them if necessary
        // Additionally flip (as necessary) all the secondary edges that stem from the seed node's adjacent nodes
        for (auto i: primary_edge_indexes.at(max_node)){
            auto& edge = biclique[i];
            size_t seed_side;

            if (gfa_graph.get_id(edge.first) == max_node){
                seed_side = 0;
            }
            else{
                seed_side = 1;
            }

            auto& seed_handle = get_side(edge, seed_side);
            auto& other_handle = get_side(edge, !seed_side);
            auto other_node = gfa_graph.get_id(other_handle);

            edge_t flipped_edge;

            if (gfa_graph.get_is_reverse(seed_handle) != seed_majority_reversal){
                flipped_edge.first = gfa_graph.flip(edge.second);
                flipped_edge.second = gfa_graph.flip(edge.first);
            }

            // Find whether the secondary node would be reversed by this operation
            bool secondary_reversal = gfa_graph.get_is_reverse(get_side(flipped_edge, seed_side));

            for (auto& secondary_edge_index: primary_edge_indexes.at(other_node)) {
                auto& secondary_edge = biclique[secondary_edge_index];

                if (secondary_edge == edge) {
                    // Don't reevaluate the edge that stems from the seed node
                    continue;
                }

                bool secondary_seed_side;
                if (gfa_graph.get_id(secondary_edge.first) == other_node) {
                    secondary_seed_side = 0;
                } else {
                    secondary_seed_side = 1;
                }

                auto& secondary_seed_handle = get_side(secondary_edge, secondary_seed_side);
                auto& secondary_other_handle = get_side(secondary_edge, !secondary_seed_side);

                if (gfa_graph.get_is_reverse(secondary_seed_handle) != secondary_reversal) {
                    edge_t secondary_flipped_edge;
                    secondary_flipped_edge.first = gfa_graph.flip(secondary_edge.second);
                    secondary_flipped_edge.second = gfa_graph.flip(secondary_edge.first);

                    biclique[secondary_edge_index] = secondary_flipped_edge;
                }

                biclique[i] = flipped_edge;
            }
        }
    }
}

}
