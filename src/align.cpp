#include "handle_to_gfa.hpp"
#include "unchop.hpp"
#include "align.hpp"

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



void convert_spoa_to_bdsg(
        const HandleGraph& gfa_graph,
        Subgraph& subgraph,
        Graph& spoa_graph){

    auto& paths = spoa_graph.sequences();
    unordered_map <uint32_t, handle_t> nodes_created;
    handle_t previous_subgraph_handle;

    for (size_t side: {0,1}){
        for (auto& item: subgraph.paths_per_handle[side]){
            auto& gfa_handle = item.first;
            PathInfo& path_info = item.second;

            // This points to the first SPOA node within the path that this sequence aligned to in the SPOA graph
            auto node = paths[path_info.spoa_id];

            size_t base_index = 0;

            while (true) {
                // Check if this node has already been copied to the BDSGraph
                auto iter = nodes_created.find(node->id);

                if (iter == nodes_created.end()) {
                    char base = gfa_graph.get_base(gfa_handle, base_index);

                    auto new_subgraph_handle = subgraph.graph.create_handle(string(1, base));
                    nodes_created.emplace(node->id, new_subgraph_handle);

                    if (base_index > 0) {
                        subgraph.graph.create_edge(previous_subgraph_handle, new_subgraph_handle);
                    }

                    previous_subgraph_handle = new_subgraph_handle;
                }
                else{
                    if (base_index > 0) {
                        subgraph.graph.create_edge(previous_subgraph_handle, iter->second);
                    }

                    previous_subgraph_handle = iter->second;
                }

                subgraph.graph.append_step(path_info.path_handle, previous_subgraph_handle);

                base_index++;

                // Check if the spoa path has ended
                if (!(node = node->Successor(path_info.spoa_id))) {
                    break;
                }
            }
        }
    }
}


void add_alignments_to_poa(
        const HandleGraph& gfa_graph,
        Subgraph& subgraph,
        Graph& spoa_graph,
        unique_ptr<AlignmentEngine>& alignment_engine,
        const vector<edge_t>& biclique){

    // Since alignment may be done twice (for iterative POA), path data might need to be cleared
    subgraph.paths_per_handle[0].clear();
    subgraph.paths_per_handle[1].clear();

    // If the graph already has some sequences in it, then start the id at that number
    uint32_t spoa_id = spoa_graph.sequences().size();

    for (auto& edge: biclique){
        if (subgraph.paths_per_handle[0].count(edge.first) == 0){

            string path_name = to_string(gfa_graph.get_id(edge.first)) + "_" + to_string(0);

            path_handle_t path_handle;

            // Paths might exist from previous alignment, but will be empty
            if (not subgraph.graph.has_path(path_name)){
                path_handle = subgraph.graph.create_path_handle(path_name);
            }
            else{
                path_handle = subgraph.graph.get_path_handle(path_name);
            }

            PathInfo path_info(path_handle, spoa_id++, 0);

            subgraph.paths_per_handle[0].emplace(edge.first, path_info);
            auto sequence = gfa_graph.get_sequence(edge.first);

            auto alignment = alignment_engine->Align(sequence, spoa_graph);
            spoa_graph.AddAlignment(alignment, sequence);
        }
        if (subgraph.paths_per_handle[1].count(edge.second) == 0){
            string path_name = to_string(gfa_graph.get_id(edge.second)) + "_" + to_string(1);

            path_handle_t path_handle;

            // Paths might exist from previous alignment, but will be empty
            if (not subgraph.graph.has_path(path_name)){
                path_handle = subgraph.graph.create_path_handle(path_name);
            }
            else{
                path_handle = subgraph.graph.get_path_handle(path_name);
            }

            PathInfo path_info(path_handle, spoa_id++, 1);

            subgraph.paths_per_handle[1].emplace(edge.second, path_info);

            auto sequence = gfa_graph.get_sequence(edge.second);

            auto alignment = alignment_engine->Align(sequence, spoa_graph);
            spoa_graph.AddAlignment(alignment, sequence);
        }
    }
}


void align_biclique_overlaps(
        size_t i,
        const HandleGraph& gfa_graph,
        const Bicliques& bicliques,
        vector <Subgraph>& subgraphs){

    // TODO: switch to fetch_add atomic
    const auto& biclique = bicliques[i];
    auto& subgraph = subgraphs[i];

    if (biclique.empty()){
        return;
    }

    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 5, -3, -3, -1);

    spoa::Graph spoa_graph{};

    add_alignments_to_poa(gfa_graph, subgraph, spoa_graph, alignment_engine, biclique);

    std::cout << '\n';

    auto consensus = spoa_graph.GenerateConsensus();

//    std::cout << ">Consensus: " << consensus << std::endl;

    spoa::Graph seeded_spoa_graph{};

    auto seeded_alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 6, -2, -4, -1);

    auto alignment = seeded_alignment_engine->Align(consensus, seeded_spoa_graph);
    seeded_spoa_graph.AddAlignment(alignment, consensus);

    // Iterate a second time on alignment, this time with consensus as the seed
    add_alignments_to_poa(gfa_graph, subgraph, spoa_graph, alignment_engine, biclique);

    {
        auto seeded_consensus = seeded_spoa_graph.GenerateConsensus();

//        string prefix = "spoa_overlap_" + to_string(i);
//        seeded_spoa_graph.PrintDot(prefix + ".dot");
//
//        string command = "dot -Tpng " + prefix + ".dot -o " + prefix + ".png";
//        run_command(command);

        auto msa = seeded_spoa_graph.GenerateMultipleSequenceAlignment();

        for (const auto& it : msa) {
            std::cout << it << std::endl;
        }
        std::cout << '\n';
    }

    convert_spoa_to_bdsg(gfa_graph, subgraph, spoa_graph);

//    if (subgraph.graph.get_node_count() < 200){
//        string test_path_prefix = "test_bluntify_subgraph_" + std::to_string(i);
//        handle_graph_to_gfa(subgraph.graph, test_path_prefix + ".gfa");
//        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
//                         + test_path_prefix + ".png";
//        run_command(command);
//    }

    unchop(&subgraph.graph);

    if (subgraph.graph.get_node_count() < 200){
        string test_path_prefix = "test_bluntify_subgraph_unchopped_" + std::to_string(i);
        handle_graph_to_gfa(subgraph.graph, test_path_prefix + ".gfa");
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

}



/// For all the edges in a biclique, reorient them by matching the majority orientation of the node with the most
/// edges. In the case where there is no such orientation, pick arbitrarily
void harmonize_biclique_orientations(HandleGraph& gfa_graph, Bicliques& bicliques){
    for (auto& biclique: bicliques.bicliques){
        if (biclique.size() < 2){
            continue;
        }

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

            // If the node has not been traversed yet, add a data container
            if (iter_left == n_edges_per_node.end()){
                array<uint64_t,2> a = {0,0};
                iter_left = n_edges_per_node.emplace(left_node, a).first;
            }
            if (iter_right == n_edges_per_node.end()){
                array<uint64_t,2> a = {0,0};
                iter_right = n_edges_per_node.emplace(right_node, a).first;
            }

            // Update the counts for F and R orientations
            iter_left->second[gfa_graph.get_is_reverse(edge.first)]++;
            iter_right->second[gfa_graph.get_is_reverse(edge.second)]++;

            uint64_t total_left = iter_left->second[0] + iter_left->second[1];
            uint64_t total_right = iter_right->second[0] + iter_right->second[1];

            // Update the max observed edges
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

//        std::cout << max_node << " " << n_forward << " " << n_reverse << '\n';

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
            bool secondary_reversal;

            if (gfa_graph.get_is_reverse(seed_handle) != seed_majority_reversal){
                flipped_edge.first = gfa_graph.flip(edge.second);
                flipped_edge.second = gfa_graph.flip(edge.first);
                secondary_reversal = gfa_graph.get_is_reverse(get_side(flipped_edge, seed_side));

            }
            else{
                flipped_edge = edge;
                secondary_reversal = gfa_graph.get_is_reverse(get_side(flipped_edge, !seed_side));
            }

            // Find whether the secondary node would be reversed by this operation

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

                if (gfa_graph.get_is_reverse(secondary_seed_handle) != secondary_reversal) {
                    edge_t secondary_flipped_edge;
                    secondary_flipped_edge.first = gfa_graph.flip(secondary_edge.second);
                    secondary_flipped_edge.second = gfa_graph.flip(secondary_edge.first);

                    biclique[secondary_edge_index] = secondary_flipped_edge;
                }
            }

            biclique[i] = flipped_edge;
        }
    }
}

}
