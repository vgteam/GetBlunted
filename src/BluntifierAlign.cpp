#include "Bluntifier.hpp"
#include "handle_to_gfa.hpp"
#include "unchop.hpp"
#include <array>
#include <map>
#include "sdsl/bits.hpp"

using handlegraph::HandleGraph;
using handlegraph::nid_t;

using std::unordered_map;
using std::unique_ptr;
using std::to_string;
using std::string;
using std::array;
using std::map;

using handlegraph::nid_t;


namespace bluntifier{


handle_t& get_side(edge_t& e, bool side){
    if (side == 0){
        return e.first;
    }
    else{
        return e.second;
    }
}



void Bluntifier::convert_spoa_to_bdsg(Graph& spoa_graph, size_t i){
    auto& paths = spoa_graph.sequences();
    unordered_map <uint32_t, handle_t> nodes_created;
    handle_t previous_subgraph_handle;

    for (size_t side: {0,1}){
        for (auto& item: subgraphs[i].paths_per_handle[side]){
            auto& gfa_handle = item.first;
            PathInfo& path_info = item.second;

            // This points to the first SPOA node within the path that this sequence aligned to in the SPOA graph
            auto node = paths[path_info.poa_id];

            size_t base_index = 0;

            while (true) {
                // Check if this node has already been copied to the BDSGraph
                auto iter = nodes_created.find(node->id);

                if (iter == nodes_created.end()) {
                    char base = gfa_graph.get_base(gfa_handle, base_index);

                    auto new_subgraph_handle = subgraphs[i].graph.create_handle(string(1, base));
                    nodes_created.emplace(node->id, new_subgraph_handle);

                    if (base_index > 0) {
                        subgraphs[i].graph.create_edge(previous_subgraph_handle, new_subgraph_handle);
                    }

                    previous_subgraph_handle = new_subgraph_handle;
                }
                else{
                    if (base_index > 0) {
                        subgraphs[i].graph.create_edge(previous_subgraph_handle, iter->second);
                    }

                    previous_subgraph_handle = iter->second;
                }

                subgraphs[i].graph.append_step(path_info.path_handle, previous_subgraph_handle);

                base_index++;

                // Check if the spoa path has ended
                if (!(node = node->Successor(path_info.poa_id))) {
                    break;
                }
            }
        }
    }
}


void Bluntifier::add_alignments_to_poa(
        const function<void(const string&,const string&)>& add_alignment,
        uint32_t first_seq_id,
        size_t i){

    // Since alignment may be done twice (for iterative POA), path data might need to be cleared
    subgraphs[i].paths_per_handle[0].clear();
    subgraphs[i].paths_per_handle[1].clear();

    // If the graph already has some sequences in it, then start the id at that number
    uint32_t poa_id = first_seq_id;

    // Treating the biclique subgraph as an object with sides, add alignments, and keep track of the paths through which
    // the left and right handles traverse so they can be used for splicing later
    for (auto& edge: bicliques[i]){
        for (int side : {0, 1}) {
            handle_t handle = side == 0 ? edge.first : edge.second;
            if (subgraphs[i].paths_per_handle[side].count(handle) == 0) {
                
                string path_name = to_string(gfa_graph.get_id(handle)) + "_" + to_string(side);
                
                path_handle_t path_handle;
                
                // Paths might exist from previous alignment, but will be empty because we
                // fill it during the downstream conversion to handle graph
                if (!subgraphs[i].graph.has_path(path_name)){
                    path_handle = subgraphs[i].graph.create_path_handle(path_name);
                }
                else{
                    path_handle = subgraphs[i].graph.get_path_handle(path_name);
                }
                
                PathInfo path_info(path_handle, poa_id++, side);
                
                subgraphs[i].paths_per_handle[side].emplace(handle, path_info);
                auto sequence = gfa_graph.get_sequence(handle);
                
                add_alignment(sequence, path_name);
            }
        }
    }
}


bool Bluntifier::biclique_overlaps_are_exact(size_t i){
    bool exact = true;

    set<size_t> sizes;

    // Check if any of the edges are NOT exact overlaps
    for (auto& edge: bicliques[i]){
        auto iter = overlaps.canonicalize_and_find(edge, gfa_graph);
//        edge = iter->first;

        if (iter == overlaps.overlaps.end()){
            throw runtime_error("ERROR: edge not found in overlaps: "
                                + to_string(gfa_graph.get_id(edge.first)) + "->"
                                + to_string(gfa_graph.get_id(edge.second)));
        }
        Alignment& alignment = iter->second;

        if (not (alignment.operations.size() == 1 and alignment.operations[0].type() == 'M')){
           exact = false;
           break;
        }
        else{
            auto query_start = gfa_graph.get_length(edge.first) - alignment.operations[0].length;
            auto ref_start = 0;

            auto explicit_cigar_operations = alignment.explicitize_mismatches(gfa_graph, edge, ref_start, query_start);

            sizes.emplace(explicit_cigar_operations[0].length);

            if (not (explicit_cigar_operations.size() == 1 and explicit_cigar_operations[0].type() == '=')){
                exact = false;
                break;
            }
        }
    }

    // Only return true for bicliques that have all the same size overlaps (could be extended to more cases later)
    if (sizes.size() != 1){
        exact = false;
    }

    return exact;
}

bool Bluntifier::biclique_overlaps_are_short(size_t i, size_t max_len) {
    
    bool are_short = true;
    for (auto& edge: bicliques[i]) {
        
        pair<size_t, size_t> lengths;
        overlaps.canonicalize_and_compute_lengths(lengths, edge, gfa_graph);
        if (lengths.first > max_len || lengths.second > max_len) {
            are_short = false;
            break;
        }
    }
    return are_short;
}


void Bluntifier::create_exact_subgraph(size_t i) {
    // For an exact biclique the alignment is trivial, just pick one of the suffixes/prefixes
    auto sequence = gfa_graph.get_sequence(bicliques[i][0].first);
    auto new_subgraph_handle = subgraphs[i].graph.create_handle(sequence);

    // Treating the biclique subgraph as an object with sides, add alignments, and keep track of the paths through which
    // the left and right handles traverse so they can be used for splicing later
    for (auto& edge: bicliques[i]) {
        for (size_t side: {0, 1}) {
            handle_t h;
            if (side == 0) {
                h = edge.first;
            } else {
                h = edge.second;
            }

            string path_name = to_string(gfa_graph.get_id(h)) + "_" + to_string(side);

            path_handle_t path_handle;

            // The same handle can branch into multiple edges within a biclique, so dont add it twice
            if (not subgraphs[i].graph.has_path(path_name)) {
                path_handle = subgraphs[i].graph.create_path_handle(path_name);
                subgraphs[i].graph.append_step(path_handle, new_subgraph_handle);

                PathInfo path_info(path_handle, 0, side);
                subgraphs[i].paths_per_handle[side].emplace(h, path_info);
            }
        }
    }
}


void Bluntifier::align_biclique_overlaps(size_t i){

    
    // TODO: switch to fetch_add atomic

    // Skip trivial bicliques
    if (bicliques[i].empty()){
        return;
    }

    if (biclique_overlaps_are_exact(i)){
        // we can infer exact overlaps without alignment (helps for de Bruijn graphs)
        create_exact_subgraph(i);
    }
    else {
        // we need to align the sequences to produce the subgraph
        
        if (biclique_overlaps_are_short(i, 250)) {
            // the alignments are short enough that we can use SPOA
            
            auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 5, -3, -3, -1);
            spoa::Graph spoa_graph{};
            
            add_alignments_to_poa([&](const string& sequence,
                                      const string& name) {
                auto alignment = alignment_engine->Align(sequence, spoa_graph);
                spoa_graph.AddAlignment(alignment, sequence);
            }, spoa_graph.sequences().size(), i);
            
            // use initial POA to get a consensus sequence
            auto consensus = spoa_graph.GenerateConsensus();
            
            // Seed the alignment with the consensus
            spoa::Graph seeded_spoa_graph{};
            auto seeded_alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 6, -2, -4, -1);
            auto cons_alignment = seeded_alignment_engine->Align(consensus, seeded_spoa_graph);
            seeded_spoa_graph.AddAlignment(cons_alignment, consensus);
            
            // Iterate a second time on alignment, this time with consensus as the seed
            add_alignments_to_poa([&](const string& sequence,
                                      const string& name) {
                auto alignment = alignment_engine->Align(sequence, seeded_spoa_graph);
                seeded_spoa_graph.AddAlignment(alignment, sequence);
            }, seeded_spoa_graph.sequences().size(), i);
            
            convert_spoa_to_bdsg(seeded_spoa_graph, i);
        }
        else {
            // we use abPOA for long alignments
            
            abpoa_t* abpoa = align_with_abpoa(i);
            
            convert_abpoa_to_bdsg(abpoa, i);
            
            abpoa_free(abpoa);
        }
        
        // make long nodes from short ones
        unchop(&subgraphs[i].graph);
    }
}

abpoa_t* Bluntifier::align_with_abpoa(size_t i) {
    
    // this is annoyingly not accessible by any header, so i'll just copy it
    static const uint8_t ab_nt4_table[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };
    
    abpoa_t* abpoa = abpoa_init();
    
    // instead of doing the alignment here, just use the lambda to queue up the
    // sequences
    vector<uint8_t*> encoded_seqs;
    vector<char*> seq_names;
    vector<int> seq_lens;
    function<void(const string&,const string&)> add_alignment_to_abpoa = [&](const string& sequence,
                                                                             const string& name) {
        // convert sequence to 0-4 vals and C array
        uint8_t* enc_seq = (uint8_t*) malloc(sequence.size() * sizeof(uint8_t));
        for (size_t j = 0; j < sequence.size(); ++j) {
            enc_seq[j] = ab_nt4_table[(uint8_t)sequence[j]];
        }
        // make a C-style string for the name
        char* seq_name = (char*) malloc((name.size() + 1) * sizeof(char));
        for (size_t j = 0; j < name.size(); ++j) {
            seq_name[j] = name[j];
        }
        seq_name[name.size()] = '\0';
        
        // remember the results
        encoded_seqs.push_back(enc_seq);
        seq_names.push_back(seq_name);
        seq_lens.push_back(sequence.size());
    };
        
    // execute the lambda
    add_alignments_to_poa(add_alignment_to_abpoa, 0, i);

    
    // and now actually perform the MSA now that we've added everything
    uint8_t** dummy_cons_seq;
    int** dummy_cons_cov;
    int* dummy_cons_len;
    int dummy_cons_num;
    abpoa_msa(abpoa,
              abpoa_params,
              encoded_seqs.size(),
              seq_names.data(),
              seq_lens.data(),
              encoded_seqs.data(),
              NULL,
              &dummy_cons_seq,
              &dummy_cons_cov,
              &dummy_cons_len,
              &dummy_cons_num,
              NULL,
              NULL);
    
    // discard the consensus sequene info that we don't actually care about
    for (int j = 0; j < dummy_cons_num; ++j) {
        free(dummy_cons_seq[j]);
        free(dummy_cons_cov[j]);
    }
    free(dummy_cons_len);
    free(dummy_cons_seq);
    free(dummy_cons_cov);
    
    // clean up C arrays
    for (auto enc_seq : encoded_seqs) {
        free(enc_seq);
    }
    for (auto name : seq_names) {
        free(name);
    }
    
    return abpoa;
}

void Bluntifier::convert_abpoa_to_bdsg(abpoa_t* abpoa, size_t i) {
    
    // heavily based on abPOA's abpoa_generate_gfa function
    
    // annoyingly not available in header, so we copy it here
    static const char ab_nt256_table[256] = {
        'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
        'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
    };
    
    unordered_map<int, handle_t> abpoa_node_to_handle;
    
    auto& graph = subgraphs[i].graph;
    abpoa_graph_t* ab_graph = abpoa->abg;
    
    // get all of the nodes
    // note: 0 and 1 are reserved for the source and sink nodes (which are only
    // for internal use, no associated sequence)
    for (int j = 2; j < ab_graph->node_n; ++j) {
        string seq(1, ab_nt256_table[ab_graph->node[j].base]);
        handle_t handle = graph.create_handle(seq);
        abpoa_node_to_handle[j] = handle;
    }
    // get all of the edges
    for (int j = 2; j < ab_graph->node_n; ++j) {
        for (int k = 0; k < ab_graph->node[j].out_edge_n; ++k) {
            int next_id = ab_graph->node[j].out_id[k];
            if (next_id == ABPOA_SRC_NODE_ID || next_id == ABPOA_SINK_NODE_ID) {
                continue;
            }
            graph.create_edge(abpoa_node_to_handle[j],
                              abpoa_node_to_handle[next_id]);
        }
    }
    
    // figure out the correspondence between abpoa's integer ids and the path
    // handles that we added during the init step
    abpoa_seq_t* ab_seqs = abpoa->abs;
    vector<path_handle_t> read_id_to_path(ab_seqs->n_seq);
    for (uint64_t j = 0; j < ab_seqs->n_seq; ++j) {
        assert(ab_seqs->name[j].l != 0);
        
        string path_name(ab_seqs->name[j].s);
        assert(graph.has_path(path_name));
        read_id_to_path[j] = graph.get_path_handle(path_name);
    }
    
    // Kahn's algorithm to compute topological order
    vector<int> in_degree(ab_graph->node_n);
    vector<int> stack;
    for (int j = 0; j < ab_graph->node_n; ++j) {
        in_degree[j] = ab_graph->node[j].in_edge_n;
        if (in_degree[j] == 0) {
            stack.push_back(j);
        }
    }
    vector<int> topological_order;
    topological_order.reserve(in_degree.size());
    while (!stack.empty()) {
        int here = stack.back();
        stack.pop_back();
        topological_order.push_back(here);
        
        for (int j = 0; j < ab_graph->node[here].out_edge_n; ++j) {
            int next = ab_graph->node[here].out_id[j];
            --in_degree[next];
            if (in_degree[next] == 0) {
                stack.push_back(next);
            }
        }
    }
    
    // add paths for the input sequences
    for (int nid : topological_order) {
        if (nid == ABPOA_SRC_NODE_ID || nid == ABPOA_SINK_NODE_ID) {
            continue;
        }
        
        // add node id to read path
        // (directly copied from the abpoa function, i have no idea how this works)
        int b = 0;
        for (size_t j = 0; j < ab_graph->node[nid].read_ids_n; ++j) {
            uint64_t num = ab_graph->node[nid].read_ids[j];
            while (num) {
                uint64_t tmp = num & -num;
                int read_id = sdsl::bits::hi(tmp); // replace ilog2_64 from abpoa
                graph.append_step(read_id_to_path[b+read_id], abpoa_node_to_handle.at(nid));
                num ^= tmp;
            }
            b += 64;
        }
    }
}

/// For all the edges in a biclique, reorient them by matching the majority orientation of the node with the most
/// edges. In the case where there is no such orientation, pick arbitrarily
void Bluntifier::harmonize_biclique_orientations(){
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
