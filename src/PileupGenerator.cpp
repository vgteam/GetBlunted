#include "PileupGenerator.hpp"
#include "handle_to_gfa.hpp"
#include "unchop.hpp"

using std::max;
using handlegraph::step_handle_t;

namespace bluntifier {


BicliqueIterator::BicliqueIterator():
    first_step(true)
{}


PoaPileupGenerator::PoaPileupGenerator(
        const BipartiteGraph& bipartite_graph,
        const IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps,
        HandleGraph& gfa_graph,
        PoaPileup& pileup,
        vector <vector <SpliceData> >& splice_sites,
        vector <mutex>& splice_site_mutexes
):
    bipartite_graph(bipartite_graph),
    id_map(id_map),
    overlaps(overlaps),
    gfa_graph(gfa_graph),
    pileup(pileup),
    splice_sites(splice_sites),
    splice_site_mutexes(splice_site_mutexes)
{}


PileupGenerator::PileupGenerator(
        const BipartiteGraph& bipartite_graph,
        const IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps,
        HandleGraph& gfa_graph,
        Pileup& pileup,
        vector <vector <SpliceData> >& splice_sites,
        vector <mutex>& splice_site_mutexes
):
    bipartite_graph(bipartite_graph),
    id_map(id_map),
    overlaps(overlaps),
    gfa_graph(gfa_graph),
    pileup(pileup),
    splice_sites(splice_sites),
    splice_site_mutexes(splice_site_mutexes)
{}


bool PileupGenerator::traverse_bipartition_edges(BicliqueIterator& iterator) {

    bool done;

    if (iterator.first_step) {
        PileupGenerator::traverse_bipartition_nodes(iterator);
        if (iterator.is_left) {
            iterator.prev_left_node = iterator.node;
            iterator.prev_left_is_valid = true;
        } else {
            iterator.prev_right_node = iterator.node;
            iterator.prev_right_is_valid = true;
        }
    }

    done = not PileupGenerator::traverse_bipartition_nodes(iterator);

    if (not done) {
        if (iterator.is_left) {
            if (iterator.prev_right_is_valid) {
                iterator.edge = {iterator.node, gfa_graph.flip(iterator.prev_right_node)};
            }
            iterator.prev_left_node = iterator.node;
            iterator.prev_left_is_valid = true;
        } else {
            if (iterator.prev_left_is_valid) {
                iterator.edge = {iterator.prev_left_node, gfa_graph.flip(iterator.node)};
            }
            iterator.prev_right_node = iterator.node;
            iterator.prev_right_is_valid = true;
        }
    }

    return not done;
}


bool PileupGenerator::traverse_bipartition_nodes(BicliqueIterator& iterator){

    if (iterator.first_step){
        // Pick an arbitrary node to start with
        handle_t start;

        std::cout << "LEFT size:\t" << bipartite_graph.left_size() << '\n';
        std::cout << "RIGHT size:\t" << bipartite_graph.right_size() << '\n';
        // Choose one from the bigger set, if possible
        if (bipartite_graph.left_size() > bipartite_graph.right_size()) {
            std::cout << "STARTING traversal from LEFT side\n";
            start = *bipartite_graph.left_begin();
        }
        else{
            start = *bipartite_graph.right_begin();
            std::cout << "STARTING traversal from RIGHT side\n";
        }

        iterator.node_stack.push(start);
        iterator.first_step = false;
    }

    while (iterator.visited.size() < bipartite_graph.left_size() + bipartite_graph.right_size()) {
        // Emulate DFS
        handle_t node = iterator.node_stack.top();
        iterator.node_stack.pop();
        if (iterator.visited.count(node) == 0) {
            iterator.node = node;
            iterator.visited.insert(node);

            // Just search the other side of the partition
            if (bipartite_graph.is_left_side(node)) {
                iterator.is_left = true;
                for (auto it = bipartite_graph.right_begin(); it != bipartite_graph.right_end(); ++it) {
                    iterator.node_stack.push(*it);
                }
            } else {
                iterator.is_left = false;
                for (auto it = bipartite_graph.left_begin(); it != bipartite_graph.left_end(); ++it) {
                    iterator.node_stack.push(*it);
                }
            }
            return true;
        }
    }
    return false;
}


void PileupGenerator::debug_print(BicliqueIterator& biclique_iterator){

    auto iter = overlaps.canonicalize_and_find(biclique_iterator.edge, gfa_graph);
    const edge_t& edge = iter->first;
    pair<size_t, size_t> lengths;
    size_t start;

    std::cout << id_map.get_name(gfa_graph.get_id(edge.first)) << "->";
    std::cout << id_map.get_name(gfa_graph.get_id(edge.second)) << '\n';

    iter->second.compute_lengths(lengths);

    start = gfa_graph.get_length(edge.first) - lengths.first;

    std::cout << "Is left: " << biclique_iterator.is_left << '\n';
    std::cout << lengths.first << " " << lengths.second << '\n';
    std::cout << iter->second.create_formatted_alignment_string(gfa_graph, edge, start, 0) << '\n';
}


void PileupGenerator::update_pseudoref(
        handle_t pseudo_reference,
        size_t pseudo_ref_length,
        size_t prev_pseudo_ref_length,
        int64_t pseudo_ref_id,
        bool is_left) {

    // Walk from the last index of the prev pseudoref up to the current index
    if (is_left) {
        for (size_t i = prev_pseudo_ref_length; i < pseudo_ref_length; i++) {
            auto h = pileup.graph.create_handle(string(1, gfa_graph.get_base(pseudo_reference, i)));
            std::cout << "created node:\t" << gfa_graph.get_id(h) << " " << gfa_graph.has_node(gfa_graph.get_id(h)) << '\n';

            if (not pileup.paths[pseudo_ref_id].empty()) {
                pileup.graph.create_edge(pileup.paths[pseudo_ref_id].back(), h);
                std::cout << "created edge:\t"
                          << pileup.graph.get_id(pileup.paths[pseudo_ref_id].back()) << "->"
                          << pileup.graph.get_id(h) << " "
                          << pileup.graph.has_edge(pileup.paths[pseudo_ref_id].back(), h) << '\n';
            }
            pileup.paths[pseudo_ref_id].push_back(h);

            std::cout << "base:\t" << gfa_graph.get_base(pseudo_reference, i) << '\n';
            std::cout << "id:\t" << gfa_graph.get_id(h) << '\n' << '\n';
        }
    }
    // Walk BACKWARDS from the last index of the prev pseudoref down to the current index
    else {
        size_t start = gfa_graph.get_length(pseudo_reference) - prev_pseudo_ref_length - 1;
        size_t stop = gfa_graph.get_length(pseudo_reference) - pseudo_ref_length - 1;

        std::cout << gfa_graph.get_length(pseudo_reference) << " " << pileup.paths[pseudo_ref_id].size() << " " << start << " " << stop << '\n';

        for (int64_t i = start; i > int64_t(stop); i--) {
            auto h = pileup.graph.create_handle(string(1, gfa_graph.get_base(pseudo_reference, i)));
            std::cout << "created node:\t" << gfa_graph.get_id(h) << " " << gfa_graph.has_node(gfa_graph.get_id(h)) << '\n';

            if (not pileup.paths[pseudo_ref_id].empty()) {
                pileup.graph.create_edge(h, pileup.paths[pseudo_ref_id].front());
                std::cout << "created edge:\t"
                          << pileup.graph.get_id(h) << "->"
                          << pileup.graph.get_id(pileup.paths[pseudo_ref_id].front()) << " "
                          << pileup.graph.has_edge(h, pileup.paths[pseudo_ref_id].front()) << '\n';

            }
            pileup.paths[pseudo_ref_id].push_front(h);
        }
    }
}


/// If the pseudoref is on the right side of the current edge, then everything
/// that depends on the expected relationship between ref/query must be flipped.
/// This method tells whether that is the case
bool PileupGenerator::pseudoref_is_reversed(
        const edge_t& canonical_edge,
        const BicliqueIterator& biclique_iterator){

    // If this is the first edge in the traversal, just use the directionality of the
    // traversal to pick the pseudoref. The "node" field in bicliqueIterator tells which
    // is the sink node
    if (pileup.edges_traversed.empty()){
        if (biclique_iterator.node == canonical_edge.second){
            return false;
        }
        else{
            return true;
        }
    }

    // Walk backward through the traversed edges and stop when a previously added node
    // matches one of the current nodes in the edge. That node is the pseudoref.
    for (int64_t i=pileup.edges_traversed.size()-1; i >= 0; i--){
        auto prev_edge = pileup.edges_traversed[i];
        if (prev_edge.first == canonical_edge.first or prev_edge.second == canonical_edge.first){
            return false;
        }
        else if (prev_edge.first == canonical_edge.second or prev_edge.second == canonical_edge.second){
            return true;
        }
        else{
            throw runtime_error("ERROR: no previously traversed node matches current nodes in edge");
        }
    }
}


void PileupGenerator::generate_from_bipartition() {

    BicliqueIterator biclique_iterator;
    pileup.paths.resize(bipartite_graph.left_size() + bipartite_graph.right_size());

    int i = 0;
    while (PileupGenerator::traverse_bipartition_edges(biclique_iterator)) {
        debug_print(biclique_iterator);
        std::cout << '\n';

        auto iter = overlaps.canonicalize_and_find(biclique_iterator.edge, gfa_graph);
        const edge_t& canonical_edge = iter->first;
        Alignment& alignment = iter->second;
        pair<size_t, size_t> lengths;

        alignment.compute_lengths(lengths);
        size_t left_start = gfa_graph.get_length(canonical_edge.first) - lengths.first;
        size_t right_start = 0;

        handle_t pseudo_reference;
        int64_t pseudo_ref_id;
        size_t pseudo_ref_length;
        size_t prev_pseudo_ref_length;

        handle_t pseudo_query;
        int64_t pseudo_query_id;

        bool pseudoref_reversed = pseudoref_is_reversed(canonical_edge, biclique_iterator);

        // Figure out which sequence is being treated as the "ref" in the overlap, call it pseudo reference
        if (pseudoref_reversed) {
            pseudo_reference = canonical_edge.second;
            pseudo_query = canonical_edge.first;
            pseudo_ref_length = lengths.second;
        }
        else{
            pseudo_reference = canonical_edge.first;
            pseudo_query = canonical_edge.second;
            pseudo_ref_length = lengths.first;
        }

        // Find out whether the pseudo ref has been added to the gfa_graph and what its length was
        if (not pileup.id_map.exists(pseudo_reference)){
            prev_pseudo_ref_length = 0;
        }
        else{
            pseudo_ref_id = pileup.id_map.get_id(pseudo_reference);
            prev_pseudo_ref_length = pileup.paths[pseudo_ref_id].size();
        }

        // If the previous alignment that used this pseudo ref did not span enough bases, the gfa_graph will need updating.
        // This also happens for the first alignment in the biclique, because no "ref" exists yet.
        if (prev_pseudo_ref_length < pseudo_ref_length){
            // This is the first sequence in the gfa_graph
            if (not pileup.id_map.exists(pseudo_reference)) {
                pseudo_ref_id = pileup.id_map.insert(pseudo_reference);
            }
            update_pseudoref(
                    pseudo_reference,
                    pseudo_ref_length,
                    prev_pseudo_ref_length,
                    pseudo_ref_id,
                    pseudoref_reversed);

        }

        {
            std::cout << "pseudoref id:\t" << gfa_graph.get_id(pseudo_reference) << '\n';
            std::cout << "pseudoref sequence:\t" << gfa_graph.get_sequence(pseudo_reference) << '\n';
            string test_path = "test_alignment_graph_" + std::to_string(i) + ".gfa";
            handle_graph_to_gfa(pileup.graph, test_path);
            i++;
        }

        std::cout << gfa_graph.get_length(canonical_edge.first) << " " << lengths.first << " " << left_start << " "
                  << gfa_graph.get_length(canonical_edge.second) << " " << lengths.second << " " << right_start << '\n';
        AlignmentIterator alignment_iterator(left_start, right_start);

        if (not pileup.id_map.exists(pseudo_query)) {
            pseudo_query_id = pileup.id_map.insert(pseudo_query);
        }
        size_t pseudo_ref_index = 0;
        int64_t offset = 0;
        size_t pseudo_query_index = 0;

        while (alignment.step_through_alignment(alignment_iterator)) {
            char pseudo_query_base;
            char pseudo_ref_base;
            bool is_pseudo_query_move;
            bool is_pseudo_ref_move;
            uint8_t code = alignment.operations[alignment_iterator.cigar_index].code;

            // Is the traversal walking forwards or backwards? Reassign the relevant flags/data accordingly
            if (pseudoref_reversed) {
                pseudo_query_base = gfa_graph.get_base(canonical_edge.first, alignment_iterator.ref_index);
                pseudo_ref_base = gfa_graph.get_base(canonical_edge.second, alignment_iterator.query_index);
                is_pseudo_query_move = Alignment::is_ref_move[code];
                is_pseudo_ref_move = Alignment::is_query_move[code];
            } else {
                pseudo_query_base = gfa_graph.get_base(canonical_edge.second, alignment_iterator.query_index);
                pseudo_ref_base = gfa_graph.get_base(canonical_edge.first, alignment_iterator.ref_index);
                is_pseudo_query_move = Alignment::is_query_move[code];
                is_pseudo_ref_move = Alignment::is_ref_move[code];
                offset = prev_pseudo_ref_length - pseudo_ref_length;
            }

            auto pseudo_ref_node = pileup.paths[pseudo_ref_id][pseudo_ref_index + max(int64_t(0), offset)];


            {
                std::cout << "offset:\t" << offset << '\n';
                std::cout << "pseudo_ref_index:\t" << pseudo_ref_index << '\n';

                std::cout <<
                          "query_base:\t" << pseudo_query_base << '\n' <<
                          "query_index:\t" << alignment_iterator.query_index << '\n' <<
                          "ref nid:\t" << gfa_graph.get_id(pseudo_ref_node) << '\n' <<
                          "ref_base:\t" << pseudo_ref_base << '\n' <<
                          "ref_index:\t" << alignment_iterator.ref_index << '\n' <<
                          "query_move:\t" << is_pseudo_query_move << '\n' <<
                          "ref_move:\t" << is_pseudo_ref_move << '\n' << '\n';

                std::cout << "Pseudoref node ids:\n";
                for (auto& item: pileup.paths[pseudo_ref_id]) {
                    std::cout << gfa_graph.get_id(item) << ",";
                }
                std::cout << '\n' << '\n';

                string test_path = "test_alignment_graph_" + std::to_string(i) + ".gfa";
                handle_graph_to_gfa(pileup.graph, test_path);
            }

            // Match or mismatch
            if (is_pseudo_query_move and is_pseudo_ref_move) {
                /// Match
                if (pseudo_ref_base == pseudo_query_base) {
                    std::cout << "MATCH\n";
                    // Connect it to whatever the previous node in this path was (which may already be connected)
                    if (not pileup.paths[pseudo_query_id].empty()) {
                        auto prev_node = pileup.paths[pseudo_query_id].back();

                        std::cout << "Pseudoquery node ids:\n";
                        for (auto& item: pileup.paths[pseudo_query_id]){
                            std::cout << pileup.graph.get_id(item) << ",";
                        }
                        std::cout << '\n' << '\n';

                        auto prev_id = pileup.graph.get_id(prev_node);
                        auto _pseudo_ref_id = pileup.graph.get_id(pseudo_ref_node);
                        std::cout << prev_id << " " << pileup.graph.has_node(prev_id) << '\n';
                        std::cout << _pseudo_ref_id << " " << pileup.graph.has_node(_pseudo_ref_id) << '\n';

                        if (not pileup.graph.has_edge(prev_node, pseudo_ref_node)) {
                            pileup.graph.create_edge(prev_node, pseudo_ref_node);
                            std::cout << "created edge:\t"
                                      << pileup.graph.get_id(prev_node) << "->"
                                      << pileup.graph.get_id(pseudo_ref_node) << " "
                                      << pileup.graph.has_edge(prev_node, pseudo_ref_node) << '\n';
                        }
                    }

                    // Simply update the path for this sequence to follow the ref
                    pileup.paths[pseudo_query_id].push_back(pseudo_ref_node);

                }
                /// Mismatch
                else {
                    std::cout << "MISMATCH\n";

                    // Create a new node and splice into the gfa_graph
                    auto node = pileup.graph.create_handle(string(1,pseudo_query_base));
                    std::cout << "created node:\t" << pileup.graph.get_id(node) << " " << pileup.graph.has_node(gfa_graph.get_id(node)) << '\n';

                    if (not pileup.paths[pseudo_query_id].empty()){
                        auto prev_node = pileup.paths[pseudo_query_id].back();
                        pileup.graph.create_edge(prev_node, node);
                        std::cout << "created edge:\t"
                                  << pileup.graph.get_id(prev_node) << "->"
                                  << pileup.graph.get_id(node) << " "
                                  << pileup.graph.has_edge(prev_node, node) << '\n';
                    }

                    pileup.paths[pseudo_query_id].push_back(node);
                    std::cout << "Pseudoquery node ids:\n";
                    for (auto& item: pileup.paths[pseudo_query_id]){
                        std::cout << pileup.graph.get_id(item) << ",";
                    }
                    std::cout << '\n' << '\n';

                    std::cout << "created node:\t" << pileup.graph.get_id(pileup.paths[pseudo_query_id].back()) << " " << pileup.graph.has_node(gfa_graph.get_id(pileup.paths[pseudo_query_id].back())) << '\n';
                }
                pseudo_query_index++;
                pseudo_ref_index++;
            }
            // Delete-like operation (relative to the pseudoref)
            else if (not is_pseudo_query_move and is_pseudo_ref_move) {
                std::cout << "DELETE\n";
                // Don't need to do anything in this situation
                pseudo_ref_index++;
            }
            // Insert-like operation (relative to the pseudoref)
            else if (is_pseudo_query_move and not is_pseudo_ref_move) {
                std::cout << "INSERT\n";

                // Similar to a mismatch
                // Create a new node and splice into the gfa_graph
                auto node = pileup.graph.create_handle(string(1,pseudo_query_base));
                std::cout << "created node:\t" << gfa_graph.get_id(node) << " " << gfa_graph.has_node(gfa_graph.get_id(node)) << '\n';

                if (not pileup.paths[pseudo_query_id].empty()){
                    auto prev_node = pileup.paths[pseudo_query_id].back();
                    pileup.graph.create_edge(prev_node, node);
                    std::cout << "created edge:\t"
                              << pileup.graph.get_id(prev_node) << "->"
                              << pileup.graph.get_id(node) << " "
                              << pileup.graph.has_edge(prev_node, node) << '\n';
                }

                pileup.paths[pseudo_query_id].push_back(node);

                pseudo_query_index++;
            }
            // Should never happen?
            else{

            }
        }

        pileup.edges_traversed.push_back(biclique_iterator.edge);
    }

    {
        string test_path = "test_alignment_graph_" + std::to_string(i) + ".gfa";
        handle_graph_to_gfa(pileup.graph, test_path);
    }

}


void PoaPileupGenerator::add_alignments_to_poa(Graph& spoa_graph, unique_ptr<AlignmentEngine>& alignment_engine){

    // If the graph already has some sequences in it, then start the id at that number
    uint32_t spoa_id = spoa_graph.sequences().size();

    for (bool is_left: {false, true}){
        for (uint64_t id=0; id < pileup.splice_site_indexes[!is_left].size(); id++) {
            // Assume the alignment data is reverse sorted by length at this point (not totally safe, but safe enough)
            auto data = get_splice_data(is_left,id,0);

            auto node_handle = pileup.id_map[!is_left].get_name(id);
            auto gfa_node_length = gfa_graph.get_length(node_handle);
            auto start_index = data.get_start_index(gfa_node_length);

            string subsequence = gfa_graph.get_subsequence(node_handle, start_index, data.length);

            std::cout << gfa_node_length << " " << start_index << " " << data.length << '\n';
            std::cout << subsequence << '\n';

            auto alignment = alignment_engine->Align(subsequence, spoa_graph);
            spoa_graph.AddAlignment(alignment, subsequence);

            data.spoa_id = spoa_id++;
        }
    }
}


void PoaPileupGenerator::convert_spoa_to_bdsg(Graph& spoa_graph){

    auto paths = spoa_graph.sequences();
    unordered_map <uint32_t, handle_t> nodes_created;
    handle_t previous_pileup_node;

    for (bool is_left: {false, true}){

        // Iterate all the sequences involved in this biclique (find them by their coordinates)
        for (uint64_t id=0; id < pileup.splice_site_indexes[!is_left].size(); id++) {
            // Get the longest overlap for this sequence
            auto alignment_data = get_splice_data(is_left, id, 0);

            // This points to the first SPOA node within the path that this sequence aligned to in the SPOA graph
            auto node = paths[alignment_data.spoa_id];

            auto gfa_node_length = gfa_graph.get_length(pileup.id_map[!is_left].get_name(id));
            uint64_t gfa_start_index = alignment_data.get_start_index(gfa_node_length);
            uint64_t base_index = 0;

            // Iterate through the path
            while (true) {
//              std::cout << node->code << " " << node->id << " " << node->Coverage() << " " << '\n';

                // Check if this node has already been copied to the BDSGraph
                auto iter = nodes_created.find(node->id);

                if (iter == nodes_created.end()) {
                    auto gfa_handle = pileup.id_map[!is_left].get_name(id);
                    char base = gfa_graph.get_base(gfa_handle, gfa_start_index + base_index);

                    auto new_pileup_node = pileup.graph.create_handle(string(1, base));
                    nodes_created.emplace(node->id, new_pileup_node);

                    if (base_index > 0) {
                        pileup.graph.create_edge(previous_pileup_node, new_pileup_node);
                    }

                    previous_pileup_node = new_pileup_node;
                } else {
                    if (base_index > 0) {
                        pileup.graph.create_edge(previous_pileup_node, iter->second);
                    }

                    previous_pileup_node = iter->second;
                }

                // If there are other (shorter) overlaps that this node participated in, then they also need paths
                for (auto& splice_data_index: pileup.splice_site_indexes[!is_left][id]) {
                    auto splice_data = get_splice_data(is_left, id, splice_data_index);
                    auto gfa_index = gfa_start_index + base_index;

                    auto other_start_index = splice_data.get_start_index(gfa_node_length);
                    auto other_stop_index = splice_data.get_stop_index(gfa_node_length);

                    if (gfa_index >= other_start_index and gfa_index <= other_stop_index) {
                        auto path_handle = pileup.graph.get_path_handle(splice_data.path_name);
                        pileup.graph.append_step(path_handle, previous_pileup_node);
                    }
                }

                base_index++;

                // Move to the next SPOA node
                if (!(node = node->Successor(alignment_data.spoa_id))) {
                    break;
                }
            }
        }
    }
}


//void PileupGenerator::generate_unique_path_name(){
//    auto side_string = to_string(is_left);
//    auto id_string = to_string(graph.get_id(node));
//    auto overlap_index_string = to_string(splice_site_indexes[!is_left][id].size());
//    auto biclique_index_string = to_string(biclique_index);
//    auto component_index_string = to_string(component_index);
//
//    path_name = component_index_string + '_' +
//                biclique_index_string + '_' +
//                id_string + '_' +
//                to_string(id) + '_' +
//                side_string + '_' +
//                overlap_index_string;
//
//
//}


SpliceData& PoaPileupGenerator::get_splice_data(bool is_left, size_t id, size_t i){
    auto gfa_node_handle = pileup.id_map[!is_left].get_name(id);
    auto node_id = gfa_graph.get_id(gfa_node_handle);
    auto splice_site_index = pileup.splice_site_indexes[!is_left][id][i];

    return splice_sites[node_id][splice_site_index];
}


void PoaPileupGenerator::sort_alignment_data_by_length(){
    for (bool is_left: {false, true}){
        for (size_t id=0; id<pileup.splice_site_indexes[is_left].size(); id++) {
            auto& indexes = pileup.splice_site_indexes[is_left][id];

            sort(indexes.begin(), indexes.end(), [&](const size_t& a, const size_t& b){
                auto a_value = get_splice_data(is_left, id, a).length;
                auto b_value = get_splice_data(is_left, id, b).length;

                return a_value > b_value;
            });
        }
    }
}


void PoaPileupGenerator::generate_from_edges(const vector<edge_t>& edges, size_t component_index) {
    
    for (auto non_canonical_edge : edges) {
        auto iter = overlaps.canonicalize_and_find(non_canonical_edge, gfa_graph);
        
        const edge_t& edge = iter->first;
        Alignment& alignment = iter->second;
        
        // Do a quick check to see if this is actually a non-overlap (0M). If so, don't process it
        if (alignment.operations.empty() or (alignment.operations[0].length == 0)){
            pileup.blunt_edges.push_back(edge);
            continue;
        }
        
        pair<size_t, size_t> lengths;

        std::cout << id_map.get_name(gfa_graph.get_id(edge.first)) << "->";
        std::cout << id_map.get_name(gfa_graph.get_id(edge.second)) << '\n';
        
        iter->second.compute_lengths(lengths);

        {
            size_t start = gfa_graph.get_length(edge.first) - lengths.first;
            std::cout << lengths.first << " " << lengths.second << '\n';
            std::cout << iter->second.create_formatted_alignment_string(gfa_graph, edge, start, 0) << '\n';
            std::cout << '\n';
        }

        alignment.compute_lengths(lengths);

        handle_t reference = edge.first;
        int64_t ref_id;
        size_t ref_overlap_length;
        bool ref_is_reverse;
        ref_overlap_length = lengths.first;
        ref_is_reverse = gfa_graph.get_is_reverse(reference);
        ref_id = gfa_graph.get_id(reference);

        handle_t query = edge.second;
        int64_t query_id;
        size_t query_overlap_length;
        bool query_is_reverse;
        query_overlap_length = lengths.second;
        query_is_reverse = gfa_graph.get_is_reverse(query);
        query_id = gfa_graph.get_id(query);

        std::cout << "REF IS_REVERSE: " << ref_is_reverse << '\n';
        std::cout << "QUERY IS_REVERSE: " << query_is_reverse << '\n';

        splice_site_mutexes[ref_id].lock();
        size_t ref_splice_site_index = splice_sites[ref_id].size();
        string ref_path_name = to_string(ref_id) + '_' + to_string(ref_splice_site_index);
        pileup.graph.create_path_handle(ref_path_name);

        pileup.update_alignment_data(true, reference, ref_splice_site_index);

        splice_sites[ref_id].emplace_back(ref_is_reverse,
                                          true,
                                          ref_overlap_length,
                                          ref_path_name,
                                          pileup.biclique_index,
                                          component_index);

        splice_site_mutexes[ref_id].unlock();

        splice_site_mutexes[query_id].lock();
        size_t query_splice_site_index = splice_sites[query_id].size();
        string query_path_name = to_string(query_id) + to_string(query_splice_site_index);
        pileup.graph.create_path_handle(query_path_name);

        pileup.update_alignment_data(false, query, query_splice_site_index);

        splice_sites[query_id].emplace_back(query_is_reverse,
                                            false,
                                            query_overlap_length,
                                            query_path_name,
                                            pileup.biclique_index,
                                            component_index);

        splice_site_mutexes[query_id].unlock();

//        {
//            size_t ref_length = gfa_graph.get_length(reference);
//            size_t query_length = gfa_graph.get_length(query);
//            std::cout << "ref id:\t" << gfa_graph.get_id(reference) << '\n';
//            std::cout << "ref sequence:\t" << gfa_graph.get_sequence(reference) << '\n';
//            std::cout << "query sequence:\t" << gfa_graph.get_sequence(query) << '\n';
//            std::cout <<
//                      "ref_overlap_length\t" << ref_overlap_length << '\n' <<
//                      "query_overlap_length\t" << query_overlap_length << '\n' <<
//                      "ref_start\t" << splice_sites[ref_id].back().get_start_index(ref_length) << '\n' <<
//                      "ref_stop\t" << splice_sites[ref_id].back().get_stop_index(ref_length) << '\n' <<
//                      "query_start\t" << splice_sites[query_id].back().get_start_index(query_length) << '\n' <<
//                      "query_stop\t" << splice_sites[query_id].back().get_stop_index(query_length) << '\n' << '\n';
//        }
    }
    
    sort_alignment_data_by_length();
    
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 5, -3, -3, -1);

    spoa::Graph spoa_graph{};
    
    add_alignments_to_poa(spoa_graph, alignment_engine);

    auto consensus = spoa_graph.GenerateConsensus();
    
//    std::cout << ">Consensus: " << consensus << std::endl;
    
    spoa::Graph seeded_spoa_graph{};
    
    auto seeded_alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 6, -2, -4, -1);
    
    auto alignment = seeded_alignment_engine->Align(consensus, seeded_spoa_graph);
    seeded_spoa_graph.AddAlignment(alignment, consensus);
    
    // Iterate a second time on alignment, this time with consensus as the seed
    add_alignments_to_poa(seeded_spoa_graph, seeded_alignment_engine);

    {
        auto seeded_consensus = seeded_spoa_graph.GenerateConsensus();

        seeded_spoa_graph.PrintDot("spoa_overlap.dot");
        auto msa = seeded_spoa_graph.GenerateMultipleSequenceAlignment();

        for (const auto& it : msa) {
            std::cout << it << std::endl;
        }
        std::cout << '\n';
    }

    convert_spoa_to_bdsg(seeded_spoa_graph);

    handle_graph_to_gfa(pileup.graph, "test_output.gfa");

    unchop(&pileup.graph);
    
    handle_graph_to_gfa(pileup.graph, "test_output_unchopped.gfa");

//    {
//        pileup.gfa_graph.for_each_path_handle([&](const path_handle_t& p) {
//            pileup.gfa_graph.for_each_step_in_path(p, [&](const step_handle_t& step) {
//                auto h = pileup.gfa_graph.get_handle_of_step(step);
//                std::cout << pileup.gfa_graph.get_sequence(h);
//            });
//            std::cout << "\n";
//        });
//    }
}


void PoaPileupGenerator::generate_from_bipartition(const bipartition& bipartition, size_t component_index) {

    vector<edge_t> edges;
    edges.reserve(bipartition.first.size() * bipartition.second.size());
    for (const auto& left: bipartition.first) {
        for (const auto& right: bipartition.second) {
            edges.emplace_back(left, gfa_graph.flip(right));
        }
    }
    
    generate_from_edges(edges, component_index);
}



}
