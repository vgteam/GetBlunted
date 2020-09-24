#include "PileupGenerator.hpp"
#include "handle_to_gfa.hpp"

using std::max;

namespace bluntifier {


BicliqueIterator::BicliqueIterator():
    first_step(true)
{}


PileupGenerator::PileupGenerator()=default;


bool PileupGenerator::traverse_bipartition_edges(
        const HandleGraph& graph,
        const OverlapMap& overlaps,
        const BipartiteGraph& bipartite_graph,
        BicliqueIterator& iterator) {

    bool done;

    if (iterator.first_step) {
        PileupGenerator::traverse_bipartition_nodes(graph, overlaps, bipartite_graph, iterator);
        if (iterator.is_left) {
            iterator.prev_left_node = iterator.node;
            iterator.prev_left_is_valid = true;
        } else {
            iterator.prev_right_node = iterator.node;
            iterator.prev_right_is_valid = true;
        }
    }

    done = not PileupGenerator::traverse_bipartition_nodes(graph, overlaps, bipartite_graph, iterator);

    if (not done) {
        if (iterator.is_left) {
            if (iterator.prev_right_is_valid) {
                iterator.edge = {iterator.node, graph.flip(iterator.prev_right_node)};
            }
            iterator.prev_left_node = iterator.node;
            iterator.prev_left_is_valid = true;
        } else {
            if (iterator.prev_left_is_valid) {
                iterator.edge = {iterator.prev_left_node, graph.flip(iterator.node)};
            }
            iterator.prev_right_node = iterator.node;
            iterator.prev_right_is_valid = true;
        }
    }

    return not done;
}


bool PileupGenerator::traverse_bipartition_nodes(
        const HandleGraph& graph,
        const OverlapMap& overlaps,
        const BipartiteGraph& bipartite_graph,
        BicliqueIterator& iterator){

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


void PileupGenerator::debug_print(const IncrementalIdMap<string>& id_map,
                                  HandleGraph& graph,
                                  OverlapMap& overlaps,
                                  BicliqueIterator& biclique_iterator){

    auto iter = overlaps.canonicalize_and_find(biclique_iterator.edge, graph);
    const edge_t& edge = iter->first;
    Alignment& alignment = iter->second;
    pair<size_t, size_t> lengths;
    size_t start;

    std::cout << id_map.get_name(graph.get_id(edge.first)) << "->";
    std::cout << id_map.get_name(graph.get_id(edge.second)) << '\n';

    iter->second.compute_lengths(lengths);

    start = graph.get_length(edge.first) - lengths.first;

    std::cout << "Is left: " << biclique_iterator.is_left << '\n';
    std::cout << lengths.first << " " << lengths.second << '\n';
    std::cout << iter->second.create_formatted_alignment_string(graph, edge, start, 0) << '\n';
}


void PileupGenerator::update_pseudoref(
        Pileup& pileup,
        HandleGraph& graph,
        handle_t pseudo_reference,
        size_t pseudo_ref_length,
        size_t prev_pseudo_ref_length,
        int64_t pseudo_ref_id,
        bool is_left) {

    // Walk from the last index of the prev pseudoref up to the current index
    if (is_left) {
        for (size_t i = prev_pseudo_ref_length; i < pseudo_ref_length; i++) {
            auto h = pileup.graph.create_handle(string(1, graph.get_base(pseudo_reference, i)));
            std::cout << "created node:\t" << graph.get_id(h) << " " << graph.has_node(graph.get_id(h)) << '\n';

            if (not pileup.paths[pseudo_ref_id].empty()) {
                pileup.graph.create_edge(pileup.paths[pseudo_ref_id].back(), h);
                std::cout << "created edge:\t"
                          << pileup.graph.get_id(pileup.paths[pseudo_ref_id].back()) << "->"
                          << pileup.graph.get_id(h) << " "
                          << pileup.graph.has_edge(pileup.paths[pseudo_ref_id].back(), h) << '\n';
            }
            pileup.paths[pseudo_ref_id].push_back(h);

            std::cout << "base:\t" << graph.get_base(pseudo_reference, i) << '\n';
            std::cout << "id:\t" << graph.get_id(h) << '\n' << '\n';
        }
    }
    // Walk BACKWARDS from the last index of the prev pseudoref down to the current index
    else {
        size_t start = graph.get_length(pseudo_reference) - prev_pseudo_ref_length - 1;
        size_t stop = graph.get_length(pseudo_reference) - pseudo_ref_length - 1;

        std::cout <<  graph.get_length(pseudo_reference) << " " << pileup.paths[pseudo_ref_id].size() << " " << start << " " << stop << '\n';

        for (int64_t i = start; i > stop; i--) {
            auto h = pileup.graph.create_handle(string(1, graph.get_base(pseudo_reference, i)));
            std::cout << "created node:\t" << graph.get_id(h) << " " << graph.has_node(graph.get_id(h)) << '\n';

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
        const Pileup& pileup,
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


void PileupGenerator::generate_from_bipartition(
        const BipartiteGraph& bipartite_graph,
        const IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps,
        HandleGraph& graph,
        Pileup& pileup) {

    BicliqueIterator biclique_iterator;
    pileup.paths.resize(bipartite_graph.left_size() + bipartite_graph.right_size());

    int i = 0;
    while (PileupGenerator::traverse_bipartition_edges(graph, overlaps, bipartite_graph, biclique_iterator)) {
        debug_print(id_map, graph, overlaps, biclique_iterator);
        std::cout << '\n';

        auto iter = overlaps.canonicalize_and_find(biclique_iterator.edge, graph);
        const edge_t& canonical_edge = iter->first;
        Alignment& alignment = iter->second;
        pair<size_t, size_t> lengths;

        alignment.compute_lengths(lengths);
        size_t left_start = graph.get_length(canonical_edge.first) - lengths.first;
        size_t right_start = 0;

        handle_t pseudo_reference;
        int64_t pseudo_ref_id;
        size_t pseudo_ref_length;
        size_t prev_pseudo_ref_length;

        handle_t pseudo_query;
        int64_t pseudo_query_id;

        bool pseudoref_is_reversed = PileupGenerator::pseudoref_is_reversed(pileup, canonical_edge, biclique_iterator);

        // Figure out which sequence is being treated as the "ref" in the overlap, call it pseudo reference
        if (pseudoref_is_reversed) {
            pseudo_reference = canonical_edge.second;
            pseudo_query = canonical_edge.first;
            pseudo_ref_length = lengths.second;
        }
        else{
            pseudo_reference = canonical_edge.first;
            pseudo_query = canonical_edge.second;
            pseudo_ref_length = lengths.first;
        }

        // Find out whether the pseudo ref has been added to the graph and what its length was
        if (not pileup.id_map.exists(pseudo_reference)){
            prev_pseudo_ref_length = 0;
        }
        else{
            pseudo_ref_id = pileup.id_map.get_id(pseudo_reference);
            prev_pseudo_ref_length = pileup.paths[pseudo_ref_id].size();
        }

        // If the previous alignment that used this pseudo ref did not span enough bases, the graph will need updating.
        // This also happens for the first alignment in the biclique, because no "ref" exists yet.
        if (prev_pseudo_ref_length < pseudo_ref_length){
            // This is the first sequence in the graph
            if (not pileup.id_map.exists(pseudo_reference)) {
                pseudo_ref_id = pileup.id_map.insert(pseudo_reference);
            }
            update_pseudoref(
                    pileup,
                    graph,
                    pseudo_reference,
                    pseudo_ref_length,
                    prev_pseudo_ref_length,
                    pseudo_ref_id,
                    pseudoref_is_reversed);

        }

        {
            std::cout << "pseudoref id:\t" << graph.get_id(pseudo_reference) << '\n';
            std::cout << "pseudoref sequence:\t" << graph.get_sequence(pseudo_reference) << '\n';
            string test_path = "test_alignment_graph_" + std::to_string(i) + ".gfa";
            handle_graph_to_gfa(pileup.graph, test_path);
            i++;
        }

        std::cout << graph.get_length(canonical_edge.first) << " " <<  lengths.first << " " << left_start << " "
                  << graph.get_length(canonical_edge.second) << " " <<  lengths.second << " " << right_start << '\n';
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
            if (pseudoref_is_reversed) {
                pseudo_query_base = graph.get_base(canonical_edge.first, alignment_iterator.ref_index);
                pseudo_ref_base = graph.get_base(canonical_edge.second, alignment_iterator.query_index);
                is_pseudo_query_move = Alignment::is_ref_move[code];
                is_pseudo_ref_move = Alignment::is_query_move[code];
            } else {
                pseudo_query_base = graph.get_base(canonical_edge.second, alignment_iterator.query_index);
                pseudo_ref_base = graph.get_base(canonical_edge.first, alignment_iterator.ref_index);
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
                          "ref nid:\t" << graph.get_id(pseudo_ref_node) << '\n' <<
                          "ref_base:\t" << pseudo_ref_base << '\n' <<
                          "ref_index:\t" << alignment_iterator.ref_index << '\n' <<
                          "query_move:\t" << is_pseudo_query_move << '\n' <<
                          "ref_move:\t" << is_pseudo_ref_move << '\n' << '\n';

                std::cout << "Pseudoref node ids:\n";
                for (auto& item: pileup.paths[pseudo_ref_id]) {
                    std::cout << graph.get_id(item) << ",";
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
                        auto pseudo_ref_id = pileup.graph.get_id(pseudo_ref_node);
                        std::cout << prev_id << " " << pileup.graph.has_node(prev_id) << '\n';
                        std::cout << pseudo_ref_id << " " << pileup.graph.has_node(pseudo_ref_id) << '\n';

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

                    // Create a new node and splice into the graph
                    auto node = pileup.graph.create_handle(string(1,pseudo_query_base));
                    std::cout << "created node:\t" << pileup.graph.get_id(node) << " " << pileup.graph.has_node(graph.get_id(node)) << '\n';

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

                    std::cout << "created node:\t" << pileup.graph.get_id(pileup.paths[pseudo_query_id].back()) << " " << pileup.graph.has_node(graph.get_id(pileup.paths[pseudo_query_id].back())) << '\n';
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
                // Create a new node and splice into the graph
                auto node = pileup.graph.create_handle(string(1,pseudo_query_base));
                std::cout << "created node:\t" << graph.get_id(node) << " " << graph.has_node(graph.get_id(node)) << '\n';

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


void PileupGenerator::generate_spoa_graph_from_bipartition(
        const BipartiteGraph& bipartite_graph,
        const IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps,
        HandleGraph& graph,
        Pileup& pileup) {

    BicliqueIterator biclique_iterator;

    // For each input sequence, what are the points at which the alignment starts/ends (may be multiple)
//    vector <vector <uint64_t> > subsequence_lengths;
//    vector <vector <uint64_t> > subsequence_starts;
//    vector <vector <uint64_t> > subsequence_stops;
//    sequence_splice_positions.resize(bipartite_graph.left_size() + bipartite_graph.right_size());
//    sequence_splice_positions.resize(bipartite_graph.left_size() + bipartite_graph.right_size());
//    sequence_splice_positions.resize(bipartite_graph.left_size() + bipartite_graph.right_size());

    uint16_t i=0;
    while (PileupGenerator::traverse_bipartition_edges(graph, overlaps, bipartite_graph, biclique_iterator)) {
        debug_print(id_map, graph, overlaps, biclique_iterator);
        std::cout << '\n';

        auto iter = overlaps.canonicalize_and_find(biclique_iterator.edge, graph);
        const edge_t& canonical_edge = iter->first;
        Alignment& alignment = iter->second;
        pair<size_t, size_t> lengths;

        alignment.compute_lengths(lengths);
        size_t left_start = graph.get_length(canonical_edge.first) - lengths.first;
        size_t right_start = 0;

        handle_t pseudo_reference;
        int64_t pseudo_ref_id;
        size_t pseudo_ref_length;
        size_t pseudo_ref_start;
        size_t pseudo_ref_stop;

        handle_t pseudo_query;
        int64_t pseudo_query_id;
        size_t pseudo_query_length;
        size_t pseudo_query_start;
        size_t pseudo_query_stop;

        if (not pileup.id_map.exists(pseudo_reference)) {
            pseudo_ref_id = pileup.id_map.insert(pseudo_reference);
        }
        else{
            pseudo_ref_id = pileup.id_map.get_id(pseudo_reference);
        }

        if (not pileup.id_map.exists(pseudo_query)) {
            pseudo_query_id = pileup.id_map.insert(pseudo_query);
        }
        else{
            pseudo_query_id = pileup.id_map.get_id(pseudo_query);
        }

        bool pseudoref_is_reversed = PileupGenerator::pseudoref_is_reversed(pileup, canonical_edge, biclique_iterator);

        // Figure out which sequence is being treated as the "query" in the overlap
        if (pseudoref_is_reversed) {
            pseudo_reference = canonical_edge.second;
            pseudo_query = canonical_edge.first;

            pseudo_ref_length = lengths.second;
            pseudo_query_length = lengths.first;

            pseudo_ref_start = right_start;
            pseudo_ref_stop = pseudo_ref_length - 1;

            pseudo_query_start = left_start;
            pseudo_query_stop = graph.get_length(pseudo_query) - 1;
        }
        else{
            pseudo_reference = canonical_edge.first;
            pseudo_query = canonical_edge.second;

            pseudo_ref_length = lengths.first;
            pseudo_query_length = lengths.second;

            pseudo_ref_start = left_start;
            pseudo_ref_stop = graph.get_length(pseudo_reference) - 1;

            pseudo_query_start = right_start;
            pseudo_query_stop = pseudo_query_length - 1;
        }


        {
            std::cout << "pseudoref id:\t" << graph.get_id(pseudo_reference) << '\n';
            std::cout << "pseudoref sequence:\t" << graph.get_sequence(pseudo_reference) << '\n';
            std::cout << "pseudoquery sequence:\t" << graph.get_sequence(pseudo_query) << '\n';
            std::cout <<
                "ref_length\t" << pseudo_ref_length << '\n' <<
                "query_length\t" << pseudo_query_length << '\n' <<
                "ref_start\t" << pseudo_ref_start << '\n' <<
                "ref_stop\t" << pseudo_ref_stop << '\n' <<
                "query_start\t" << pseudo_query_start << '\n' <<
                "query_stop\t" << pseudo_query_stop << '\n' << '\n';
        }


//        std::cout << pseudo_ref_id << " " << sequence_splice_positions.size();
//        if (i == 0){
//            sequence_splice_positions[pseudo_ref_id].push_back({pseudo_ref_start, pseudo_ref_stop});
//        }
//        sequence_splice_positions[pseudo_query_id].push_back({pseudo_query_start, pseudo_query_stop});

        pileup.edges_traversed.push_back(biclique_iterator.edge);
        i++;
    }
}



}

/**

        pair<size_t, size_t> lengths;
        size_t start;

        std::cout << id_map.get_name(graph.get_id(edge.first)) << "->";
        std::cout << id_map.get_name(graph.get_id(edge.second)) << '\n';

        iter->second.compute_lengths(lengths);

        start = graph.get_length(edge.first) - lengths.first;

        std::cout << "Is left: " << biclique_iterator.is_left << '\n';
        std::cout << lengths.first << " " << lengths.second << '\n';
        std::cout << iter->second.create_formatted_alignment_string(graph, edge, start, 0) << '\n';
**/