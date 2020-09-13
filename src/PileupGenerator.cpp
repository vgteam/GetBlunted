#include "PileupGenerator.hpp"


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

        // Choose one from the bigger set, if possible
        if (bipartite_graph.left_size() > bipartite_graph.right_size()) {
            start = *bipartite_graph.left_begin();
        }
        else{
            start = *bipartite_graph.right_begin();
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


void PileupGenerator::generate_from_bipartition(
        const BipartiteGraph& bipartite_graph,
        const IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps,
        HandleGraph& graph,
        Pileup& pileup) {

    BicliqueIterator biclique_iterator;
    bool first_step = true;
    pileup.paths.resize(bipartite_graph.left_size() + bipartite_graph.right_size());

    while (PileupGenerator::traverse_bipartition_edges(graph, overlaps, bipartite_graph, biclique_iterator)) {
        debug_print(id_map, graph, overlaps, biclique_iterator);
        std::cout << '\n';

        auto iter = overlaps.canonicalize_and_find(biclique_iterator.edge, graph);
        const edge_t& edge = iter->first;
        Alignment& alignment = iter->second;
        pair<size_t, size_t> lengths;

        alignment.compute_lengths(lengths);
        size_t left_start = graph.get_length(edge.first) - lengths.first;
        size_t right_start = 0;

        handle_t pseudo_reference;
        int64_t pseudo_ref_id;
        size_t pseudo_ref_length;
        size_t prev_pseudo_ref_length;

        // Figure out which sequence is being treated as the "ref" in the overlap, call it pseudo_ref
        if (biclique_iterator.is_left) {
            pseudo_reference = edge.second;
            pseudo_ref_length = lengths.second;
        }
        else{
            pseudo_reference = edge.first;
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

        // If the previous alignment that used this pseudo ref did not span enough bases, the graph will need updating
        if (prev_pseudo_ref_length < pseudo_ref_length){
            // Seq of interest ("pseudoref") is on the RIGHT of the overlap
            if (biclique_iterator.is_left){
                // If this pseudoref is not in the pileup graph at all, it must be the first sequence
                // In this case, just initialize the graph with it.
                if (not pileup.id_map.exists(pseudo_reference)){
                    pseudo_ref_id = pileup.id_map.insert(pseudo_reference);

                    for (size_t i=0; i<pseudo_ref_length; i++){
                        auto h = pileup.graph.create_handle(string(1,graph.get_base(pseudo_reference,i)));
                        if (i > 0){
                            pileup.graph.create_edge(pileup.paths[pseudo_ref_id][i-1], h);
                        }
                        pileup.paths[pseudo_ref_id].push_back(h);
                    }
                }
                // If the current pseudoref is longer than the previous one, then stitch the correct side
                else{
                    for (size_t i=prev_pseudo_ref_length; i<pseudo_ref_length; i++){
                        auto h = pileup.graph.create_handle(string(1,graph.get_base(pseudo_reference,i)));
                        pileup.graph.create_edge(pileup.paths[pseudo_ref_id][i-1], h);
                        pileup.paths[pseudo_ref_id].push_back(h);
                    }
                }
            }
            // Seq of interest ("pseudoref") is on the LEFT of the overlap
            else{
                // If this pseudoref is not in the pileup graph at all, it must be the first sequence
                // In this case, just initialize the graph with it.
                if (not pileup.id_map.exists(pseudo_reference)){
                    pseudo_ref_id = pileup.id_map.insert(pseudo_reference);

                    for (size_t i=0; i<pseudo_ref_length; i++){
                        auto h = pileup.graph.create_handle(string(1,graph.get_base(pseudo_reference,i)));
                        if (i > 0){
                            pileup.graph.create_edge(pileup.paths[pseudo_ref_id][i-1], h);
                        }
                        pileup.paths[pseudo_ref_id].push_back(h);
                    }
                }
                // If the current pseudoref is longer than the previous one, then stitch the correct side
                else{
//                    size_t start = graph.get_length(edge.first) - prev_pseudo_ref_length - 1;
//                    size_t stop = graph.get_length(edge.first) - pseudo_ref_length;
//                    for (int64_t i=start; i>=stop; i--){
//                        auto h = pileup.graph.create_handle(string(1,graph.get_base(pseudo_reference,i)));
//                        pileup.graph.create_edge(h, pileup.paths[pseudo_ref_id][i+1]);
//                        pileup.paths[pseudo_ref_id].push_back(h);
                    }
                }
            }
        }


        AlignmentIterator alignment_iterator(left_start, right_start);
        size_t index = 0;

        while (alignment.step_through_alignment(alignment_iterator)){
            char base;
            bool is_move;
            uint8_t code = alignment.operations[alignment_iterator.cigar_index].code;

            // Is the traversal walking backwards? If so pick the sequence accordingly
            if (biclique_iterator.is_left) {
                base = graph.get_base(edge.first, alignment_iterator.query_index);
                is_move = Alignment::is_query_move[code];
            }else{
                base = graph.get_base(edge.second, alignment_iterator.ref_index);
                is_move = Alignment::is_ref_move[code];
            }

            if (is_move) {
            }

            index++;
        }
        first_step = false;

        string s;
        pileup.to_string(s);
        std::cout << s << '\n';
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