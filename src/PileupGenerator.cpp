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


//void PileupGenerator::generate_from_bipartition(
//        bipartition& partition,
//        OverlapMap& overlaps,
//        HandleGraph& graph,
//        Pileup& pileup) {
//
//}



}