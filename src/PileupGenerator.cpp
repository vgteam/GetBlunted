#include "PileupGenerator.hpp"


namespace bluntifier {


BicliqueIterator::BicliqueIterator():
    first_step(true)
{}


//void BicliqueIterator::update(edge_t edge){
//    this->visited.insert(edge.first);
//    this->visited.insert(edge.second);
//    this->edge = edge;
//}


PileupGenerator::PileupGenerator()=default;


bool PileupGenerator::traverse_bipartition(
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
            start = graph.flip(*bipartite_graph.right_begin());
        }

        iterator.nodes.push(start);
        iterator.first_step = false;
    }

    while (iterator.visited.size() < bipartite_graph.left_size() + bipartite_graph.right_size()) {
        // Emulate DFS
        handle_t node = iterator.nodes.top();
        iterator.nodes.pop();
        if (iterator.visited.count(node) == 0) {
            iterator.node = node;
            iterator.visited.insert(node);

            // Just search the other side of the partition
            if (bipartite_graph.is_left_side(node)) {
                for (auto it = bipartite_graph.right_begin(); it != bipartite_graph.right_end(); ++it) {
                    iterator.nodes.push(*it);
                }
            } else {
                for (auto it = bipartite_graph.left_begin(); it != bipartite_graph.left_end(); ++it) {
                    iterator.nodes.push(*it);
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