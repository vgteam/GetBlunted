#include "PileupGenerator.hpp"


namespace bluntifier {


void BicliqueIterator::update(edge_t edge){
    this->visited.insert(edge.first);
    this->visited.insert(edge.second);
    this->edge = edge;
}


bool PileupGenerator::traverse_bipartition(
        HandleGraph& graph,
        OverlapMap& overlaps,
        bipartition& partition,
        BicliqueIterator& iterator){

    bool done = false;

    if (iterator.first_step){
        // Pick an arbitrary node to start with
        handle_t start;

        if (partition.first.size() > partition.second.size()) {
            start = *partition.first.begin();
        }
        else{
            start = graph.flip(*partition.second.begin());
        }

        iterator.nodes.push(start);
        iterator.visited.insert(start);
    }

    // Emulate DFS
    handle_t node = iterator.nodes.top();
    if (iterator.visited.count(node) == 0){
        iterator.visited.insert(node);

        // Just search the other side of the partition
        if (partition.first.count(node) > 0){
            for (auto& neighbor: partition.second){
                iterator.nodes.push(neighbor);
            }
        }
        if (partition.second.count(node) > 0){
            for (auto& neighbor: partition.first){
                iterator.nodes.push(neighbor);
            }
        }
    }


    if (iterator.visited.size() == partition.first.size() + partition.second.size()) {
        done = true;
    }

    return not done;
}


void PileupGenerator::generate_from_bipartition(
        bipartition& partition,
        OverlapMap& overlaps,
        HandleGraph& graph,
        Pileup& pileup) {

}



}