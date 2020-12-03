#include "OverlappingOverlap.hpp"

using std::cerr;

namespace bluntifier{


OverlappingChild::OverlappingChild(handle_t handle, size_t biclique_index, bool side):
        handle(handle),
        biclique_index(biclique_index),
        side(side)
{}


OverlappingChild::OverlappingChild():
        handle(),
        biclique_index(0)
{}

void OverlappingChild::print(HandleGraph& gfa_graph) const{
    std::cerr << "ID=" << gfa_graph.get_id(handle) << " Biclique=" << biclique_index << " Length=" << gfa_graph.get_length(handle) << '\n';
}


OverlappingNodeInfo::OverlappingNodeInfo(nid_t parent_node):
        parent_node(parent_node)
{}


void OverlappingNodeInfo::print(HandleGraph& graph){
    for (auto side: {0,1}) {
        cerr << "Side: " << side << '\n';
        cerr << "Overlapping:\n";
        for (auto& item1: overlapping_children[side]){
            auto& overlapping_child = item1.second;
            cerr << '\t' << item1.first << " ";
            overlapping_child.print(graph);
        }
        cerr << "Normal:\n";
        for (auto& item1: normal_children[side]){
            auto& normal_child = item1.second;
            cerr << '\t' << item1.first << " ";
            normal_child.print(graph);
        }
        cerr << '\n';
    }
    cerr << '\n';

}



}