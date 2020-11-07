#include "OverlappingOverlap.hpp"



namespace bluntifier{


OverlappingChild::OverlappingChild(handle_t handle, size_t biclique_index):
        handle(handle),
        biclique_index(biclique_index)
{}


void OverlappingChild::print(HandleGraph& gfa_graph) const{
    std::cout << gfa_graph.get_id(handle) << " Biclique=" << biclique_index << " Length=" << gfa_graph.get_length(handle) << '\n';
}


OverlappingNodeInfo::OverlappingNodeInfo(nid_t parent_node):
        parent_node(parent_node)
{}

}