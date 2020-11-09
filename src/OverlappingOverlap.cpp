#include "OverlappingOverlap.hpp"

using std::cout;

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


void OverlappingNodeInfo::print(HandleGraph& graph){
    for (auto side: {0,1}) {
        cout << "Side: " << side << '\n';
        cout << "Overlapping:\n";
        for (auto& item1: overlapping_children[side]){
            auto& overlapping_child = item1.second;
            cout << '\t';
            overlapping_child.print(graph);
        }
        cout << "Normal:\n";
        for (auto& item1: normal_children[side]){
            auto& normal_child = item1.second;
            cout << '\t';
            normal_child.print(graph);
        }
        cout << '\n';
    }
    cout << '\n';

}


void OverlappingNodeInfo::find_overlapping_children_by_base_index(size_t index, bool side, deque<OverlappingChild>& children){
    if (index == 0){
        for (auto& item: overlapping_children[side]){
            children.emplace_back(item.second);
        }
    }
    else {
        while (true){
            auto result = overlapping_children[side].upper_bound(index - 1);

            if (result != overlapping_children[side].end()){
                children.emplace_back(result->second);
                index = result->first + 1;
            }
            else{
                break;
            }
        }
    }
}


void OverlappingNodeInfo::find_normal_children_by_base_index(size_t index, bool side, deque<OverlappingChild>& children){
    if (index == 0){
        for (auto& item: normal_children[side]){
            children.emplace_back(item.second);
        }
    }
    else {
        while (true){
            auto result = normal_children[side].upper_bound(index - 1);

            if (result != normal_children[side].end()){
                children.emplace_back(result->second);
                index = result->first + 1;
            }
            else{
                break;
            }
        }
    }
}



}