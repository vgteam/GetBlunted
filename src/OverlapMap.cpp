#include "OverlapMap.hpp"

using std::make_pair;
using std::to_string;


namespace bluntifier{

OverlapMap::OverlapMap()=default;


Alignment OverlapMap::insert(const gfak::edge_elem& e, handle_t source, handle_t sink){
    Alignment alignment(e.alignment);
    overlaps.emplace(make_pair(source, sink), alignment);

    return alignment;
}


void OverlapMap::insert(Alignment& alignment, handle_t source, handle_t sink){
    overlaps.emplace(make_pair(source, sink), alignment);
}


Alignment OverlapMap::insert(const gfak::edge_elem& e, const edge_t& edge_handle){
    Alignment alignment(e.alignment);
    overlaps.emplace(edge_handle, alignment);

    return alignment;
}


void OverlapMap::update_edge(const edge_t& old_edge, const edge_t& new_edge){
    Alignment alignment = overlaps.find(old_edge)->second;

    overlaps.erase(old_edge);
    overlaps.emplace(new_edge, alignment);
}


unordered_map<edge_t,Alignment>::iterator OverlapMap::at(const handle_t source, handle_t sink){
    return overlaps.find(make_pair(source, sink));
}


unordered_map<edge_t,Alignment>::iterator OverlapMap::at(const edge_t& edge_handle){
    return overlaps.find(edge_handle);
}


unordered_map<edge_t,Alignment>::const_iterator OverlapMap::at(const handle_t source, handle_t sink)const{
    return overlaps.find(make_pair(source, sink));
}


unordered_map<edge_t,Alignment>::const_iterator OverlapMap::at(const edge_t& edge_handle)const{
    return overlaps.find(edge_handle);
}


unordered_map<edge_t,Alignment>::iterator OverlapMap::canonicalize_and_find(const edge_t& edge, const HandleGraph& graph){
    // Check if edge is found in the overlaps in its provided orientation
    auto iter = overlaps.find(edge);
    if (iter == overlaps.end()){
        // If not found, then flip it and search again
        handle_t right = graph.flip(edge.first);
        handle_t left = graph.flip(edge.second);

        iter = overlaps.find({left, right});

        // If it still isn't found then throw an error
        if (iter == overlaps.end()) {
            std::cerr << "ERROR: edge not found in overlaps:\n\t("
                                << graph.get_id(edge.first) << '-'
                                << graph.get_is_reverse(edge.first) << ")->("
                                << graph.get_id(edge.second) << '-'
                                << graph.get_is_reverse(edge.second) << ')';
        }
    }

    return iter;
}


unordered_map<edge_t,Alignment>::const_iterator OverlapMap::canonicalize_and_find(
        const edge_t& edge,
        const HandleGraph& graph) const{

    // Check if edge is found in the overlaps in its provided orientation
    auto iter = overlaps.find(edge);
    if (iter == overlaps.end()){
        // If not found, then flip it and search again
        handle_t right = graph.flip(edge.first);
        handle_t left = graph.flip(edge.second);

        iter = overlaps.find({left,right});

        // If it still isn't found then throw an error
        if (iter == overlaps.end()) {
            throw runtime_error("ERROR: edge not found in overlaps:\n\t("
                                + to_string(graph.get_id(edge.first)) + '-'
                                + to_string(graph.get_is_reverse(edge.first)) + ")->("
                                + to_string(graph.get_id(edge.second)) + '-'
                                + to_string(graph.get_is_reverse(edge.second)) + ')');
        }
    }

    return iter;
}


void OverlapMap::canonicalize_and_compute_lengths(pair<size_t,size_t>& lengths, edge_t& edge, const HandleGraph& graph){
    auto iter = canonicalize_and_find(edge, graph);
    iter->second.compute_lengths(lengths);
}





}
