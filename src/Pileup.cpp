#include "Pileup.hpp"
#include <algorithm>



using std::sort;


namespace bluntifier{


SpliceData::SpliceData(uint64_t start, uint64_t stop, string& path_name):
        sequence_start_index(start),
        sequence_stop_index(stop),
        path_name(path_name)
{}


SpliceData::SpliceData(
        bool is_reverse,
        bool is_left,
        uint64_t start,
        uint64_t stop,
        string& path_name,
        size_t pileup_index,
        size_t component_index):
        is_reverse(is_reverse),
        is_left(is_left),
        sequence_start_index(start),
        sequence_stop_index(stop),
        path_name(path_name),
        biclique_index(pileup_index),
        component_index(component_index)
{}


/// Tell the user which of the sequence indexes is the middlemost index, so that the leftness
/// attribute doesn't need to be considered
size_t SpliceData::get_coordinate(){
    if (is_left){
        return sequence_start_index;
    }
    else{
        // Because divide handle cuts AFTER the index provided, right sites need to be offset by 1
        return sequence_stop_index - 1;
    }
}


/// Tell the user which of the sequence indexes is the middlemost index, so that leftness
/// doesn't need to be considered AND flip the coord if the node was reverse at the time
/// this site was created
size_t SpliceData::get_forward_coordinate(HandleGraph& gfa_graph, size_t node_id){
    auto h = gfa_graph.get_handle(node_id, is_reverse);

    size_t coordinate;

    // If the SpliceData tells us that it is "left", that means the node is on the left of an overlap,
    // which indicates that the overlap is happening on the right end of the node
    if (is_left){
        if (not is_reverse){
            // Because divide handle cuts AFTER the index provided, right sites need to be offset by -1
            coordinate = sequence_start_index - 1;
        }
        else{
            size_t forward_start_index = gfa_graph.get_length(h) - sequence_start_index - 1;
            coordinate = forward_start_index;
        }
    }
    else{
        if (not is_reverse){
            coordinate = sequence_stop_index;
        }
        else{
            size_t forward_start_index = gfa_graph.get_length(h) - sequence_stop_index - 1;

            // Because divide handle cuts AFTER the index provided, right sites need to be offset by -1
            coordinate = forward_start_index - 1;
        }
    }

    return coordinate;
}


bool SpliceData::forward_splice_is_left(){
    bool splice_is_left;

    // If the SpliceData tells us that it is "left", that means the node is on the left of an overlap,
    // which indicates that the overlap is happening on the right end of the node
    if (is_left){
        if (not is_reverse){
            splice_is_left = false;
        }
        else{
            splice_is_left = true;
        }
    }
    else{
        // Because divide handle cuts AFTER the index provided, right sites need to be offset by -1
        if (not is_reverse){
            splice_is_left = true;
        }
        else{
            splice_is_left = false;
        }
    }

    return splice_is_left;
}


ostream& operator<<(ostream& os, SpliceData& alignment_data){
    os << "is_reverse:\t" << alignment_data.is_reverse << '\n';
    os << "is_left:\t" << alignment_data.is_left << '\n';
    os << "start:\t" << alignment_data.sequence_start_index << '\n';
    os << "stop:\t" << alignment_data.sequence_stop_index << '\n';
    os << "path_name:\t" << alignment_data.path_name << '\n';
    os << "biclique_index:\t" << alignment_data.biclique_index << '\n';
    os << "component_index:\t" << alignment_data.component_index << '\n';

    return os;
}


bool SpliceData::operator<(const SpliceData& other) const{
    auto a = sequence_stop_index - sequence_start_index + 1;
    auto b = other.sequence_stop_index - other.sequence_start_index + 1;

    return a < b;
}


void PoaPileup::print_paths(){
    graph.for_each_path_handle([&](path_handle_t path){
        graph.for_each_step_in_path(path, [&](step_handle_t step){
            auto h = graph.get_handle_of_step(step);
            std::cout << graph.get_sequence(h);
        });

        std::cout << '\n';
    });
}


void PoaPileup::update_alignment_data(
        bool is_left,
        handle_t node,
        uint64_t start_index,
        uint64_t stop_index,
        string& path_name,
        size_t component_index){

    int64_t id;

    if (not id_map[!is_left].exists(node)) {
        id = id_map[!is_left].insert(node);
        alignment_data[!is_left].push_back({});
    } else {
        id = id_map[!is_left].get_id(node);
    }

    auto side_string = to_string(is_left);
    auto id_string = to_string(graph.get_id(node));
    auto overlap_index_string = to_string(alignment_data[!is_left][id].size());
    auto biclique_index_string = to_string(biclique_index);
    auto component_index_string = to_string(component_index);

    path_name = component_index_string + '_' +
            biclique_index_string + '_' +
            id_string + '_' +
            to_string(id) + '_' +
            side_string + '_' +
            overlap_index_string;

    auto path_handle = graph.create_path_handle(path_name, false);

    alignment_data[!is_left][id].emplace_back(start_index, stop_index, path_name);
}


void PoaPileup::sort_alignment_data_by_length(){
    for (auto& side: alignment_data){
        for (auto& data: side) {
            sort(data.rbegin(), data.rend());
        }
    }
}


Pileup::Pileup():
        id_map(true)
{}


PoaPileup::PoaPileup():
    id_map({true, true})
{}


}