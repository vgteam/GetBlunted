#include "Pileup.hpp"
#include <algorithm>

using std::sort;


namespace bluntifier{


AlignmentData::AlignmentData(uint64_t start, uint64_t stop, string& path_name):
        sequence_start_index(start),
        sequence_stop_index(stop),
        path_name(path_name)
{}


AlignmentData::AlignmentData(
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
    pileup_index(pileup_index),
    component_index(component_index)
{}


ostream& operator<<(ostream& os, AlignmentData& alignment_data){
    os << "is_reverse:\t" << alignment_data.is_reverse << '\n';
    os << "is_left:\t" << alignment_data.is_left << '\n';
    os << "start:\t" << alignment_data.sequence_start_index << '\n';
    os << "stop:\t" << alignment_data.sequence_stop_index << '\n';
    os << "path_name:\t" << alignment_data.path_name << '\n';
    os << "pileup_index:\t" << alignment_data.pileup_index << '\n';
    os << "component_index:\t" << alignment_data.component_index << '\n';

    return os;
}


bool AlignmentData::operator<(const AlignmentData& other) const{
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
    auto pileup_index_string = to_string(index);
    auto component_index_string = to_string(component_index);

    path_name = component_index_string + '_' +
            pileup_index_string + '_' +
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