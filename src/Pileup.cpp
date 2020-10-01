#include "Pileup.hpp"
#include <algorithm>

using std::sort;


namespace bluntifier{


bool AlignmentData::operator<(const AlignmentData& other) const{
    auto a = sequence_stop_index - sequence_start_index + 1;
    auto b = other.sequence_stop_index - other.sequence_start_index + 1;

    return a < b;
}


void PoaPileup::update_alignment_data(bool is_left, handle_t node, uint64_t start_index, uint64_t stop_index){
    int64_t id;

    if (not id_map[!is_left].exists(node)) {
        id = id_map[!is_left].insert(node);
        alignment_data[!is_left].push_back({});
    } else {
        id = id_map[!is_left].get_id(node);
    }

    alignment_data[!is_left][id].emplace_back(start_index, stop_index);

}


void PoaPileup::sort_alignment_data_by_length(){
    for (auto& side: alignment_data){
        sort(side.rbegin(), side.rend());
    }
}


Pileup::Pileup():
        id_map(true)
{}


PoaPileup::PoaPileup():
    id_map({true, true})
{}


}