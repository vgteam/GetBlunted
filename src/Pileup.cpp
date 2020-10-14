#include "Pileup.hpp"
#include <algorithm>



using std::sort;


namespace bluntifier{


SpliceData::SpliceData(size_t length, string& path_name):
        length(length),
        path_name(path_name)
{}


SpliceData::SpliceData(
        bool is_reverse,
        bool is_left,
        size_t length,
        string& path_name,
        size_t pileup_index,
        size_t component_index):
        is_reverse(is_reverse),
        is_left(is_left),
        length(length),
        path_name(path_name),
        biclique_index(pileup_index),
        component_index(component_index)
{}


size_t SpliceData::get_start_index(size_t node_length) const{
    size_t index;

    if (is_left){
        if (length > node_length){
            throw runtime_error("ERROR: node length is shorter than overlap length");
        }
        index = node_length - length;
    }
    else{
        index = 0;
    }

    return index;
}


size_t SpliceData::get_stop_index(size_t node_length) const{
    size_t index;

    if (is_left){
        index = node_length - 1;
    }
    else{
        if (length > node_length){
            throw runtime_error("ERROR: node length is shorter than overlap length");
        }
        index = length - 1;
    }

    return index;
}


/// Tell the user which of the sequence indexes is the middlemost index, so that the leftness
/// attribute doesn't need to be considered
size_t SpliceData::get_coordinate(size_t node_length){
    if (is_left){
        return get_start_index(node_length);
    }
    else{
        // Because divide handle cuts AFTER the index provided, right sites need to be offset by 1
        return get_stop_index(node_length) - 1;
    }
}


/// Tell the user which of the sequence indexes is the middlemost index, so that leftness
/// doesn't need to be considered AND flip the coord if the node was reverse at the time
/// this site was created
size_t SpliceData::get_forward_coordinate(size_t node_length) const{
    size_t coordinate;

    // If the SpliceData tells us that it is "left", that means the node is on the left of an overlap,
    // which indicates that the overlap is happening on the right end of the node
    if (is_left){
        if (not is_reverse){
            // Because divide handle cuts AFTER the index provided, right sites need to be offset by -1
            coordinate = node_length - length - 1;
        }
        else{
            coordinate = length - 1;
        }
    }
    else{
        if (not is_reverse){
            coordinate = length - 1;
        }
        else{
            // Because divide handle cuts AFTER the index provided, right sites need to be offset by -1
            coordinate = node_length - length - 1;
        }
    }

    return coordinate;
}


bool SpliceData::forward_splice_is_left() const{
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
    os << "length:\t" << alignment_data.length << '\n';
    os << "path_name:\t" << alignment_data.path_name << '\n';
    os << "biclique_index:\t" << alignment_data.biclique_index << '\n';
    os << "component_index:\t" << alignment_data.component_index << '\n';

    return os;
}


bool SpliceData::operator<(const SpliceData& other) const{
    return length < other.length;
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
        size_t splice_site_index){

    int64_t id;

    if (not id_map[!is_left].exists(node)) {
        id = id_map[!is_left].insert(node);
        splice_site_indexes[!is_left].emplace_back();
    } else {
        id = id_map[!is_left].get_id(node);
    }

    splice_site_indexes[!is_left][id].emplace_back(splice_site_index);
}


Pileup::Pileup():
        id_map(true)
{}


PoaPileup::PoaPileup():
    id_map({true, true})
{}


}