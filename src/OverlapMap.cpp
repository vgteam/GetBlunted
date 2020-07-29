#include "OverlapMap.hpp"

using std::make_pair;

namespace bluntifier{

OverlapMap::OverlapMap()=default;


void OverlapMap::insert(const gfak::edge_elem& e, handle_t source, handle_t sink){
    Cigar cigar(e.alignment);

    overlaps.insert({make_pair(source, sink), cigar});
}


void OverlapMap::insert(const gfak::edge_elem& e, edge_t& edge_handle){
    Cigar cigar(e.alignment);

    overlaps.insert({edge_handle, cigar});
}


Cigar OverlapMap::at(handle_t source, handle_t sink){
    return overlaps.at(make_pair(source, sink));
}


Cigar OverlapMap::at(edge_t& edge_handle){
    return overlaps.at(edge_handle);
}


void OverlapMap::get_lengths(pair<size_t, size_t>& lengths, edge_t& edge_handle){
    // Reset lengths
    lengths.first = 0;
    lengths.second = 0;

    for (const auto& c: overlaps.at(edge_handle).operations){
        const uint8_t cigar_code = Cigar::cigar_code[c.second];

        // Assume the right side (sink) node is treated as the "reference" in the cigar
        if (Cigar::is_ref_move[cigar_code]){
            lengths.second++;
        }

        // Assume the left side (source) node is treated as the "query" in the cigar
        if (Cigar::is_query_move[cigar_code]){
            lengths.first++;
        }
    }
}


}