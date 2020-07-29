#ifndef BLUNTIFIER_OVERLAPMAP_HPP
#define BLUNTIFIER_OVERLAPMAP_HPP

#include "handlegraph/handle_graph.hpp"
#include "gfakluge.hpp"
#include "utility.hpp"
#include "Cigar.hpp"
#include <unordered_map>
#include <string>
#include <utility>

using handlegraph::handle_t;
using handlegraph::edge_t;
using bluntifier::Cigar;
using std::unordered_map;
using std::string;
using std::pair;

namespace bluntifier {

class OverlapMap {
public:
    /// Attributes ///
    unordered_map <pair <handle_t, handle_t>, Cigar> overlaps;

    /// Methods ///
    OverlapMap();

    void insert(const gfak::edge_elem& e, handle_t source, handle_t sink);
    void insert(const gfak::edge_elem& e, edge_t& edge_handle);
    Cigar at(handle_t source, handle_t sink);
    Cigar at(edge_t& edge_handle);
    void get_lengths(pair<size_t, size_t>& lengths, edge_t& edge_handle);
};

}

#endif //BLUNTIFIER_OVERLAPMAP_HPP
