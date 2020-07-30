#ifndef BLUNTIFIER_OVERLAPMAP_HPP
#define BLUNTIFIER_OVERLAPMAP_HPP

#include "handlegraph/handle_graph.hpp"
#include "gfakluge.hpp"
#include "utility.hpp"
#include "Cigar.hpp"
#include <unordered_map>
#include <string>
#include <utility>

using handlegraph::HandleGraph;
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
    void insert(const gfak::edge_elem& e, const edge_t& edge_handle);
    unordered_map<edge_t,Cigar>::iterator at(handle_t source, handle_t sink);
    unordered_map<edge_t,Cigar>::iterator at(edge_t& edge_handle);
    unordered_map<edge_t,Cigar>::iterator canonicalize_and_find(edge_t& edge, const HandleGraph& graph);
    void compute_lengths(pair<size_t,size_t>& lengths, unordered_map<edge_t,Cigar>::iterator& iter);
    void canonicalize_and_compute_lengths(pair<size_t,size_t>& lengths, edge_t& edge, const HandleGraph& graph);
};

}

#endif //BLUNTIFIER_OVERLAPMAP_HPP
