#ifndef BLUNTIFIER_PILEUP_HPP
#define BLUNTIFIER_PILEUP_HPP

#include "handlegraph/handle_graph.hpp"
#include "bdsg/packed_graph.hpp"
#include "IncrementalIdMap.hpp"
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <deque>
#include <string>
#include <array>

using bluntifier::IncrementalIdMap;
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;
using handlegraph::handle_t;
using handlegraph::edge_t;
using bdsg::PackedGraph;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::string;
using std::deque;
using std::pair;
using std::array;


namespace bluntifier{



class AlignmentData{
public:
    bool is_reverse;
    bool is_left;
    uint64_t sequence_start_index;
    uint64_t sequence_stop_index;
    string path_name;
    size_t pileup_index;
    size_t component_index;
    uint32_t spoa_id;

    AlignmentData()=default;
    AlignmentData(uint64_t start, uint64_t stop, string& path_name);
    AlignmentData(
            bool is_forward,
            bool is_left,
            uint64_t start,
            uint64_t stop,
            string& path_name,
            size_t pileup_index,
            size_t component_index);

    bool operator<(const AlignmentData& other) const;
};


class PoaPileup {
public:
    /// Attributes ///

    // A new graph to represent the multiple sequence alignment
    PackedGraph graph;

    // Convert nodes to integer IDs
    array <IncrementalIdMap<handle_t>, 2> id_map;

    // For each node in the original gfa graph (stored here in order of their appearance) what are the nodes in the
    // pileup graph (from left to right) that will be spliced back into the gfa
    array <vector <vector <AlignmentData> >, 2> alignment_data;

    // Keep track of any edges that are non-overlapping, and cant be represented in the POA graph
    vector <edge_t> blunt_edges;

    // Index to help create unique names for paths in the context of the main GFA graph
    size_t index = 0;


    /// Methods ///
    void print_paths();
    void sort_alignment_data_by_length();
    void update_alignment_data(
            bool is_left,
            handle_t node,
            uint64_t start_index,
            uint64_t stop_index,
            string& path_name,
            size_t component_index);

    PoaPileup();
    PoaPileup(PoaPileup&&) = default;

private:
    const size_t LEFT = 0;
    const size_t RIGHT = 1;
};


class Pileup {
public:
    /// Attributes ///

    // A new graph to represent the multiple sequence alignment
    PackedGraph graph;

    // In terms of nodes in the graph, what path does each alignment take?
    vector <deque <handle_t> > paths;

    // List of edges that have been traversed in order of traversal
    vector <edge_t> edges_traversed;

    // Convert nodes to integer IDs
    IncrementalIdMap<handle_t> id_map;

    /// Methods ///
    Pileup();
};

}

#endif //BLUNTIFIER_PILEUP_HPP
