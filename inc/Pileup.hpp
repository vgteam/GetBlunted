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
using handlegraph::HandleGraph;
using handlegraph::handle_t;
using handlegraph::edge_t;
using bdsg::PackedGraph;
using std::unordered_map;
using std::unordered_set;
using std::ostream;
using std::vector;
using std::string;
using std::deque;
using std::pair;
using std::array;


namespace bluntifier{


class SpliceData{
public:
    /// Attributes ///

    // Was the splice site calculated from a node that was flipped in the handle graph
    bool is_reverse;

    // Was the node to the left or right of the biclique (in the GFA-canonical direction)
    bool is_left;

    // At which coordinate is the splice event (orientation sensitive)
    size_t sequence_start_index;
    size_t sequence_stop_index;

    // A globally unique string to point to a path in the overlap POA graph which describes where
    // this node's sequence enters the graph
    string path_name;

    // To which biclique does this site splice into (globally unique index)
    size_t biclique_index;

    // From which adjacency component was this site
    size_t component_index;

    // Only used internally for keeping track of this sequence's place in the SPOA overlap
    uint32_t spoa_id;


    /// Methods ///

    SpliceData()=default;
    SpliceData(uint64_t start, uint64_t stop, string& path_name);
    SpliceData(
            bool is_reverse,
            bool is_left,
            uint64_t start,
            uint64_t stop,
            string& path_name,
            size_t pileup_index,
            size_t component_index);

    // Simplify finding which coord is the relevant one when splicing into the gfa (variable from left/right nodes)
    size_t get_coordinate();

    // Remove left/right AND forward/reverse ambiguity when finding the splice coordinate
    size_t get_forward_coordinate(HandleGraph& gfa_graph, size_t node_id) const;

    // Assume an index is given in forward orientation, then re-zero this splice sites coords based on that index
    void offset_splice_coordinate(HandleGraph& gfa_graph, size_t node_id, size_t coordinate);

    // Save some sanity by telling the user whether the splice site is on the left end of the canonical node
    bool forward_splice_is_left() const;

    bool operator<(const SpliceData& other) const;
};


ostream& operator<<(ostream& os, SpliceData& alignment_data);


class PoaPileup {
public:
    /// Attributes ///

    // A new graph to represent the multiple sequence alignment
    PackedGraph graph;

    // Convert nodes to integer IDs
    array <IncrementalIdMap<handle_t>, 2> id_map;

    // For each node in the original gfa graph (stored here in order of their appearance) what are the nodes in the
    // pileup graph (from left to right) that will be spliced back into the gfa
    array <vector <vector <SpliceData> >, 2> alignment_data;

    // Keep track of any edges that are non-overlapping, and cant be represented in the POA graph
    vector <edge_t> blunt_edges;

    // Index to help create unique names for paths in the context of the main GFA graph
    size_t biclique_index = 0;


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
