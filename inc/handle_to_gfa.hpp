#ifndef BLUNTIFIER_HANDLE_TO_GFA_HPP
#define BLUNTIFIER_HANDLE_TO_GFA_HPP

#include "handlegraph/handle_graph.hpp"
#include "gfa_to_handle.hpp"
#include "Filesystem.hpp"

#include <fstream>

using handlegraph::HandleGraph;
using handlegraph::handle_t;
using handlegraph::edge_t;
using ghc::filesystem::path;

using std::runtime_error;
using std::ostream;
using std::string;


namespace bluntifier {


char get_reversal_character(const HandleGraph& graph, const handle_t& node);


void write_node_to_gfa(const HandleGraph& graph, const handle_t& node, ostream& output_file);


void write_edge_to_gfa(const HandleGraph& graph, const edge_t& edge, ostream& output_file);


/// With no consideration for directionality, just dump all the edges/nodes into GFA format
void handle_graph_to_gfa(const HandleGraph& graph, ostream& output_gfa);


// TODO write this method to use the overlaps and id map to write the linkages/sequences in the canonical direction
// using the canonical names as well, wherever possible
void handle_graph_to_canonical_gfa(const HandleGraph& graph, const string& output_path);

void write_bfs_to_gfa(path gfa_path, string start_name, size_t radius);


}

#endif //BLUNTIFIER_HANDLE_TO_GFA_HPP
