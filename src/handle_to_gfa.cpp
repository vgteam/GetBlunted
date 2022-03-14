#include "IncrementalIdMap.hpp"
#include "handle_to_gfa.hpp"
#include "OverlapMap.hpp"
#include "bdsg/hash_graph.hpp"

#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <queue>


using handlegraph::nid_t;
using bdsg::HashGraph;

using std::runtime_error;
using std::ofstream;
using std::function;
using std::unordered_set;
using std::unordered_map;
using std::queue;

namespace bluntifier {


char get_reversal_character(const HandleGraph& graph, const handle_t& node){
    bool reversed = graph.get_is_reverse(node);

    if (reversed){
        return '-';
    }
    else{
        return '+';
    }
}


void write_node_to_gfa(const HandleGraph& graph, const handle_t& node, ostream& output_file){
    output_file << "S\t" << graph.get_id(node) << '\t' << graph.get_sequence(node) << '\n';
}


void write_node_to_gfa(const HandleGraph& graph, IncrementalIdMap<string>& id_map, const handle_t& node, ostream& output_file){
    output_file << "S\t" << id_map.get_name(graph.get_id(node)) << '\t' << graph.get_sequence(node) << '\n';
}


void write_edge_to_gfa(const HandleGraph& graph, const edge_t& edge, ostream& output_file){
    output_file << "L\t" << graph.get_id(edge.first) << '\t' << get_reversal_character(graph, edge.first) << '\t'
                << graph.get_id(edge.second) << '\t' << get_reversal_character(graph, edge.second) << '\t'
                << "0M" << '\n';
}


void write_edge_to_gfa(const HandleGraph& graph, IncrementalIdMap<string>& id_map, const edge_t& edge, OverlapMap& overlap_map, ostream& output_file){
    auto result = overlap_map.canonicalize_and_find(edge, graph);
    const edge_t& e = result->first;

    output_file << "L\t" << id_map.get_name(graph.get_id(e.first)) << '\t' << get_reversal_character(graph, e.first) << '\t'
                << id_map.get_name(graph.get_id(e.second)) << '\t' << get_reversal_character(graph, e.second) << '\t'
                << result->second.get_cigar_string() << '\n';
}


/// With no consideration for directionality, just dump all the edges/nodes into GFA format
void handle_graph_to_gfa(const HandleGraph& graph, ostream& output_gfa){

    output_gfa << "H\tHVN:Z:1.0\n";

    graph.for_each_handle([&](const handle_t& node){
        write_node_to_gfa(graph, node, output_gfa);
    });

    graph.for_each_edge([&](const edge_t& edge){
        write_edge_to_gfa(graph, edge, output_gfa);
    });
}


// TODO write this method to use the overlaps and id map to write the linkages/sequences in the canonical direction
// using the canonical names as well, wherever possible
void handle_graph_to_canonical_gfa(const HandleGraph& graph, const string& output_path){

    throw runtime_error("ERROR: called unfinished method 'handle_graph_to_canonical_gfa'");
//    ofstream output_gfa(output_path);
//    graph.for_each_handle([&](handle_t& node){
//
//    });
//
//    graph.for_each_edge([&](edge_t& edge){
//
//    });
}


void for_node_in_bfs(
        HandleGraph& graph,
        nid_t start_node,
        size_t radius,
        const function<void(const handle_t& h)>& f){
    unordered_set<nid_t> visited;
    queue<nid_t> q;

    q.emplace(start_node);
    visited.emplace(start_node);

    size_t distance = 0;

    while (not q.empty() and distance < radius) {
        nid_t n = q.front();
        q.pop();

        auto h = graph.get_handle(n);
        f(h);

        graph.follow_edges(h, false, [&](const handle_t& other_handle) {
            auto other_node = graph.get_id(other_handle);

            // Attempt to add this to the set of visited nodes
            auto pass = visited.emplace(other_node).second;

            // Check that this has NOT been visited before queuing it
            if (pass) {
                q.emplace(other_node);
            }
        });

        graph.follow_edges(h, true, [&](const handle_t& other_handle) {
            auto other_node = graph.get_id(other_handle);

            // Attempt to add this to the set of visited nodes
            auto pass = visited.emplace(other_node).second;

            // Check that this has NOT been visited before queuing it
            if (pass) {
                q.emplace(other_node);
            }
        });
    }
}


void for_edge_and_node_in_bfs(
        HandleGraph& graph,
        nid_t start_node,
        size_t radius,
        const function<void(const handle_t& h)>& f_node,
        const function<void(const handle_t& h1, const handle_t& h2)>& f_edge
){
    unordered_set <nid_t> visited_nodes;
    unordered_map <handle_t, unordered_set<handle_t> > visited_edges;

    queue <pair <handle_t, size_t> > q;

    q.emplace(graph.get_handle(start_node), 0);

    bool begin = true;

    while (not q.empty()) {
        handle_t h;
        size_t distance;
        tie(h,distance) = q.front();

        f_node(h);

        q.pop();

        vector<nid_t> other_nodes;

        graph.follow_edges(h, false, [&](const handle_t& other_handle) {
            auto other_node = graph.get_id(other_handle);

            // See if this node has been visited
            bool visited = visited_nodes.count(other_node);

            // Check that this has NOT been visited before queuing it
            if (not visited and distance < radius) {
                q.emplace(other_handle, (distance+1));
                other_nodes.emplace_back(other_node);
            }

            auto edge_result = visited_edges[h].insert(other_handle);

            // If this node is at the periphery, don't report any edges that extend outside the radius.
            // Otherwise, report all edges
            if (edge_result.second) {
                if (distance < radius) {
                    f_edge(h, other_handle);
                }
                else {
                    if (visited){
                        f_edge(h, other_handle);
                    }
                }
            }
        });

        graph.follow_edges(h, true, [&](const handle_t& other_handle) {
            auto other_node = graph.get_id(other_handle);

            // See if this node has been visited
            bool visited = visited_nodes.count(other_node);

            // Check that this has NOT been visited before queuing it
            if (not visited and distance < radius) {
                q.emplace(other_handle, (distance+1));
                other_nodes.emplace_back(other_node);
            }

            auto edge_result = visited_edges[other_handle].emplace(h);

            // If this node is at the periphery, don't report any edges that extend outside the radius.
            // Otherwise, report all edges
            if (edge_result.second) {
                if (distance < radius) {
                    f_edge(other_handle, h);
                }
                else {
                    if (visited){
                        f_edge(other_handle, h);
                    }
                }
            }

        });

        for (auto& n: other_nodes){
            visited_nodes.emplace(n);
        }

        if (begin){
            begin = false;
        }
    }
}


void write_bfs_to_gfa(HandleGraph& gfa_graph, IncrementalIdMap<string>& id_map, OverlapMap& overlap_map, nid_t start_node, size_t radius, ofstream& output_gfa){

    auto f_node = [&](const handle_t& h) {
        write_node_to_gfa(gfa_graph, id_map, h, output_gfa);
    };

    auto f_edge = [&](const handle_t& h1, const handle_t& h2) {
        write_edge_to_gfa(gfa_graph, id_map, std::make_pair(h1, h2), overlap_map, output_gfa);
    };

    for_edge_and_node_in_bfs(gfa_graph, start_node, radius, f_node, f_edge);
}


void write_bfs_to_gfa(path gfa_path, string start_name, size_t radius){
    HashGraph gfa_graph;
    IncrementalIdMap<string> id_map;
    OverlapMap overlaps;

    gfa_to_handle_graph(gfa_path, gfa_graph, id_map, overlaps);

    path output_path = gfa_path;
    output_path.replace_extension("subgraph_" + start_name + "_" + to_string(radius) + ".gfa");

    ofstream output_gfa(output_path);

    nid_t start_node = id_map.get_id(start_name);

    write_bfs_to_gfa(gfa_graph, id_map, overlaps, start_node, radius, output_gfa);
}




}
