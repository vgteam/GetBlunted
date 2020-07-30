#include "bdsg/packed_graph.hpp"
#include "gfa_to_handle.hpp"
#include "utility.hpp"
#include "IncrementalIdMap.hpp"
#include "OverlapMap.hpp"


using std::ifstream;
using std::unordered_map;

using bluntifier::gfa_to_path_handle_graph_in_memory;
using bluntifier::gfa_to_path_handle_graph;
using bluntifier::gfa_to_handle_graph;
using bluntifier::parent_path;
using bluntifier::join_paths;
using bluntifier::IncrementalIdMap;
using bluntifier::OverlapMap;
using bluntifier::Cigar;
using handlegraph::handle_t;
using bdsg::PackedGraph;
using bdsg::MutablePathMutableHandleGraph;


int main(){
    string script_path = __FILE__;
    string project_directory = parent_path(script_path, 3);

    // Get test GFA path
    string relative_gfa_path = "/data/test_cigar.gfa";
    const string absolute_gfa_path = join_paths(project_directory, relative_gfa_path);
    ifstream file(absolute_gfa_path);

    PackedGraph g;
    IncrementalIdMap id_map;
    OverlapMap overlaps;

    gfa_to_handle_graph(absolute_gfa_path, g, id_map, overlaps);

    string left_name;
    string right_name;

    for (auto& item: overlaps.overlaps){
        left_name = id_map.get_name(g.get_id(item.first.first));
        right_name = id_map.get_name(g.get_id(item.first.second));

        cerr << left_name << '-' << g.get_is_reverse(item.first.first) << ','
             << right_name << '-' << g.get_is_reverse(item.first.second) << '\n';
        cerr << item.second << '\n' << '\n';
    }

    pair<size_t,size_t> lengths;
    Cigar cigar("");

    g.for_each_edge([&](edge_t edge) {
        overlaps.canonicalize_and_compute_lengths(lengths, edge, g);

        left_name = id_map.get_name(g.get_id(edge.first));
        right_name = id_map.get_name(g.get_id(edge.second));

        cerr << left_name << '-' << g.get_is_reverse(edge.first) << ','
             << right_name << '-' << g.get_is_reverse(edge.second) << '\n';

        cerr << lengths.first << ',' << lengths.second << '\n';
    });


    return 0;
}
