#include "bdsg/packed_graph.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "utility.hpp"
#include "IncrementalIdMap.hpp"
#include "OverlapMap.hpp"

using handlegraph::handle_t;
using bdsg::PackedGraph;
using std::ifstream;
using std::ofstream;

using bluntifier::gfa_to_path_handle_graph_in_memory;
using bluntifier::gfa_to_path_handle_graph;
using bluntifier::gfa_to_handle_graph;
using bluntifier::handle_graph_to_gfa;
using bluntifier::parent_path;
using bluntifier::join_paths;
using bluntifier::IncrementalIdMap;
using bluntifier::OverlapMap;
using bluntifier::Alignment;


int main(){
    string script_path = __FILE__;
    string project_directory = parent_path(script_path, 3);

    // Get test GFA path
    string relative_gfa_path = "/data/test_gfa1.gfa";
    string absolute_gfa_path = join_paths(project_directory, relative_gfa_path);

    ifstream file(absolute_gfa_path);

    PackedGraph g;
    IncrementalIdMap<string> id_map;
    OverlapMap overlaps;

    gfa_to_handle_graph(absolute_gfa_path, g, id_map, overlaps);

    string output_dir = join_paths(project_directory, "data/");
    string output_path = join_paths(output_dir, "handle_to_gfa_test_output.gfa");
    ofstream out(output_path);
    handle_graph_to_gfa(g, out);

    return 0;
}

