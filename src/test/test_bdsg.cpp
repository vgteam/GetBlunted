#include "bdsg/packed_graph.hpp"
#include "gfa_to_handle.hpp"
#include "utility.hpp"
#include "topological_sort.hpp"
#include "IncrementalIdMap.hpp"
#include "OverlapMap.hpp"


using std::ifstream;
using std::unordered_map;

using bluntifier::gfa_to_path_handle_graph_in_memory;
using bluntifier::gfa_to_path_handle_graph;
using bluntifier::gfa_to_handle_graph;
using bluntifier::parent_path;
using bluntifier::join_paths;
using bluntifier::orient_nodes_forward;
using bluntifier::IncrementalIdMap;
using bluntifier::OverlapMap;
using bluntifier::Alignment;
using handlegraph::handle_t;
using bdsg::PackedGraph;
using bdsg::MutablePathMutableHandleGraph;


int main(){
    string script_path = __FILE__;
    string project_directory = parent_path(script_path, 3);

    // Get test GFA path
    string relative_gfa_path = "/data/diploid_case_d.gfa";
    const string absolute_gfa_path = join_paths(project_directory, relative_gfa_path);
    ifstream file(absolute_gfa_path);

    PackedGraph g;
    IncrementalIdMap<string> id_map;
    OverlapMap overlaps;

    gfa_to_handle_graph(absolute_gfa_path, g, id_map, overlaps);
    g.clear();

    gfa_to_path_handle_graph_in_memory(file, g, id_map, overlaps);

    orient_nodes_forward(g);



    return 0;
}

