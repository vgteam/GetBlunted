#include "bdsg/packed_graph.hpp"
#include "gfa_to_handle.hpp"
#include "utility.hpp"


using std::ifstream;

using bluntifier::gfa_to_path_handle_graph_in_memory;
using bluntifier::gfa_to_path_handle_graph;
using bluntifier::gfa_to_handle_graph;
using bluntifier::parent_path;
using bluntifier::join_paths;
using handlegraph::handle_t;
using bdsg::PackedGraph;
using bdsg::MutablePathMutableHandleGraph;


int main(){
    string script_path = __FILE__;
    string project_directory = parent_path(script_path, 3);

    // Get test GFA path
    string relative_gfa_path = "/data/test_gfa1.gfa";
    const string absolute_gfa_path = join_paths(project_directory, relative_gfa_path);
    ifstream file(absolute_gfa_path);

    PackedGraph g;

    gfa_to_handle_graph(absolute_gfa_path, g);

    g.clear();

    gfa_to_path_handle_graph_in_memory(file, g);





    return 0;
}

