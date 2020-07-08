#include <experimental/filesystem>
#include <fstream>

#include "bdsg/packed_graph.hpp"
#include "handlegraph/handle_graph.hpp"


using std::experimental::filesystem::path;
using std::ifstream;

using bdsg::PackedGraph;
using handlegraph::handle_t;


int main(){
    PackedGraph g;

    handle_t h = g.create_handle("GATTAGA");

    return 0;
}

