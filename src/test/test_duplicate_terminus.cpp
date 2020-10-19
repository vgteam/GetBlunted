
#include "duplicate_terminus.hpp"
#include "handle_to_gfa.hpp"
#include "bdsg/hash_graph.hpp"

using bluntifier::duplicate_terminus;
using bluntifier::handle_graph_to_gfa;
using bdsg::HashGraph;

using std::string;
using std::deque;


int main(){
    HashGraph graph;

    auto h = graph.create_handle("GATTACA");

    deque<handle_t> children;
    queue<size_t> sizes;

    sizes.emplace(3);
    sizes.emplace(2);
    sizes.emplace(1);

    // TODO: disallow unsorted values?

    {
        string test_path = "test_duplicated_terminus_" + std::to_string(0) + ".gfa";
        handle_graph_to_gfa(graph, test_path);
    }

    duplicate_terminus(graph, sizes, children, h);

    {
        string test_path = "test_duplicated_terminus_" + std::to_string(1) + ".gfa";
        handle_graph_to_gfa(graph, test_path);
    }

    return 0;
}



