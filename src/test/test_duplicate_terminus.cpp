
#include "duplicate_terminus.hpp"
#include "traverse.hpp"
#include "handle_to_gfa.hpp"
#include "bdsg/hash_graph.hpp"

using bluntifier::duplicate_terminus;
using bluntifier::handle_graph_to_gfa;
using bluntifier::source_sink_paths;
using bdsg::HashGraph;

using std::string;
using std::deque;


int main(){
    HashGraph graph;
    string original_sequence = "GATTACA";
    auto h = graph.create_handle(original_sequence);

    deque<handle_t> children;
    deque<size_t> sizes;

    sizes.emplace_back(3);
    sizes.emplace_back(2);
    sizes.emplace_back(2);
    sizes.emplace_back(1);

    {
        string test_path = "test_duplicated_terminus_" + std::to_string(0) + ".gfa";
        handle_graph_to_gfa(graph, test_path);
    }

    duplicate_terminus(graph, sizes, children, h);

    {
        string test_path = "test_duplicated_terminus_" + std::to_string(1) + ".gfa";
        handle_graph_to_gfa(graph, test_path);
    }

    vector<vector<handle_t>> paths = source_sink_paths(graph);

    for (auto& item: children){
        std::cout << graph.get_id(item) << " " << graph.get_sequence(item) << '\n';
    }

    for (auto& path: paths){
        string sequence;

        for (auto& item: path){
            sequence += graph.get_sequence(item);
        }

        if (sequence != original_sequence){
            throw runtime_error("FAIL: traversal sequence does not equal parent node sequence: "
                                + sequence + " != " + original_sequence);
        }
    }

    return 0;
}



