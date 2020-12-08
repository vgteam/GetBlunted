#include "utility.hpp"
#include "handle_to_gfa.hpp"
#include "bdsg/hash_graph.hpp"

#include <iostream>

using bluntifier::handle_graph_to_gfa;
using bluntifier::run_command;
using bdsg::HashGraph;

using std::string;
using std::vector;
using std::ofstream;

int main(){

    for (size_t o: {0,1,2,3,4}){
        HashGraph graph;
        string original_sequence = "GATT";
        auto h = graph.create_handle(original_sequence);

        vector<size_t> offsets = {o};
        graph.divide_handle(h, offsets);

        {
            string test_path_prefix = "test_divide_handle" + std::to_string(o);
            ofstream out(test_path_prefix + ".gfa");
            handle_graph_to_gfa(graph, out);
            string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                             + test_path_prefix + ".png";
            run_command(command);
        }
    }

    return 0;
}

