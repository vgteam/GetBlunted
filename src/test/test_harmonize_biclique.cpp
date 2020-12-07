#include "BluntifierAlign.hpp"
#include "handle_to_gfa.hpp"
#include "utility.hpp"
#include "bdsg/hash_graph.hpp"

using std::cout;
using std::string;
using bdsg::HashGraph;
using handlegraph::handle_t;
using handlegraph::edge_t;
using bluntifier::Bicliques;
using bluntifier::run_command;
using bluntifier::handle_graph_to_gfa;
using bluntifier::harmonize_biclique_orientations;


void print_biclique(HandleGraph& graph, Bicliques& bicliques){
    for (auto& edge: bicliques[0]){
        cout << "(" << graph.get_id(edge.first);
        cout << (graph.get_is_reverse(edge.first) ? "-" : "+");
        cout << ") -> (" << graph.get_id(edge.second);
        cout << (graph.get_is_reverse(edge.second) ? "-" : "+") << ") ";
        cout << graph.get_sequence(edge.first) << " " << graph.get_sequence(edge.second) << '\n';
    }
}


void test_a(){
    HashGraph graph;

    auto a = graph.create_handle("CAA");
    auto f = graph.create_handle("CAA");
    auto b = graph.create_handle("CAA");
    auto c = graph.create_handle("CAA");
    auto d = graph.create_handle("TTG");
    auto e = graph.create_handle("TTG");

    Bicliques bicliques;
    bicliques.bicliques.emplace_back();

    bicliques.bicliques[0] = {
            {a,b},
            {a,c},
            {graph.flip(d),a},
            {graph.flip(a),e},
            {f,b},
            {f,c},
            {graph.flip(d),f},
            {graph.flip(f),e}
    };

    for (auto& edge: bicliques[0]) {
        graph.create_edge(edge.first, edge.second);
    }

    {
        string test_path_prefix = "test_harmonize_biclique_a";
        handle_graph_to_gfa(graph, test_path_prefix + ".gfa");
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    print_biclique(graph, bicliques);
    cout << '\n';

    harmonize_biclique_orientations(graph, bicliques);

    print_biclique(graph, bicliques);
    cout << '\n';


    for (auto& edge: bicliques[0]){
        if (graph.get_sequence(edge.first) != "CAA" or graph.get_sequence(edge.second) != "CAA"){
            throw runtime_error("FAIL");
        }
    }


}

void test_b(){
//    Biclique 5
//    (a+) -> (b-)
//    (a+) -> (c+)
//    (b+) -> (d+)
//    (c-) -> (d+)

//    Biclique 5 (harmonized)
//    (a+) -> (b-)
//    (a+) -> (c+)
//    (b+) -> (d+)
//    (c-) -> (d+)

    HashGraph graph;

    auto a = graph.create_handle("CAA");
    auto d = graph.create_handle("TTG");
    auto b = graph.create_handle("TTG");
    auto c = graph.create_handle("CAA");

    Bicliques bicliques;
    bicliques.bicliques.emplace_back();

    bicliques.bicliques[0] = {
            {a,graph.flip(b)},
            {a,c},
            {b,d},
            {graph.flip(c),d}
    };

    for (auto& edge: bicliques[0]) {
        graph.create_edge(edge.first, edge.second);
    }

    {
        string test_path_prefix = "test_harmonize_biclique_b";
        handle_graph_to_gfa(graph, test_path_prefix + ".gfa");
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    print_biclique(graph, bicliques);
    cout << '\n';

    harmonize_biclique_orientations(graph, bicliques);

    print_biclique(graph, bicliques);
    cout << '\n';

    for (auto& edge: bicliques[0]){
        if (graph.get_sequence(edge.first) != "CAA" or graph.get_sequence(edge.second) != "CAA"){
            throw runtime_error("FAIL");
        }
    }
}


int main(){
    test_a();
    cout << '\n';
    test_b();

    cout << "PASS\n";
    return 0;
}