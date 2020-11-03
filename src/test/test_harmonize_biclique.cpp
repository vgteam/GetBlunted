#include "OverlapAligner.hpp"
#include "bdsg/hash_graph.hpp"

using std::cout;
using bdsg::HashGraph;
using handlegraph::handle_t;
using handlegraph::edge_t;
using bluntifier::Bicliques;
using bluntifier::harmonize_biclique_orientations;


void print_biclique(HandleGraph& graph, Bicliques& bicliques){
    for (auto& edge: bicliques[0]){
        cout << "(" << graph.get_id(edge.first);
        cout << (graph.get_is_reverse(edge.first) ? "-" : "+");
        cout << ") -> (" << graph.get_id(edge.second);
        cout << (graph.get_is_reverse(edge.second) ? "-" : "+") << ")" << '\n';
    }
}


int main(){
    HashGraph graph;

    auto a = graph.create_handle("CAA");
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
            {graph.flip(a),e}
    };

    for (auto& edge: bicliques[0]) {
        graph.create_edge(edge.first, edge.second);
    }

    print_biclique(graph, bicliques);
    cout << '\n';

    harmonize_biclique_orientations(graph, bicliques);

    for (auto& item: bicliques[0]){
        graph.create_edge(item.first, item.second);
    }

    print_biclique(graph, bicliques);
    cout << '\n';


    return 0;
}