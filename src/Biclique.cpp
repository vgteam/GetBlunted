#include "Biclique.hpp"

using std::cerr;

namespace bluntifier {


BicliqueEdgeIndex::BicliqueEdgeIndex(size_t biclique, size_t edge) :
        biclique_index(biclique),
        edge_index(edge) {}


size_t Bicliques::size() const {
    return bicliques.size();
}


edge_t& Bicliques::operator[](BicliqueEdgeIndex i) {
    return bicliques[i.biclique_index][i.edge_index];
}


const edge_t& Bicliques::operator[](BicliqueEdgeIndex i) const {
    return bicliques[i.biclique_index][i.edge_index];
}


vector<edge_t>& Bicliques::operator[](size_t i) {
    return bicliques[i];
}


const vector<edge_t>& Bicliques::operator[](size_t i) const {
    return bicliques[i];
}

void Bicliques::print(HandleGraph& graph){
    for (size_t i=0; i<bicliques.size(); i++){
        cerr << "Biclique: " << i << '\n';
        for (auto& e: bicliques[i]){
            cerr << '\t' << graph.get_id(e.first) << "<-->" << graph.get_id(e.second) << '\n';
        }
    }
}

}