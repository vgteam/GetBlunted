#include "Biclique.hpp"

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

}