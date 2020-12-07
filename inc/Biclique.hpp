#ifndef BLUNTIFIER_BICLIQUE_HPP
#define BLUNTIFIER_BICLIQUE_HPP

#include "handle_graph.hpp"

#include <vector>

using handlegraph::handle_t;
using handlegraph::edge_t;

using std::vector;


namespace bluntifier {


class BicliqueEdgeIndex {
public:
    size_t biclique_index;
    size_t edge_index;

    BicliqueEdgeIndex(size_t biclique, size_t edge);
};


class Bicliques {
public:
    vector <vector <edge_t> > bicliques;

    edge_t& operator[](BicliqueEdgeIndex i);

    const edge_t& operator[](BicliqueEdgeIndex i) const;

    vector <edge_t>& operator[](size_t i);

    const vector <edge_t>& operator[](size_t i) const;

    size_t size() const;
};


}

#endif //BLUNTIFIER_BICLIQUE_HPP


