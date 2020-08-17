#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>

#include "bdsg/hash_graph.hpp"
#include "BicliqueCover.hpp"

using bluntifier::adjacency_components;
using bluntifier::AdjacencyComponent;
using bluntifier::bipartition;
using bluntifier::ordered_bipartition;
using bluntifier::BipartiteGraph;
using bluntifier::GaloisLattice;
using bluntifier::CenteredGaloisTree;

using bdsg::HashGraph;
using handlegraph::handle_t;
using handlegraph::HandleGraph;
using std::vector;
using std::set;
using std::sort;
using std::cerr;
using std::endl;
using std::function;


int main(){
    
    // dominos can be detected
    {
        HashGraph graph;
        
        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");
        handle_t h4 = graph.create_handle("A");
        handle_t h5 = graph.create_handle("A");
        
        // initialize with an incomplete domino
        graph.create_edge(h0, h3);
        graph.create_edge(h0, h4);
        graph.create_edge(h1, h3);
        graph.create_edge(h1, h4);
        graph.create_edge(h1, h5);
        graph.create_edge(h2, h4);
        
        bipartition partition({h0, h1, h2},
                              {graph.flip(h3), graph.flip(h4), graph.flip(h5)});
        
        BipartiteGraph bigraph(graph, partition);
        
        {
            GaloisLattice lattice(bigraph);
            if (!lattice.is_domino_free()) {
                return 1;
            }
        }
        
        // add edge to make a domino
        graph.create_edge(h2, h5);
        
        {
            GaloisLattice lattice(bigraph);
            if (lattice.is_domino_free()) {
                return 1;
            }
        }
        
        // break the domino
        graph.create_edge(h0, h5);
        
        {
            GaloisLattice lattice(bigraph);
            if (!lattice.is_domino_free()) {
                return 1;
            }
        }
    }
    
    {
         
    }
    
    cerr << "biclique cover tests successful" << endl;
    return 0;
}

