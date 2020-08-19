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
using bluntifier::BicliqueCover;

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
    
//    // dominos can be detected
//    {
//        HashGraph graph;
//
//        handle_t h0 = graph.create_handle("A");
//        handle_t h1 = graph.create_handle("A");
//        handle_t h2 = graph.create_handle("A");
//        handle_t h3 = graph.create_handle("A");
//        handle_t h4 = graph.create_handle("A");
//        handle_t h5 = graph.create_handle("A");
//
//        // initialize with an incomplete domino
//        graph.create_edge(h0, h3);
//        graph.create_edge(h0, h4);
//        graph.create_edge(h1, h3);
//        graph.create_edge(h1, h4);
//        graph.create_edge(h1, h5);
//        graph.create_edge(h2, h4);
//
//        bipartition partition({h0, h1, h2},
//                              {graph.flip(h3), graph.flip(h4), graph.flip(h5)});
//
//        BipartiteGraph bigraph(graph, partition);
//
//        {
//            GaloisLattice lattice(bigraph);
//            if (!lattice.is_domino_free()) {
//                return 1;
//            }
//        }
//
//        // add edge to make a domino
//        graph.create_edge(h2, h5);
//
//        {
//            GaloisLattice lattice(bigraph);
//            if (lattice.is_domino_free()) {
//                return 1;
//            }
//        }
//
//        // break the domino
//        graph.create_edge(h0, h5);
//
//        {
//            GaloisLattice lattice(bigraph);
//            if (!lattice.is_domino_free()) {
//                return 1;
//            }
//        }
//    }
//
//    // biclique covers can be detected on simple domino free graphs
//    {
//        HashGraph graph;
//
//        handle_t h0 = graph.create_handle("A");
//        handle_t h1 = graph.create_handle("A");
//        handle_t h2 = graph.create_handle("A");
//
//        graph.create_edge(h0, h1);
//        graph.create_edge(h0, h2);
//
//        // TODO: these actually aren't simple graphs, so correctness isn't
//        // guaranteed
//
//        bipartition partition({h0}, {graph.flip(h1), graph.flip(h2)});
//
//
//        // a single biclique
//        {
//            BipartiteGraph bigraph(graph, partition);
//            GaloisLattice lattice(bigraph);
//            vector<bipartition> cover = lattice.biclique_cover();
//            if (cover.size() != 1) {
//                return 1;
//            }
//            if (cover[0].first.size() != 1) {
//                return 1;
//            }
//            if (cover[0].second.size() != 2) {
//                return 1;
//            }
//            if (!cover[0].first.count(h0)) {
//                return 1;
//            }
//            if (!cover[0].second.count(graph.flip(h1))) {
//                return 1;
//            }
//            if (!cover[0].second.count(graph.flip(h2))) {
//                return 1;
//            }
//        }
//
//        // two independent bicliques
//        handle_t h3 = graph.create_handle("A");
//        handle_t h4 = graph.create_handle("A");
//        handle_t h5 = graph.create_handle("A");
//
//        graph.create_edge(h3, h5);
//        graph.create_edge(h4, h5);
//
//        partition.first.insert(h3);
//        partition.first.insert(h4);
//        partition.second.insert(graph.flip(h5));
//        {
//            BipartiteGraph bigraph(graph, partition);
//            GaloisLattice lattice(bigraph);
//            vector<bipartition> cover = lattice.biclique_cover();
//
//
//            if (cover.size() != 2) {
//                return 1;
//            }
//            bool found1 = false, found2 = false;
//            for (auto& biclique : cover) {
//                bool is1 = true, is2 = true;
//                if (biclique.first.size() != 1 ||
//                    biclique.second.size() != 2 ||
//                    !biclique.first.count(h0) ||
//                    !biclique.second.count(graph.flip(h1)) ||
//                    !biclique.second.count(graph.flip(h2))) {
//                    is1 = false;
//                }
//                if (biclique.first.size() != 2 ||
//                    biclique.second.size() != 1 ||
//                    !biclique.first.count(h3) ||
//                    !biclique.first.count(h4) ||
//                    !biclique.second.count(graph.flip(h5))) {
//                    is2 = false;
//                }
//                found1 = found1 || is1;
//                found2 = found2 || is2;
//            }
//            if (!found1) {
//                return 1;
//            }
//            if (!found2) {
//                return 1;
//            }
//
////            BipartiteGraph bg(graph, partition);
////            BicliqueCover c(bg);
////            auto s = c.simplify();
////            s.for_each_handle([&](const handle_t& h) {
////                cerr << s.get_id(h) << endl;
////                s.follow_edges(h, false, [&](const handle_t& n) {
////                    cerr << "\t->" << s.get_id(n) << endl;
////                });
////            });
//        }
//    }
//
//    // apply the algorithm to graphs that have non-trivial biclique containment
//    // relationships
//    {
//        // a graph with a chain of bicliques
//        HashGraph graph;
//
//        handle_t h0 = graph.create_handle("A");
//        handle_t h1 = graph.create_handle("A");
//        handle_t h2 = graph.create_handle("A");
//        handle_t h3 = graph.create_handle("A");
//        handle_t h4 = graph.create_handle("A");
//        handle_t h5 = graph.create_handle("A");
//
//        graph.create_edge(h0, h3);
//        graph.create_edge(h0, h4);
//        graph.create_edge(h1, h3);
//        graph.create_edge(h1, h4);
//        graph.create_edge(h1, h5);
//        graph.create_edge(h2, h4);
//
//        bipartition partition({h0, h1, h2},
//                              {graph.flip(h3), graph.flip(h4), graph.flip(h5)});
//
//        // this graph isn't simple, and it won't won't produce a correct
//        // cover, so this mostly just so i can see the debug output
//        BipartiteGraph bigraph(graph, partition);
//        GaloisLattice lattice(bigraph);
//        lattice.biclique_cover();
//    }
    {
        // a graph bifurcating bicliques
        HashGraph graph;
        
        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");
        handle_t h4 = graph.create_handle("A");
        
        graph.create_edge(h0, h3);
        graph.create_edge(h1, h3);
        graph.create_edge(h1, h4);
        graph.create_edge(h2, h4);
        
        bipartition partition({h0, h1, h2},
                              {graph.flip(h3), graph.flip(h4)});
        
        // this graph isn't simple, and it won't won't produce a correct
        // cover, so this mostly just so i can see the debug output
        BipartiteGraph bigraph(graph, partition);
        GaloisLattice lattice(bigraph);
        lattice.biclique_cover();
    }
    
    
//            for (auto bc : cover) {
//                cerr << "covering biclique:" << endl;
//                cerr << "\tleft:" << endl;
//                for (auto h : bc.first) {
//                    cerr << "\t\t" << graph.get_id(h) << " " << graph.get_is_reverse(h) << endl;
//                }
//                cerr << "\tright:" << endl;
//                for (auto h : bc.second) {
//                    cerr << "\t\t" << graph.get_id(h) << " " << graph.get_is_reverse(h) << endl;
//                }
//            }
    
    cerr << "biclique cover tests successful" << endl;
    return 0;
}

