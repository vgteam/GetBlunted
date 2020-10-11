#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>

#include "bdsg/hash_graph.hpp"
#include "BicliqueCover.hpp"

using bluntifier::bipartition;
using bluntifier::ordered_bipartition;
using bluntifier::BipartiteGraph;
using bluntifier::GaloisLattice;
using bluntifier::CenteredGaloisTree;
using bluntifier::BicliqueCover;
using bluntifier::ReducedDualGraph;

using bdsg::HashGraph;
using handlegraph::handle_t;
using handlegraph::HandleGraph;
using std::vector;
using std::set;
using std::sort;
using std::cerr;
using std::endl;
using std::function;
using std::pair;
using std::max_element;

class TestBicliqueCover : public BicliqueCover{
public:
    TestBicliqueCover(const BipartiteGraph& graph) : BicliqueCover(graph) {}
    using BicliqueCover::biclique_cover_apx;
    using BicliqueCover::lattice_polish;
};

bool verify_biclique_cover(const vector<bipartition>& cover,
                           const bipartition& partition,
                           const BipartiteGraph& bigraph) {
    bool complete = true;
    for (auto side : partition.first) {
        bigraph.for_each_adjacent_side(side, [&](handle_t adj) {
            bool found = false;
            for (auto& biclique : cover) {
                if (biclique.first.count(side) && biclique.second.count(adj)) {
                    found = true;
                }
            }
            if (!found) {
                complete = false;
            }
        });
    }
    bool consistent = true;
    for (auto& biclique : cover) {
        for (auto side : biclique.first) {
            bool found = false;
            for (auto other_side : biclique.second) {
                bigraph.for_each_adjacent_side(side, [&](handle_t adj) {
                    if (adj == other_side) {
                        found = true;
                    }
                });
            }
            if (!found) {
                consistent = false;
            }
        }
    }
    return consistent && complete;
}

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


        {
            BipartiteGraph bigraph(graph, partition);
            GaloisLattice lattice(bigraph);
            if (!lattice.is_domino_free()) {
                return 1;
            }
        }

        // add edge to make a domino
        graph.create_edge(h2, h5);

        {
            BipartiteGraph bigraph(graph, partition);
            GaloisLattice lattice(bigraph);
            if (lattice.is_domino_free()) {
                return 1;
            }
        }

        // break the domino
        graph.create_edge(h0, h5);

        {
            BipartiteGraph bigraph(graph, partition);
            GaloisLattice lattice(bigraph);
            if (!lattice.is_domino_free()) {
                return 1;
            }
        }
    }

    // biclique covers can be detected on simple domino free graphs
    {
        HashGraph graph;

        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");

        graph.create_edge(h0, h1);
        graph.create_edge(h0, h2);

        // TODO: these actually aren't simple graphs, so correctness isn't
        // guaranteed

        bipartition partition({h0}, {graph.flip(h1), graph.flip(h2)});


        // a single biclique
        {
            BipartiteGraph bigraph(graph, partition);
            GaloisLattice lattice(bigraph);
            vector<bipartition> separator = lattice.biclique_separator();
            if (separator.size() != 1) {
                return 1;
            }
            if (separator[0].first.size() != 1) {
                return 1;
            }
            if (separator[0].second.size() != 2) {
                return 1;
            }
            if (!separator[0].first.count(h0)) {
                return 1;
            }
            if (!separator[0].second.count(graph.flip(h1))) {
                return 1;
            }
            if (!separator[0].second.count(graph.flip(h2))) {
                return 1;
            }
        }

        // two independent bicliques
        handle_t h3 = graph.create_handle("A");
        handle_t h4 = graph.create_handle("A");
        handle_t h5 = graph.create_handle("A");

        graph.create_edge(h3, h5);
        graph.create_edge(h4, h5);

        partition.first.insert(h3);
        partition.first.insert(h4);
        partition.second.insert(graph.flip(h5));
        {
            BipartiteGraph bigraph(graph, partition);
            GaloisLattice lattice(bigraph);
            vector<bipartition> separator = lattice.biclique_separator();


            if (separator.size() != 2) {
                return 1;
            }
            bool found1 = false, found2 = false;
            for (auto& biclique : separator) {
                bool is1 = true, is2 = true;
                if (biclique.first.size() != 1 ||
                    biclique.second.size() != 2 ||
                    !biclique.first.count(h0) ||
                    !biclique.second.count(graph.flip(h1)) ||
                    !biclique.second.count(graph.flip(h2))) {
                    is1 = false;
                }
                if (biclique.first.size() != 2 ||
                    biclique.second.size() != 1 ||
                    !biclique.first.count(h3) ||
                    !biclique.first.count(h4) ||
                    !biclique.second.count(graph.flip(h5))) {
                    is2 = false;
                }
                found1 = found1 || is1;
                found2 = found2 || is2;
            }
            if (!found1) {
                return 1;
            }
            if (!found2) {
                return 1;
            }

//            BipartiteGraph bg(graph, partition);
//            BicliqueCover c(bg);
//            auto s = c.simplify();
//            s.for_each_handle([&](const handle_t& h) {
//                cerr << s.get_id(h) << endl;
//                s.follow_edges(h, false, [&](const handle_t& n) {
//                    cerr << "\t->" << s.get_id(n) << endl;
//                });
//            });
        }
    }

    // apply the algorithm to graphs that have non-trivial biclique containment
    // relationships
    {
        // a graph with a chain of bicliques
        HashGraph graph;

        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");
        handle_t h4 = graph.create_handle("A");
        handle_t h5 = graph.create_handle("A");

        graph.create_edge(h0, h3);
        graph.create_edge(h0, h4);
        graph.create_edge(h1, h3);
        graph.create_edge(h1, h4);
        graph.create_edge(h1, h5);
        graph.create_edge(h2, h4);

        bipartition partition({h0, h1, h2},
                              {graph.flip(h3), graph.flip(h4), graph.flip(h5)});

        // this graph isn't simple, and it won't won't produce a correct
        // cover, so this mostly just so i can see the debug output
        BipartiteGraph bigraph(graph, partition);
        GaloisLattice lattice(bigraph);
        lattice.biclique_separator();
    }
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
        lattice.biclique_separator();
    }
    
    // simplifying domino free graphs produces the expected results
    {

        HashGraph graph;

        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");

        graph.create_edge(h0, h2);
        graph.create_edge(h1, h3);

        bipartition partition({h0, h1}, {graph.flip(h2), graph.flip(h3)});

        // simplification doesn't do anything to a simple graph

        {
            BipartiteGraph bigraph(graph, partition);
            vector<pair<handle_t, vector<handle_t>>> simplifications;
            BipartiteGraph simple = bigraph.simplify(simplifications);

            if (!simplifications.empty()) {
                return 1;
            }
            int edge_count = 0;
            bool found1 = false, found2 = false;
            for (auto it = simple.left_begin(); it != simple.left_end(); ++it) {
                if (*it == h0) {
                    simple.for_each_adjacent_side(*it, [&](handle_t adj) {
                        edge_count++;
                        if (adj == graph.flip(h2)) {
                            found1 = true;
                        }
                    });
                }
                else if (*it == h1) {
                    simple.for_each_adjacent_side(*it, [&](handle_t adj) {
                        edge_count++;
                        if (adj == graph.flip(h3)) {
                            found2 = true;
                        }
                    });
                }
                else {
                    return 1;
                }
            }
            if (!found1 || !found2 || edge_count != 2) {
                return 1;
            }
        }

        // if we make a successor, the same graph is returned
        graph.create_edge(h0, h3);

        {
            BipartiteGraph bigraph(graph, partition);
            vector<pair<handle_t, vector<handle_t>>> simplifications;
            BipartiteGraph simple = bigraph.simplify(simplifications);

            if (simplifications.size() != 1) {
                return 1;
            }
            if (simplifications[0].first != h1) {
                return 1;
            }
            if (simplifications[0].second[0] != h0) {
                return 1;
            }
            int edge_count = 0;
            bool found1 = false, found2 = false;
            for (auto it = simple.left_begin(); it != simple.left_end(); ++it) {
                if (*it == h0) {
                    simple.for_each_adjacent_side(*it, [&](handle_t adj) {
                        edge_count++;
                        if (adj == graph.flip(h2)) {
                            found1 = true;
                        }
                    });
                }
                else if (*it == h1) {
                    simple.for_each_adjacent_side(*it, [&](handle_t adj) {
                        edge_count++;
                        if (adj == graph.flip(h3)) {
                            found2 = true;
                        }
                    });
                }
                else {
                    return 1;
                }
            }
            if (!found1 || !found2 || edge_count != 2) {
                return 1;
            }
        }
    }
    
    // a more complicated example
    {
        HashGraph graph;
        
        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");
        handle_t h4 = graph.create_handle("A");
        handle_t h5 = graph.create_handle("A");
        handle_t h6 = graph.create_handle("A");
        
        graph.create_edge(h0, h3);
        graph.create_edge(h0, h4);
        graph.create_edge(h0, h5);
        graph.create_edge(h0, h6);
        graph.create_edge(h1, h3);
        graph.create_edge(h1, h4);
        graph.create_edge(h2, h3);
        graph.create_edge(h2, h6);
        
        bipartition partition({h0, h1, h2},
                              {graph.flip(h3), graph.flip(h4), graph.flip(h5), graph.flip(h6)});
        
        BipartiteGraph bigraph(graph, partition);
        vector<pair<handle_t, vector<handle_t>>> simplifications;
        BipartiteGraph simple = bigraph.simplify(simplifications);

        int edge_count = 0;
        bool found1 = false, found2 = false, found3 = false;
        for (auto it = simple.left_begin(); it != simple.left_end(); ++it) {
            if (*it == h0) {
                simple.for_each_adjacent_side(*it, [&](handle_t adj) {
                    edge_count++;
                    if (adj == graph.flip(h5)) {
                        found1 = true;
                    }
                });
            }
            else if (*it == h1) {
                simple.for_each_adjacent_side(*it, [&](handle_t adj) {
                    edge_count++;
                    if (adj == graph.flip(h4)) {
                        found2 = true;
                    }
                });
            }
            else if (*it == h2) {
                simple.for_each_adjacent_side(*it, [&](handle_t adj) {
                    edge_count++;
                    if (adj == graph.flip(h3) || adj == graph.flip(h6)) {
                        found3 = true;
                    }
                });
            }
            else {
                return 1;
            }
        }
        if (!found1 || !found2 || !found3 || edge_count != 3) {
            return 1;
        }
    }
    
    // full end to end run for domino free cover
    {
        HashGraph graph;
        
        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");
        handle_t h4 = graph.create_handle("A");
        handle_t h5 = graph.create_handle("A");
        
        graph.create_edge(h0, h3);
        graph.create_edge(h0, h4);
        graph.create_edge(h0, h5);
        graph.create_edge(h1, h3);
        graph.create_edge(h1, h4);
        graph.create_edge(h2, h3);
        
        bipartition partition({h0, h1, h2},
                              {graph.flip(h3), graph.flip(h4), graph.flip(h5)});
        
        BipartiteGraph bigraph(graph, partition);
        
        auto cover = BicliqueCover(bigraph).get();
        
        if (cover.size() != 3) {
            return 1;
        }
        
        if (!verify_biclique_cover(cover, partition, bigraph)) {
            return 1;
        }
    }
    
    // the dual graph produces the correct answer on a graph it can
    // completely decompose
    {
        HashGraph graph;
        
        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");
        
        graph.create_edge(h0, h2);
        graph.create_edge(h1, h2);
        graph.create_edge(h1, h3);
        
        bipartition partition({h0, h1},
                              {graph.flip(h2), graph.flip(h3)});
        
        BipartiteGraph bigraph(graph, partition);
        
        ReducedDualGraph dual(bigraph);
        
        bool is_exact;
        auto cover = dual.biclique_cover(is_exact);
        
        if (!is_exact) {
            return 1;
        }
        
        if (cover.size() != 2) {
            return 1;
        }
        
        if (!verify_biclique_cover(cover, partition, bigraph)) {
            return 1;
        }
    }
    
    // the dual graph produces a correct biclique cover on an
    // irreducible graph
    {
        HashGraph graph;
        
        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");
        handle_t h4 = graph.create_handle("A");
        handle_t h5 = graph.create_handle("A");
        
        graph.create_edge(h0, h3);
        graph.create_edge(h0, h5);
        graph.create_edge(h1, h3);
        graph.create_edge(h1, h4);
        graph.create_edge(h2, h4);
        graph.create_edge(h2, h5);
        
        bipartition partition({h0, h1, h2},
                              {graph.flip(h3), graph.flip(h4), graph.flip(h5)});
        
        BipartiteGraph bigraph(graph, partition);
        
        ReducedDualGraph dual(bigraph);
        
        bool is_exact;
        auto cover = dual.biclique_cover(is_exact);
        
        if (cover.size() != 3) {
            return 1;
        }
        
        if (!verify_biclique_cover(cover, partition, bigraph)) {
            return 1;
        }
        
    }
    
    // test out the fast heuristic for biclique cover
    {
        HashGraph graph;
        
        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");
        handle_t h4 = graph.create_handle("A");
        handle_t h5 = graph.create_handle("A");
        
        graph.create_edge(h0, h3);
        graph.create_edge(h0, h4);
        graph.create_edge(h1, h3);
        graph.create_edge(h1, h4);
        graph.create_edge(h1, h5);
        graph.create_edge(h2, h4);
        graph.create_edge(h2, h5);
        
        bipartition partition({h0, h1, h2},
                              {graph.flip(h3), graph.flip(h4), graph.flip(h5)});
        
        BipartiteGraph bigraph(graph, partition);
                
        auto cover = TestBicliqueCover(bigraph).biclique_cover_apx();
                
        if (!verify_biclique_cover(cover, partition, bigraph)) {
            return 1;
        }

    }
    // another test that will make it start from the right
    {
        HashGraph graph;
        
        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        
        graph.create_edge(h0, h2);
        graph.create_edge(h1, h2);
        
        bipartition partition({h0, h1},
                              {graph.flip(h2)});
        
        BipartiteGraph bigraph(graph, partition);
        
        auto cover = TestBicliqueCover(bigraph).biclique_cover_apx();
        
        if (!verify_biclique_cover(cover, partition, bigraph)) {
            return 1;
        }

    }
    
    // test the lattice polishing
    {
        
        HashGraph graph;
        
        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");
        handle_t h4 = graph.create_handle("A");
        handle_t h5 = graph.create_handle("A");
        handle_t h6 = graph.create_handle("A");
        
        graph.create_edge(h0, h3);
        graph.create_edge(h0, h4);
        graph.create_edge(h1, h3);
        graph.create_edge(h1, h4);
        graph.create_edge(h1, h5);
        graph.create_edge(h1, h6);
        graph.create_edge(h2, h5);
        graph.create_edge(h2, h6);
        
        bipartition partition({h0, h1, h2},
                              {graph.flip(h3), graph.flip(h4), graph.flip(h5), graph.flip(h6)});
        
        BipartiteGraph bigraph(graph, partition);
        
        vector<bipartition> cover(3);
        cover[0].first.insert(h0);
        cover[0].second.insert(graph.flip(h3));
        cover[0].second.insert(graph.flip(h4));
        cover[1].first.insert(h1);
        cover[1].second.insert(graph.flip(h3));
        cover[1].second.insert(graph.flip(h4));
        cover[1].second.insert(graph.flip(h5));
        cover[1].second.insert(graph.flip(h6));
        cover[2].first.insert(h2);
        cover[2].second.insert(graph.flip(h5));
        cover[2].second.insert(graph.flip(h6));
        
        if (!verify_biclique_cover(cover, partition, bigraph)) {
            return 1;
        }
        
        TestBicliqueCover biclique_cover(bigraph);
        biclique_cover.lattice_polish(cover);
        
        if (!verify_biclique_cover(cover, partition, bigraph)) {
            return 1;
        }
        
        if (cover.size() != 2) {
            return 1;
        }
    }
    
    // lattice polishing works leftward as well
    {
        
        HashGraph graph;
        
        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        
        graph.create_edge(h0, h1);
        graph.create_edge(h0, h2);
        
        bipartition partition({h0},
                              {graph.flip(h1), graph.flip(h2)});
        
        BipartiteGraph bigraph(graph, partition);
        
        vector<bipartition> cover(2);
        cover[0].first.insert(h0);
        cover[0].second.insert(graph.flip(h1));
        cover[1].first.insert(h0);
        cover[1].second.insert(graph.flip(h2));
        
        if (!verify_biclique_cover(cover, partition, bigraph)) {
            return 1;
        }
        
        TestBicliqueCover biclique_cover(bigraph);
        biclique_cover.lattice_polish(cover);
        
        if (!verify_biclique_cover(cover, partition, bigraph)) {
            return 1;
        }
        
        if (cover.size() != 1) {
            return 1;
        }
        
//        for (auto bc : cover) {
//            cerr << "covering biclique:" << endl;
//            cerr << "\tleft:" << endl;
//            for (auto h : bc.first) {
//                cerr << "\t\t" << graph.get_id(h) << " " << graph.get_is_reverse(h) << endl;
//            }
//            cerr << "\tright:" << endl;
//            for (auto h : bc.second) {
//                cerr << "\t\t" << graph.get_id(h) << " " << graph.get_is_reverse(h) << endl;
//            }
//        }
    }
    
    // a simplification case that revealed a bug at one point
    {
        HashGraph graph;
        
        handle_t h0 = graph.create_handle("A");
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");
        handle_t h4 = graph.create_handle("A");
        
        graph.create_edge(h0, h3);
        graph.create_edge(h0, h4);
        graph.create_edge(h1, h3);
        graph.create_edge(h1, h4);
        graph.create_edge(h2, h3);
        graph.create_edge(h2, h4);
        
        bipartition partition({h0, h1, h2},
                              {graph.flip(h3), graph.flip(h4)});
        
        BipartiteGraph bigraph(graph, partition);
        vector<pair<handle_t, vector<handle_t>>> simplifications;
        BipartiteGraph simple = bigraph.simplify(simplifications);
        
        bool found_unsimplified_node = false;
        for (auto h : partition.first) {
            if (!found_unsimplified_node) {
                set<handle_t> sides;
                simple.for_each_adjacent_side(h, [&](const handle_t& side){
                    sides.insert(side);
                });
                if (!sides.empty()) {
                    assert(sides.size() == 1);
                    found_unsimplified_node = true;
                }
            }
            else {
                simple.for_each_adjacent_side(h, [&](const handle_t& side){
                    assert(false);
                });
            }
        }
    }
    
    cerr << "biclique cover tests successful" << endl;
    return 0;
}

