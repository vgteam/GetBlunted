#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>

#include "bdsg/hash_graph.hpp"
#include "adjacency_components.hpp"

using bluntifier::adjacency_components;
using bluntifier::AdjacencyComponent;

using bdsg::HashGraph;
using handlegraph::handle_t;
using handlegraph::HandleGraph;
using std::vector;
using std::set;
using std::sort;
using std::cerr;
using std::endl;
using std::function;

int count_edges_across(const AdjacencyComponent& adj_comp,
                       const bipartition& partition) {
    int count = 0;
    for (auto side : adj_comp) {
        adj_comp.for_each_adjacent_side(side, [&](handle_t adj_side) {
            if (partition.first.count(side) == partition.second.count(adj_side)
                && !(side < adj_side)) {
                ++count;
            }
            return true;
        });
    }
    return count;
}

bool check_bipartite(const HandleGraph& graph, const bipartition& partition) {
    
    bool looks_bipartite = true;
    for (auto it = partition.first.begin(); it != partition.first.end() && looks_bipartite; ++it) {
        graph.follow_edges(*it, false, [&](const handle_t& adj_node) {
            handle_t adj_side = graph.flip(adj_node);
            if (partition.first.count(adj_side)) {
                looks_bipartite = false;
                return false;
            }
            if (!partition.second.count(adj_side)) {
                looks_bipartite = false;
                return false;
            }
            return true;
        });
    }
    for (auto it = partition.second.begin(); it != partition.second.end() && looks_bipartite; ++it) {
        graph.follow_edges(*it, false, [&](const handle_t& adj_node) {
            handle_t adj_side = graph.flip(adj_node);
            if (partition.second.count(adj_side)) {
                looks_bipartite = false;
                return false;
            }
            if (!partition.first.count(adj_side)) {
                looks_bipartite = false;
                return false;
            }
            return true;
        });
    }
    return looks_bipartite;
}

int main(){
    
    // tests using a simple SNP graph
    {
        HashGraph graph;
        
        handle_t h1 = graph.create_handle("GAT");
        handle_t h2 = graph.create_handle("T");
        handle_t h3 = graph.create_handle("T");
        handle_t h4 = graph.create_handle("ACA");
        
        graph.create_edge(h1, h2);
        graph.create_edge(h1, graph.flip(h3));
        graph.create_edge(h2, h4);
        graph.create_edge(graph.flip(h3), h4);
        
        // can we correctly find the adjacency components?
        {
            // manually construct the adjacency components
            set<vector<handle_t>> true_adj_comps;
            
            vector<handle_t> comp = {graph.flip(h1)};
            true_adj_comps.insert(comp);
            
            comp = {h1, graph.flip(h2), h3};
            sort(comp.begin(), comp.end());
            true_adj_comps.insert(comp);
            
            comp = {h2, graph.flip(h3), graph.flip(h4)};
            sort(comp.begin(), comp.end());
            true_adj_comps.insert(comp);
            
            comp = {h4};
            true_adj_comps.insert(comp);
            
            // get the adjacency components algorithmically
            auto found_adj_comps = adjacency_components(graph);
            
            // make sure we found all and only the same ones
            // as we constructed manually
            for (auto& found_comp : found_adj_comps) {
                vector<handle_t> comp_sides(found_comp.begin(), found_comp.end());
                sort(comp_sides.begin(), comp_sides.end());
                if (true_adj_comps.count(comp_sides)) {
                    true_adj_comps.erase(comp_sides);
                }
                else {
                    cerr << "found an erroneous adjacency component containing:" << endl;
                    return 1;
                }
            }
            
            if (!true_adj_comps.empty()) {
                cerr << "failed to find an adjacency component" << endl;
                return 1;
            }
        }
        
        // can we correctly identify adjacency components as bipartite?
        {
            for (auto& adj_comp : adjacency_components(graph)) {
                if (!adj_comp.is_bipartite()) {
                    cerr << "erroneously identified a component as non-bipartite" << endl;
                    return 1;
                }
            }
        }
    }
    
    // tests using graph with a non bipartite adjacency component
    {
        
        HashGraph graph;
        
        handle_t h1 = graph.create_handle("GAT");
        handle_t h2 = graph.create_handle("TA");
        handle_t h3 = graph.create_handle("CA");
        
        graph.create_edge(h1, graph.flip(h2));
        graph.create_edge(h1, graph.flip(h3));
        graph.create_edge(h2, graph.flip(h3));
        
        for (auto& adj_comp : adjacency_components(graph)) {
            if (adj_comp.size() > 1) {
                if (adj_comp.is_bipartite()) {
                    cerr << "erroneously identified a component as bipartite" << endl;
                    return 1;
                }
            }
            else {
                if (!adj_comp.is_bipartite()) {
                    cerr << "erroneously identified a component as non-bipartite" << endl;
                    return 1;
                }
            }
        }
    }
    
    // tests using graph with a clear maximum bipartition
    {
        
        HashGraph graph;
        
        handle_t h1 = graph.create_handle("GAT");
        handle_t h2 = graph.create_handle("TA");
        handle_t h3 = graph.create_handle("CA");
        handle_t h4 = graph.create_handle("CAT");
        
        // a biclique (1, 2) -> (3, 4)
        graph.create_edge(h1, h3);
        graph.create_edge(h1, h4);
        graph.create_edge(h2, h3);
        graph.create_edge(h2, h4);
        
        // a single other edge
        graph.create_edge(h1, graph.flip(h2));
        
        for (auto& adj_comp : adjacency_components(graph)) {
            if (adj_comp.size() == 1) {
                // skip over the trivial components
                continue;
            }
            
            {
                auto partition = adj_comp.exhaustive_maximum_bipartite_partition();
                
                if (count_edges_across(adj_comp, partition) != 4) {
                    cerr << "did not identify maximum partition in exhaustive search" << endl;
                    return 1;
                }
            }
            {
                auto partition = adj_comp.maximum_bipartite_partition_apx_1_2();
                adj_comp.refine_apx_partition(partition);
                
                if (count_edges_across(adj_comp, partition) != 4) {
                    cerr << "did not identify maximum partition in local search" << endl;
                    return 1;
                }
            }
        }
    }
    
    // tests using a more complicated graph
    {
        
        HashGraph graph;
        
        handle_t h1 = graph.create_handle("GAT");
        handle_t h2 = graph.create_handle("TA");
        handle_t h3 = graph.create_handle("CA");
        handle_t h4 = graph.create_handle("CAT");
        handle_t h5 = graph.create_handle("GAT");
        handle_t h6 = graph.create_handle("TA");
        handle_t h7 = graph.create_handle("CA");
        handle_t h8 = graph.create_handle("CAT");
        
        // a dense bipartite graph (1, 2, 3, 4) -> (5, 6, 7, 8)
        graph.create_edge(h1, h5);
        graph.create_edge(h1, h6);
        graph.create_edge(h1, h8);
        graph.create_edge(h2, h5);
        graph.create_edge(h2, h6);
        graph.create_edge(h2, h7);
        graph.create_edge(h3, h6);
        graph.create_edge(h3, h7);
        graph.create_edge(h3, h8);
        graph.create_edge(h4, h5);
        graph.create_edge(h4, h7);
        graph.create_edge(h4, h8);
        
        // a few other edges
        graph.create_edge(h1, graph.flip(h4));
        graph.create_edge(h2, graph.flip(h3));
        graph.create_edge(graph.flip(h7), h8);
        
        for (auto& adj_comp : adjacency_components(graph)) {
            if (adj_comp.size() == 1) {
                // skip over the trivial components
                continue;
            }
            
            {
                auto partition = adj_comp.exhaustive_maximum_bipartite_partition();
                
                if (count_edges_across(adj_comp, partition) != 12) {
                    cerr << "did not identify maximum partition in large exhaustive search" << endl;
                    return 1;
                }
            }
            {
                auto partition = adj_comp.maximum_bipartite_partition_apx_1_2(3);
                adj_comp.refine_apx_partition(partition);
                
                if (count_edges_across(adj_comp, partition) != 12) {
                    // TODO: this actually isn't a guarantee, should design a better test
                    cerr << "did not identify maximum partition in large local search" << endl;
                    return 1;
                }
            }
        }
    }
    
    // test using a partition that can only be improved locally using an edge swap
    {
        HashGraph graph;
        
        handle_t h1 = graph.create_handle("GAT");
        handle_t h2 = graph.create_handle("TA");
        handle_t h3 = graph.create_handle("CA");
        handle_t h4 = graph.create_handle("CAT");
        handle_t h5 = graph.create_handle("GAT");
        handle_t h6 = graph.create_handle("TA");
        handle_t h7 = graph.create_handle("CA");
        handle_t h8 = graph.create_handle("CAT");
        
        // the focal edge is 2 -> 7
        graph.create_edge(h1, graph.flip(h2));
        graph.create_edge(h2, h5);
        graph.create_edge(h2, h7);
        graph.create_edge(h4, h7);
        graph.create_edge(graph.flip(h6), h7);
        graph.create_edge(graph.flip(h7), h8);
        
        // add all other edges across to keep the partition rigid
        graph.create_edge(h1, h5);
        graph.create_edge(h1, h6);
        graph.create_edge(h1, h8);
        graph.create_edge(h3, h5);
        graph.create_edge(h3, h6);
        graph.create_edge(h3, h8);
        graph.create_edge(h4, h5);
        graph.create_edge(h4, h6);
        graph.create_edge(h4, h8);
        
        for (auto& adj_comp : adjacency_components(graph)) {
            if (adj_comp.size() == 1) {
                // skip over the trivial components
                continue;
            }
            
            bipartition partition;
            partition.first.insert(h1);
            partition.first.insert(h2);
            partition.first.insert(h3);
            partition.first.insert(h4);
            partition.second.insert(graph.flip(h5));
            partition.second.insert(graph.flip(h6));
            partition.second.insert(graph.flip(h7));
            partition.second.insert(graph.flip(h8));
            
            adj_comp.refine_apx_partition(partition);
            
            if (count_edges_across(adj_comp, partition) != 13) {
                cerr << "did not identify maximum partition in edge swap" << endl;
                return 1;
            }
        }
    }
    
    // check for ability to decompose into bipartite blocks
    {
        HashGraph graph;
        handle_t h1 = graph.create_handle("A");
        handle_t h2 = graph.create_handle("A");
        handle_t h3 = graph.create_handle("A");
        
        // fully connected graph (without reversing self-loops)
        graph.create_edge(h1, h1);
        graph.create_edge(h1, h2);
        graph.create_edge(h1, h3);
        graph.create_edge(h1, graph.flip(h2));
        graph.create_edge(h1, graph.flip(h3));
        graph.create_edge(h2, h2);
        graph.create_edge(h2, h3);
        graph.create_edge(h2, graph.flip(h3));
        graph.create_edge(h3, h3);
        
        bool success = true;
        for (auto& adj_comp : adjacency_components(graph)) {
            adj_comp.decompose_into_bipartite_blocks([&](const HandleGraph& g,
                                                         const bipartition& p) {
                success = success && check_bipartite(g, p);
            });
        }
        if (!success) {
            cerr << "failed to decompose into bipartite blocks" << endl;
            return 1;
        }
    }
    
    cerr << "adjacency components tests successful" << endl;
    return 0;
}

