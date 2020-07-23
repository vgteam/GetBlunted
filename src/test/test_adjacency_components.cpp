#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

#include "bdsg/hash_graph.hpp"
#include "adjacency_components.hpp"

using bluntifier::adjacency_components;

using bdsg::HashGraph;
using handlegraph::handle_t;
using std::vector;
using std::set;
using std::sort;
using std::cerr;
using std::endl;

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
    
    // tests using graph with a non bipartite adjacency
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
                cerr << "STARTING NON BIPARTITE COMPONENT" << endl;
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
    
    
    
    cerr << "adjacency components tests successful" << endl;
    return 0;
}

