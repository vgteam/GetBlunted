#include "AdjacencyComponent.hpp"
#include "BicliqueCover.hpp"
#include "OverlapMap.hpp"
#include "gfa_to_handle.hpp"
#include "copy_graph.hpp"
#include "utility.hpp"

#include "bdsg/hash_graph.hpp"

using bluntifier::gfa_to_path_handle_graph_in_memory;
using bluntifier::gfa_to_path_handle_graph;
using bluntifier::gfa_to_handle_graph;
using bluntifier::parent_path;
using bluntifier::join_paths;
using bluntifier::IncrementalIdMap;
using bluntifier::OverlapMap;
using bluntifier::Alignment;
using bluntifier::for_each_adjacency_component;
using bluntifier::AdjacencyComponent;
using bluntifier::BipartiteGraph;
using bluntifier::BicliqueCover;
using bluntifier::bipartition;
using bluntifier::copy_path_handle_graph;

using handlegraph::handle_t;
using bdsg::HashGraph;


void deduplicate_and_canonicalize_edges_in_bicliques(
        vector<bipartition>& biclique_cover,
        vector<edge_t>& deduplicated_biclique,
        const HandleGraph& gfa_graph,
        const OverlapMap& overlaps){

    // sort the bicliques in descending order by size (to get any repeated edges
    // into larger POAs -- likely to be more compact this way)
    sort(biclique_cover.begin(), biclique_cover.end(),
         [&](const bipartition& a, const bipartition& b) {
             return a.first.size() * a.second.size() > b.first.size() * b.second.size();
         });

    unordered_set<edge_t> edges_processed;
    for (const bipartition& biclique : biclique_cover) {
        // get the edges that haven't been handled in a previous biclique
        for (handle_t left : biclique.first) {
            for (handle_t right : biclique.second) {
                edge_t edge(left, gfa_graph.flip(right));
                auto iter = overlaps.canonicalize_and_find(edge, gfa_graph);

                if (!edges_processed.count(edge)) {
                    edges_processed.emplace(edge);
                    deduplicated_biclique.emplace_back(iter->first);
                }
            }
        }
    }
}


void compute_all_bicliques(
        size_t i,
        const HashGraph& gfa_graph,
        const OverlapMap& overlaps,
        vector<AdjacencyComponent>& adjacency_components,
        vector <vector <edge_t> >& bicliques,
        mutex& biclique_mutex){

    auto& adjacency_component = adjacency_components[i];

    // Skip trivial adjacency components (dead ends)
    if (adjacency_component.size() == 1) {
        return;
    }

    adjacency_component.decompose_into_bipartite_blocks([&](const BipartiteGraph& bipartite_graph){
        biclique_mutex.lock();
        bicliques.emplace_back();
        biclique_mutex.unlock();

        vector<bipartition> biclique_cover = BicliqueCover(bipartite_graph).get();

        deduplicate_and_canonicalize_edges_in_bicliques(biclique_cover, bicliques.back(), gfa_graph, overlaps);
    });
}


void find_middlemost_overlaps(
        const vector <vector <pair <size_t, size_t> > >& node_to_biclique_index,
        const vector <vector <edge_t> >& bicliques,
        HandleGraph& handle_graph,
        size_t node_id){

    // Map each biclique to its start index w.r.t. to the forward parent node
    map <size_t, size_t> left_overlaps;
    map <size_t, size_t> right_overlaps;


    for (auto& item: node_to_biclique_index[node_id]){
        // Do stuff
    }
}



void duplicate_termini(
        const vector <vector <pair <size_t, size_t> > >& node_to_biclique_index,
        const vector <vector <edge_t> >& bicliques,
        HandleGraph& handle_graph){

    for (size_t node_id=1; node_id<node_to_biclique_index.size(); node_id++){

    }
}


void map_splice_sites_by_node(
        const HandleGraph& gfa_graph,
        const vector <vector <edge_t> >& bicliques,
        vector <vector <pair <size_t, size_t> > >& node_to_biclique_index){

    // Create a mapping from all the nodes to their participating edges in each biclique, where the mapping
    // just keeps track of the biclique index and the intra-biclique index for each edge in the
    // "bicliques" vector of vectors, using a pair of indexes {bc_index, ibc_index}
    for (size_t i=0; i<bicliques.size(); i++){
        for (size_t j=0; j<bicliques[i].size(); j++){
            const auto& edge = bicliques[i][j];

            nid_t left_node_id;
            nid_t right_node_id;

            left_node_id = gfa_graph.get_id(edge.first);
            right_node_id = gfa_graph.get_id(edge.second);

            node_to_biclique_index[left_node_id].emplace_back(i,j);
            node_to_biclique_index[right_node_id].emplace_back(i,j);
        }
    }
}


void bluntify(string gfa_path){
    ifstream file(gfa_path);

    HashGraph gfa_graph;
    IncrementalIdMap<string> id_map;
    OverlapMap overlaps;

    gfa_to_handle_graph(gfa_path, gfa_graph, id_map, overlaps);

    // Where all the ACs go
    vector<AdjacencyComponent> adjacency_components;

    // Compute Adjacency Components and store in vector
    compute_all_adjacency_components(gfa_graph, adjacency_components);

    // Where all the Bicliques go (once we have these, no longer need Adjacency Components)
    vector <vector <edge_t> > bicliques;
    mutex biclique_mutex;

    auto size = gfa_graph.get_node_count() + 1;
    vector <vector <pair <size_t, size_t> > > node_to_biclique_index(size);

    std::cout << "Total adjacency components:\t" << adjacency_components.size() << '\n' << '\n';

    for (size_t i = 0; i<adjacency_components.size(); i++){
//        {
//            cout << "Component " << i << " of size " << adjacency_components[i].size() << '\n' << std::flush;
//            cout << "NODES IN ADJACENCY COMPONENT:\n";
//            for (auto& handle: adjacency_components[i]) {
//                std::cout << id_map.get_name(gfa_graph.get_id(handle)) << (gfa_graph.get_is_reverse(handle) ? "-" : "+")
//                          << '\n';
//            }
//            cout << '\n';
//        }

        compute_all_bicliques(i, gfa_graph, overlaps, adjacency_components, bicliques, biclique_mutex);
    }

    map_splice_sites_by_node(gfa_graph, bicliques, node_to_biclique_index);


}


int main(int argc, char **argv){
    string gfa_path;

    if (argc == 1){
        throw runtime_error("No input gfa path provided");
    }
    else if (argc == 2){
        gfa_path = argv[1];
    }
    else{
        throw runtime_error("Too many arguments. Specify 1 input gfa path.");
    }

    bluntify(gfa_path);

    return 0;
}

