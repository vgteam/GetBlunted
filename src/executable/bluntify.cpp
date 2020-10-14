#include "AdjacencyComponent.hpp"
#include "PileupGenerator.hpp"
#include "BicliqueCover.hpp"
#include "OverlapMap.hpp"
#include "Splicer.hpp"
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
using bluntifier::BicliqueIterator;
using bluntifier::PoaPileupGenerator;
using bluntifier::PoaPileup;
using bluntifier::OverlapMap;
using bluntifier::Alignment;
using bluntifier::SpliceData;
using bluntifier::Splicer;
using bluntifier::for_each_adjacency_component;
using bluntifier::AdjacencyComponent;
using bluntifier::BipartiteGraph;
using bluntifier::BicliqueCover;
using bluntifier::copy_path_handle_graph;

using handlegraph::handle_t;
using bdsg::HashGraph;



// TODO: Indexes are passed by reference for future threading (will need to be atomic)
void process_adjacency_component(
        size_t& i,
        size_t& biclique_index,
        HashGraph& gfa_graph,
        IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps,
        vector<AdjacencyComponent>& adjacency_components,
        vector <vector <PoaPileup> >& subgraphs,
        vector <vector <SpliceData> >& splice_sites,
        vector <mutex>& splice_site_mutexes){

    auto& adjacency_component = adjacency_components[i];

    cout << "Component " << i << " of size " << adjacency_component.size() << '\n' << std::flush;

    cout << "NODES IN ADJACENCY COMPONENT:\n";
    for (auto& handle: adjacency_component){
        std::cout << id_map.get_name(gfa_graph.get_id(handle)) << (gfa_graph.get_is_reverse(handle) ? "-" : "+") << '\n';
    }
    cout << '\n';

    if (adjacency_component.size() == 1) {
        return;
    }

    adjacency_component.decompose_into_bipartite_blocks([&](const BipartiteGraph& bipartite_graph){

        auto biclique_cover = BicliqueCover(bipartite_graph).get();

        // sort the bicliques in descending order by size (to get any repeated edges
        // into larger POAs -- likely to be more compact this way)
        sort(biclique_cover.begin(), biclique_cover.end(),
             [&](const bipartition& a, const bipartition& b) {
                 return a.first.size() * a.second.size() > b.first.size() * b.second.size();
             });

        unordered_set<edge_t> edges_processed;
        for (const bipartition& biclique : biclique_cover) {
            // get the edges that haven't been handled in a previous biclique
            vector<edge_t> new_edges;
            for (handle_t left : biclique.first) {
                for (handle_t right : biclique.second) {
                    edge_t edge(left, gfa_graph.flip(right));
                    if (!edges_processed.count(edge)) {
                        new_edges.push_back(edge);
                    }
                }
            }

            // Keep track of the number of pileups, used to make unique names for paths
            subgraphs[i].emplace_back();
            subgraphs[i].back().biclique_index = biclique_index++;

            // Iterate all alignments and build a set of splice sites and pileups
            PoaPileup pileup;
            PoaPileupGenerator pileup_generator(
                    bipartite_graph,
                    id_map,
                    overlaps,
                    gfa_graph,
                    pileup,
                    splice_sites,
                    splice_site_mutexes);

            pileup_generator.generate_from_edges(new_edges,i);
        }
    });
}



void bluntify(string gfa_path){
    ifstream file(gfa_path);

    HashGraph gfa_graph;
    IncrementalIdMap<string> id_map;
    OverlapMap overlaps;

    gfa_to_handle_graph(gfa_path, gfa_graph, id_map, overlaps);

    // Where all the ACs go
    vector<AdjacencyComponent> adjacency_components;

    // Compute adjacency components and store in vector
    compute_adjacency_components(gfa_graph, adjacency_components);

    // Vector of all the pileups for each adjacency component, with a 1:1 ratio of pileup:biclique
    vector <vector <PoaPileup> > subgraphs(adjacency_components.size());

    // Node ids are 1-based so the first entry is kept empty
    auto size = gfa_graph.get_node_count() + 1;

    // Vector to store all the splice sites for each node.
    vector <vector <SpliceData> > splice_sites(size);
    vector <mutex> splice_site_mutexes(size);

    std::cout << "Total adjacency components:\t" << adjacency_components.size() << '\n';

    gfa_graph.for_each_handle([&](const handle_t& h){
        cout << gfa_graph.get_is_reverse(h) << '\n';
    });
    cout << '\n';
    gfa_graph.for_each_edge([&](const edge_t& e){
        auto a = e;
        auto iter = overlaps.canonicalize_and_find(a, gfa_graph);
        cout << gfa_graph.get_is_reverse(iter->first.first) << '\n';
        cout << gfa_graph.get_is_reverse(iter->first.second) << '\n' << '\n';
    });

    // TODO: thread this function
    size_t biclique_index = 0;
    for (size_t i = 0; i<adjacency_components.size(); i++){
        process_adjacency_component(
                i,
                biclique_index,
                gfa_graph,
                id_map,
                overlaps,
                adjacency_components,
                subgraphs,
                splice_sites,
                splice_site_mutexes);
    }

    for (size_t node_id=1; node_id<splice_sites.size(); node_id++){
        cout << id_map.get_name(node_id) << '\n';
        auto h = gfa_graph.get_handle(node_id, false);
        size_t node_length = gfa_graph.get_length(h);

        for (auto& site: splice_sites[node_id]){
            cout << site.get_forward_coordinate(node_length) << " "
                 << site.forward_splice_is_left() << " "
                 << site.biclique_index << '\n';
        }
        cout << '\n';

        Splicer splicer(gfa_graph, splice_sites, node_id);

        // TODO: Actually do something with the overlapping overlaps
        splicer.find_overlapping_overlaps();

        // Find the middlemost splice site for each biclique
        splicer.duplicate_all_termini();

        // Split the node and duplicate prefixes/suffixes

        // For each prefix/suffix, perform a similar splitting process

//        auto h = gfa_graph.get_handle(node_id, false);
//        gfa_graph.divide_handle(h, all_sites);
    }

    // Add each subgraph to the GFA graph (as an island, initially)
    for (auto& pileup: subgraphs.back()) {
        copy_path_handle_graph(&pileup.graph, &gfa_graph);
    }
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

