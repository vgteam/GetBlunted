#include "AdjacencyComponent.hpp"
#include "PileupGenerator.hpp"
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
using bluntifier::BicliqueIterator;
using bluntifier::PileupGenerator;
using bluntifier::PoaPileup;
using bluntifier::OverlapMap;
using bluntifier::Alignment;
using bluntifier::AlignmentData;
using bluntifier::for_each_adjacency_component;
using bluntifier::AdjacencyComponent;
using bluntifier::BipartiteGraph;
using bluntifier::BicliqueCover;
using bluntifier::copy_path_handle_graph;

using handlegraph::handle_t;
using bdsg::HashGraph;


void process_adjacency_component(
        size_t i,
        HashGraph& gfa_graph,
        IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps,
        vector<AdjacencyComponent>& adjacency_components,
        vector <vector <PoaPileup> >& subgraphs,
        vector <vector <AlignmentData> >& splice_sites,
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
        uint64_t pileup_index = 0;
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
            subgraphs[i].back().index = pileup_index++;

            // Iterate all alignments and build a set of splice sites and pileups
            PoaPileup pileup;
            PileupGenerator::generate_spoa_graph_from_edges(
                    new_edges,
                    id_map,
                    overlaps,
                    gfa_graph,
                    subgraphs[i].back(),
                    splice_sites,
                    splice_site_mutexes,
                    i);
        }
    });

}


/// Gather all the indexes that the node needs to be split at (initially, for duplication purposes)
/// TODO: adapt for multiple bicliques
void aggregate_splice_sites(vector<size_t>& left_sites, vector<size_t>& right_sites, vector<size_t>& all_sites){
    set<size_t> all_nonduplicate_sites;

    all_nonduplicate_sites.emplace(left_sites.back());

    // Because divide handle cuts AFTER the index provided, right sites need to be offset by 1
    all_nonduplicate_sites.emplace(right_sites.front() - 1);

    for (auto& item: all_nonduplicate_sites){
        cout << item << '\n';
        all_sites.emplace_back(item);
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

    // Compute adjacency components and store in vector
    compute_adjacency_components(gfa_graph, adjacency_components);

    // Vector of all the pileups for each adjacency component, with a 1:1 ratio of pileup:biclique
    vector <vector <PoaPileup> > subgraphs(adjacency_components.size());

    // Node ids are 1-based so the first entry is kept empty
    auto size = gfa_graph.get_node_count() + 1;

    // Vector to store all the splice sites for each node.
    vector <vector <AlignmentData> > splice_sites(size);
    vector <mutex> splice_site_mutexes(size);

    std::cout << "Total adjacency components:\t" << adjacency_components.size() << '\n';

    // TODO: thread this function
    for (size_t i = 0; i<adjacency_components.size(); i++){
        process_adjacency_component(
                i,
                gfa_graph,
                id_map,
                overlaps,
                adjacency_components,
                subgraphs,
                splice_sites,
                splice_site_mutexes);
    }

    for (size_t node_id=1; node_id<splice_sites.size(); node_id++){
        vector <size_t> left_sites;
        vector <size_t> right_sites;
        vector <size_t> all_sites;

        for (auto& site: splice_sites[node_id]){
//            cout << site << '\n' << '\n';

            // TODO split on a biclique basis and duplicate starting at the middlemost index for each biclique
            auto h = gfa_graph.get_handle(node_id, site.is_reverse);

//          If the AlignmentData tells us that it is "left", that means the node is on the left of an overlap,
//          which indicates that the overlap is happening on the right end of the node...
            if (site.is_left){
                if (not site.is_reverse){
                    right_sites.emplace_back(site.sequence_start_index);
                }
                else{
                    size_t forward_start_index = gfa_graph.get_length(h) - site.sequence_start_index - 1;
                    left_sites.emplace_back(forward_start_index);
                }
            }
            else{
                if (not site.is_reverse){
                    left_sites.emplace_back(site.sequence_stop_index);
                }
                else{
                    size_t forward_start_index = gfa_graph.get_length(h) - site.sequence_stop_index - 1;
                    right_sites.emplace_back(forward_start_index);
                }
            }
        }

        sort(left_sites.begin(), left_sites.end());
        sort(right_sites.begin(), right_sites.end());

        cout << "Splice sites for node: " << id_map.get_name(node_id) << '\n';

        // Find overlapping overlaps
        if (not right_sites.empty() and not left_sites.empty()) {
            for (auto& item: left_sites) {
                if (item >= right_sites[0]) {
                    std::cout << "overlapping overlap!\n";
                }
            }
        }

        cout << "Left sites:\n";
        for (auto& item: left_sites){
            cout << item << ',';
        }
        cout << '\n';
        cout << "Right sites:\n";
        for (auto& item: right_sites){
            cout << item << ',';
        }
        cout << '\n' << '\n';

        aggregate_splice_sites(left_sites, right_sites, all_sites);
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
        throw runtime_error("No input gfa path found");
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

