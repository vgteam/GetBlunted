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
using bluntifier::for_each_adjacency_component;
using bluntifier::AdjacencyComponent;
using bluntifier::BipartiteGraph;
using bluntifier::BicliqueCover;
using bluntifier::copy_path_handle_graph;

using handlegraph::handle_t;
using bdsg::HashGraph;


void bluntify(string gfa_path){
    ifstream file(gfa_path);

    HashGraph gfa_graph;
    IncrementalIdMap<string> id_map;
    OverlapMap overlaps;

    gfa_to_handle_graph(gfa_path, gfa_graph, id_map, overlaps);

    // Vector of all the pileups for each adjacency component, with a 1:1 ratio of pileup:biclique
    vector <vector <PoaPileup> > subgraphs;

    for_each_adjacency_component(gfa_graph, [&](AdjacencyComponent& adjacency_component){
        if (adjacency_component.size() == 1) {
            return true;
        }

        subgraphs.emplace_back();

        adjacency_component.decompose_into_bipartite_blocks([&](const BipartiteGraph& bipartite_graph){
            
            auto biclique_cover = BicliqueCover(bipartite_graph).get();
            
            // sort the bicliques in descending order by size (to get any repeated edges
            // into larger POAs -- likely to be more compact this way)
            sort(biclique_cover.begin(). biclique_cover.end(),
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
                
                subgraphs.back().emplace_back();
                
                // Keep track of the number of pileups, used to make unique names for paths
                subgraphs.back().back().index = pileup_index++;
                
                // Iterate all alignments and build a set of alleles for each coordinate
                PoaPileup pileup;
                PileupGenerator::generate_spoa_graph_from_edges(
                      new_edges,
                      id_map,
                      overlaps,
                      gfa_graph,
                      subgraphs.back().back());
            }

            // Add each subgraph to the GFA graph (as an island, initially)
            for (auto& pileup: subgraphs.back()) {
                copy_path_handle_graph(&pileup.graph, &gfa_graph);
            }

            // Iterate LEFT nodes and splice
            for (auto left_iter = bipartite_graph.left_begin(); left_iter != bipartite_graph.left_end(); ++left_iter) {
                for (auto& pileup: subgraphs.back()){

//                    int64_t id;
//                    if (pileup.id_map[0].exists(*left_iter)) {
//                        id = pileup.id_map[0].get_id(*left_iter);
//                    }
//                    else if (pileup.id_map[0].exists(gfa_graph.flip(*left_iter))){
//                        id = pileup.id_map[0].get_id(gfa_graph.flip(*left_iter));
//                    }
//                    else{
//                        throw runtime_error("ERROR: could not find gfa node in pileup data");
//                    }

//                    for (auto& alignment_data: pileup.alignment_data[0][id]){
//
//                    }
                }
            }

            // Iterate RIGHT nodes and splice
            for (auto right_iter = bipartite_graph.right_begin(); right_iter != bipartite_graph.right_end(); ++right_iter) {
                for (auto& pileup: subgraphs.back()){
//                    int64_t id;
//                    if (pileup.id_map[0].exists(*right_iter)) {
//                        id = pileup.id_map[0].get_id(*right_iter);
//                    }
//                    else if (pileup.id_map[0].exists(gfa_graph.flip(*right_iter))){
//                        id = pileup.id_map[0].get_id(gfa_graph.flip(*right_iter));
//                    }
//                    else{
//                        throw runtime_error("ERROR: could not find gfa node in pileup data");
//                    }
                }
            }

        });
    });
}


int main(){
    string script_path = __FILE__;
    string project_directory = parent_path(script_path, 3);

    // Get test GFA path
    string relative_gfa_path = "/data/diploid_case_c.gfa";
//    string relative_gfa_path = "/data/unbalanced_bipartition.gfa";
//    string relative_gfa_path = "/data/staggered_overlap.gfa";
//    string relative_gfa_path = "/data/guppy_360_hg002_messy_small.gfa";
//    string relative_gfa_path = "/data/guppy_360_hg002_mess.gfa";
//    string relative_gfa_path = "/data/1q.gfa";

    const string absolute_gfa_path = join_paths(project_directory, relative_gfa_path);

    bluntify(absolute_gfa_path);

    return 0;
}

