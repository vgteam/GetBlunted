#include "AdjacencyComponent.hpp"
#include "PileupGenerator.hpp"
#include "BicliqueCover.hpp"
#include "OverlapMap.hpp"
#include "gfa_to_handle.hpp"
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

using handlegraph::handle_t;
using bdsg::HashGraph;


void bluntify(string gfa_path){
    ifstream file(gfa_path);

    HashGraph graph;
    IncrementalIdMap<string> id_map;
    OverlapMap overlaps;

    gfa_to_handle_graph(gfa_path, graph, id_map, overlaps);

    // Vector of all the pileups for each adjacency component, with a 1:1 ratio of pileup:biclique
    vector <vector <PoaPileup> > subgraphs;

    for_each_adjacency_component(graph, [&](AdjacencyComponent& adjacency_component){
        if (adjacency_component.size() == 1) {
            return true;
        }

        subgraphs.emplace_back();

        adjacency_component.decompose_into_bipartite_blocks([&](const BipartiteGraph& bipartite_graph){
            subgraphs.back().emplace_back();

            // Iterate all alignments and build a set of alleles for each coordinate
            PoaPileup pileup;
            PileupGenerator::generate_spoa_graph_from_bipartition(
                    bipartite_graph,
                    id_map,
                    overlaps,
                    graph,
                    subgraphs.back().back());

            // Construct a new graph containing the correct alleles

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

