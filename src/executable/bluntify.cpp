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
using bluntifier::Pileup;
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
    IncrementalIdMap id_map;
    OverlapMap overlaps;

    pair<size_t,size_t> lengths;
    edge_t edge;
    size_t start;

    gfa_to_handle_graph(gfa_path, graph, id_map, overlaps);

//    for_each_adjacency_component(graph, [&](AdjacencyComponent& adjacency_component) {
//        if (adjacency_component.size() == 1) {
//            return true;
//        }
//
//        adjacency_component.decompose_into_bipartite_blocks([&](const BipartiteGraph& bipartite_graph) {
//            for (auto a = bipartite_graph.left_begin(); a != bipartite_graph.left_end(); ++a) {
//                for (auto b = bipartite_graph.right_begin(); b != bipartite_graph.right_end(); ++b) {
//                    edge.first = *a;
//                    edge.second = graph.flip(*b);
//
//                    auto iter = overlaps.canonicalize_and_find(edge, graph);
//
//                    cout << id_map.get_name(graph.get_id(edge.first)) << "->"
//                         << id_map.get_name(graph.get_id(edge.second)) << '\n';
//
//                    iter->second.compute_lengths(lengths);
//
//                    start = graph.get_length(edge.first) - lengths.first;
//
//                    cout << lengths.first << " " << lengths.second << '\n';
//                    cout << start << " " << graph.get_length(edge.first) << " " << lengths.first << '\n';
//                    cout << iter->second.create_formatted_alignment_string(graph, edge, start, 0) << '\n';
//
//                }
//            }
//        });
//    });


    for_each_adjacency_component(graph, [&](AdjacencyComponent& adjacency_component){
        if (adjacency_component.size() == 1) {
            return true;
        }

        adjacency_component.decompose_into_bipartite_blocks([&](const BipartiteGraph& bipartite_graph){
            // Iterate all alignments and build a set of alleles for each coordinate
            Pileup pileup;
            BicliqueIterator iterator;

            while (PileupGenerator::traverse_bipartition(graph, overlaps, bipartite_graph, iterator)){
                cout << id_map.get_name(graph.get_id(iterator.node)) << '\n';
                cout << graph.get_sequence(iterator.node) << '\n';
            }

            // Construct a new graph containing the correct alleles

        });
    });
}


int main(){
    string script_path = __FILE__;
    string project_directory = parent_path(script_path, 3);

    // Get test GFA path
    string relative_gfa_path = "/data/staggered_overlap.gfa";
    const string absolute_gfa_path = join_paths(project_directory, relative_gfa_path);

    bluntify(absolute_gfa_path);

    return 0;
}

