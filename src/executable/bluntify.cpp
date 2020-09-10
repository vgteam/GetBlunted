#include "AdjacencyComponent.hpp"
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
using bluntifier::OverlapMap;
using bluntifier::Alignment;
using bluntifier::adjacency_components;
using bluntifier::AdjacencyComponent;
using bluntifier::bipartition;
using bluntifier::ordered_bipartition;
using bluntifier::BipartiteGraph;
using bluntifier::GaloisLattice;

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

    for (auto& adj_comp : adjacency_components(graph)) {
        if (adj_comp.size() == 1) {
            // TODO: The trivial case is already safe to convert overlap to a graph
            continue;
        }

        auto partition = adj_comp.exhaustive_maximum_bipartite_partition();

        BipartiteGraph bigraph(graph, partition);
        GaloisLattice lattice(bigraph);
        vector<bipartition> separator = lattice.biclique_separator();


        // Iterate all alignments and build a set of alleles for each coordinate


        // Construct a new graph containing the correct alleles











//        for (auto& biclique : separator) {
//            for (auto& a: biclique.first){
//                for (auto& b: biclique.second){
//                    edge.first = a;
//                    edge.second = graph.flip(b);
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
//                }
//            }
//        }
    }
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

