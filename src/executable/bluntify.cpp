#include "AdjacencyComponent.hpp"
#include "BicliqueCover.hpp"
#include "Biclique.hpp"
#include "OverlapMap.hpp"
#include "Duplicator.hpp"
#include "Subgraph.hpp"
#include "align.hpp"
#include "OverlappingOverlap.hpp"
#include "OverlappingOverlapSplicer.hpp"
#include "gfa_to_handle.hpp"
#include "utility.hpp"

#include "bdsg/hash_graph.hpp"

#include "unchop.hpp"

using bluntifier::gfa_to_path_handle_graph_in_memory;
using bluntifier::gfa_to_path_handle_graph;
using bluntifier::gfa_to_handle_graph;
using bluntifier::handle_graph_to_gfa;
using bluntifier::parent_path;
using bluntifier::join_paths;
using bluntifier::IncrementalIdMap;
using bluntifier::OverlapMap;
using bluntifier::Duplicator;
using bluntifier::Alignment;
using bluntifier::for_each_adjacency_component;
using bluntifier::AdjacencyComponent;
using bluntifier::BipartiteGraph;
using bluntifier::BicliqueCover;
using bluntifier::BicliqueEdgeIndex;
using bluntifier::Bicliques;
using bluntifier::NodeInfo;
using bluntifier::Subgraph;
using bluntifier::OverlappingOverlapSplicer;
using bluntifier::OverlappingNodeInfo;
using bluntifier::OverlappingChild;
using bluntifier::bipartition;
using bluntifier::copy_path_handle_graph;
using bluntifier::duplicate_prefix;
using bluntifier::duplicate_suffix;
using bluntifier::run_command;
using bluntifier::unchop;
using bluntifier::harmonize_biclique_orientations;

using handlegraph::MutablePathDeletableHandleGraph;
using handlegraph::as_integer;
using handlegraph::handle_t;
using bdsg::HashGraph;


void deduplicate_and_canonicalize_biclique_cover(
        vector <bipartition>& biclique_cover,
        vector <vector <edge_t> >& deduplicated_biclique_cover,
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
        deduplicated_biclique_cover.emplace_back();

        // get the edges that haven't been handled in a previous biclique
        for (handle_t left : biclique.first) {
            for (handle_t right : biclique.second) {
                edge_t edge(left, gfa_graph.flip(right));
                auto iter = overlaps.canonicalize_and_find(edge, gfa_graph);

                if (!edges_processed.count(edge)) {
                    edges_processed.emplace(edge);
                    deduplicated_biclique_cover.back().emplace_back(iter->first);
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
        Bicliques& bicliques,
        mutex& biclique_mutex){

    // TODO: switch to fetch_add atomic
    auto& adjacency_component = adjacency_components[i];

    // Skip trivial adjacency components (dead ends)
    if (adjacency_component.size() == 1) {
        return;
    }

    adjacency_component.decompose_into_bipartite_blocks([&](const BipartiteGraph& bipartite_graph){
        vector <bipartition> biclique_cover = BicliqueCover(bipartite_graph).get();
        vector <vector <edge_t> > deduplicated_biclique_cover;

        // TODO: find a lock-minimal thread safe way to prevent copying each biclique cover during duplication
        // TODO: Maybe just move the deduplication step outside of thread fn?
        deduplicate_and_canonicalize_biclique_cover(
                biclique_cover,
                deduplicated_biclique_cover,
                gfa_graph,
                overlaps);

        for (auto& biclique: deduplicated_biclique_cover) {
            biclique_mutex.lock();
            bicliques.bicliques.emplace_back(biclique);
            biclique_mutex.unlock();
        }
    });
}


void map_splice_sites_by_node(
        const HandleGraph& gfa_graph,
        const Bicliques& bicliques,
        vector <vector <BicliqueEdgeIndex> >& node_to_biclique_edge){

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

            node_to_biclique_edge[left_node_id].emplace_back(i,j);

            // Don't make 2 mappings ot he same edge if it is a self-loop
            if (right_node_id != left_node_id) {
                node_to_biclique_edge[right_node_id].emplace_back(i, j);
            }
        }
    }
}


void print_adjacency_components_stats(
        size_t i,
        vector<AdjacencyComponent>& adjacency_components,
        IncrementalIdMap<string>& id_map,
        HandleGraph& gfa_graph){

    cout << "Component " << i << " of size " << adjacency_components[i].size() << '\n' << std::flush;
    cout << "NODES IN ADJACENCY COMPONENT:\n";
    for (auto& handle: adjacency_components[i]) {
        std::cout << id_map.get_name(gfa_graph.get_id(handle)) << (gfa_graph.get_is_reverse(handle) ? "-" : "+")
                  << '\n';
    }
    cout << '\n';
}


void splice_subgraphs(
        HashGraph& gfa_graph,
        vector <Subgraph>& subgraphs,
        map <nid_t, nid_t>& child_to_parent,
        map<nid_t, OverlappingNodeInfo>& overlapping_overlap_nodes){

    size_t i = 0;
    for (auto& subgraph: subgraphs) {

        // First, copy the subgraph into the GFA graph
        subgraph.graph.increment_node_ids(gfa_graph.max_node_id());
        copy_path_handle_graph(&subgraph.graph, &gfa_graph);

        if (gfa_graph.get_node_count() < 30){
            string test_path_prefix = "test_bluntify_splice_" + std::to_string(i) + "_b";
            handle_graph_to_gfa(gfa_graph, test_path_prefix + ".gfa");
            string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                             + test_path_prefix + ".png";
            run_command(command);
        }

        i++;

        // Iterate the suffixes/prefixes that participated in this biclique
        for (bool side: {0, 1}) {
            for (auto& item: subgraph.paths_per_handle[side]) {
                auto& handle = item.first;
                auto& path_info = item.second;
                auto node_id = gfa_graph.get_id(handle);

                cout << "Splicing node " << node_id << ", side " << side << '\n';

                // Check if this is an Overlapping Overlap node
                auto original_gfa_node = child_to_parent[node_id];
                auto result = overlapping_overlap_nodes.find(original_gfa_node);
                if (result != overlapping_overlap_nodes.end()) {
                    // Do re-chopping and extra splicing ?
                    // or just skip for now?
                    continue;
                }

                // Find the path handle for the path that was copied into the GFA graph
                auto path_name = subgraph.graph.get_path_name(path_info.path_handle);
                auto path_handle = gfa_graph.get_path_handle(path_name);

                cout << "\tPath sequence:\t";
                for (const auto& h: gfa_graph.scan_path(path_handle)) {
                    cout << gfa_graph.get_sequence(h);
                }
                cout << '\n';
                cout << "\tNode sequence:\t" << gfa_graph.get_sequence(handle) << '\n';

                set<handle_t> parent_handles;
                gfa_graph.follow_edges(handle, 1 - side, [&](const handle_t& h) {
                    parent_handles.emplace(h);
                });

                if (parent_handles.empty()) {
                    throw runtime_error("ERROR: biclique terminus does not have any parent: " + to_string(node_id));
                }

                for (auto& parent_handle: parent_handles) {
                    // Depending on which side of the biclique this node is on, its path in the POA will be spliced
                    // differently
                    if (path_info.biclique_side == 0) {
                        auto& left = parent_handle;
                        auto right = gfa_graph.get_handle_of_step(gfa_graph.path_begin(path_handle));

                        gfa_graph.create_edge(left, right);
                    }
                    else {
                        auto left = gfa_graph.get_handle_of_step(gfa_graph.path_back(path_handle));
                        auto& right = parent_handle;

                        gfa_graph.create_edge(left, right);
                    }
                }

                cout << "Destroying: " << gfa_graph.get_id(handle) << '\n';

                if (subgraph.paths_per_handle[1-side].count(handle) == 0
                    and subgraph.paths_per_handle[1-side].count(gfa_graph.flip(handle)) == 0) {
                    gfa_graph.destroy_handle(handle);
                }
            }
        }
    }
}


void bluntify(string gfa_path){
    ifstream file(gfa_path);

    HashGraph gfa_graph;
    IncrementalIdMap<string> id_map;
    OverlapMap overlaps;

    gfa_to_handle_graph(gfa_path, gfa_graph, id_map, overlaps);

//    {
//        size_t id = 1;
//        for (auto& item: id_map.names) {
//            cout << id++ << " " << item << '\n';
//        }
//    }

    // Where all the ACs go
    vector<AdjacencyComponent> adjacency_components;

    // Compute Adjacency Components and store in vector
    compute_all_adjacency_components(gfa_graph, adjacency_components);

    // Where all the Bicliques go (once we have these, no longer need Adjacency Components)
    Bicliques bicliques;
    mutex biclique_mutex;

    auto size = gfa_graph.get_node_count() + 1;
    vector <vector <BicliqueEdgeIndex> > node_to_biclique_edge(size);

    std::cout << "Total adjacency components:\t" << adjacency_components.size() << '\n' << '\n';

    for (size_t i = 0; i<adjacency_components.size(); i++){
        print_adjacency_components_stats(i,adjacency_components,id_map,gfa_graph);
        compute_all_bicliques(i, gfa_graph, overlaps, adjacency_components, bicliques, biclique_mutex);
    }

    {
        size_t i = 0;
        for (auto& biclique: bicliques.bicliques) {
            cout << "Biclique " << i++ << '\n';
            for (auto& edge: biclique) {
                cout << "(" << gfa_graph.get_id(edge.first);
                cout << (gfa_graph.get_is_reverse(edge.first) ? "-" : "+");
                cout << ") -> (" << gfa_graph.get_id(edge.second);
                cout << (gfa_graph.get_is_reverse(edge.second) ? "-" : "+") << ")" << '\n';
            }
        }
        cout << '\n' << '\n';
    }

    // TODO: delete adjacency components vector if unneeded

    map_splice_sites_by_node(gfa_graph, bicliques, node_to_biclique_edge);

    Duplicator super_duper(node_to_biclique_edge,bicliques,overlaps);

    if (gfa_graph.get_node_count() < 30){
        string test_path_prefix = "test_bluntify_" + std::to_string(0);
        handle_graph_to_gfa(gfa_graph, test_path_prefix + ".gfa");
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    super_duper.duplicate_all_node_termini(gfa_graph);

    if (gfa_graph.get_node_count() < 30){
        string test_path_prefix = "test_bluntify_" + std::to_string(1);
        handle_graph_to_gfa(gfa_graph, test_path_prefix + ".gfa");
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    harmonize_biclique_orientations(gfa_graph, bicliques);

    vector <Subgraph> biclique_subgraphs(bicliques.size());

    for (size_t i=0; i<bicliques.size(); i++){
        // TODO: finish refactor?
        align_biclique_overlaps(i, gfa_graph, bicliques, biclique_subgraphs);
    }

//    for (auto& item: super_duper.parent_to_children){
//        cout << "Children of node " << item.first << " with original name: " << id_map.get_name(item.first) << '\n';
//        for (auto& child: item.second){
//            cout << child << " " << gfa_graph.get_sequence(gfa_graph.get_handle(child, false)) << '\n';
//        }
//    }
//    cout << '\n';

    splice_subgraphs(gfa_graph, biclique_subgraphs, super_duper.child_to_parent, super_duper.overlapping_overlap_nodes);

    if (gfa_graph.get_node_count() < 30){
        string test_path_prefix = "test_bluntify_spliced_" + std::to_string(1);
        handle_graph_to_gfa(gfa_graph, test_path_prefix + ".gfa");
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    OverlappingOverlapSplicer oo_splicer(super_duper.overlapping_overlap_nodes, biclique_subgraphs);

    oo_splicer.splice_overlapping_overlaps(gfa_graph);
//    splice_overlapping_overlaps(gfa_graph, biclique_subgraphs, super_duper.overlapping_overlap_nodes);
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

