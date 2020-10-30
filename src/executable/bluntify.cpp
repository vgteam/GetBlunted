#include "AdjacencyComponent.hpp"
#include "BicliqueCover.hpp"
#include "Biclique.hpp"
#include "OverlapMap.hpp"
#include "Duplicator.hpp"
#include "gfa_to_handle.hpp"
#include "utility.hpp"

#include "bdsg/hash_graph.hpp"

#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"
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
using bluntifier::bipartition;
using bluntifier::copy_path_handle_graph;
using bluntifier::duplicate_prefix;
using bluntifier::duplicate_suffix;
using bluntifier::run_command;
using bluntifier::unchop;

using handlegraph::MutablePathDeletableHandleGraph;
using handlegraph::as_integer;
using handlegraph::handle_t;
using bdsg::HashGraph;

using spoa::AlignmentEngine;
using spoa::AlignmentType;
using spoa::Graph;


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


void convert_spoa_to_bdsg(
        const HashGraph& gfa_graph,
        HashGraph& subgraph,
        Graph& spoa_graph,
        array <map <handle_t, uint32_t>, 2>& handle_paths_per_side){

    auto& paths = spoa_graph.sequences();
    unordered_map <uint32_t, handle_t> nodes_created;
    handle_t previous_subgraph_handle;

    gfa_graph.for_each_handle([&](const handle_t& h){
        cout << as_integer(h) << '\n';
    });
    cout << '\n';

    for (size_t side: {0,1}){
        for (auto& item: handle_paths_per_side[side]){
            auto& gfa_handle = item.first;
            auto spoa_id = item.second;

            // This points to the first SPOA node within the path that this sequence aligned to in the SPOA graph
            auto node = paths[spoa_id];

            size_t base_index = 0;

            string subgraph_path_name = to_string(gfa_graph.get_id(gfa_handle)) + "_" + to_string(side);
            auto subgraph_path_handle = subgraph.create_path_handle(subgraph_path_name);

            cout << gfa_graph.get_sequence(gfa_handle) << '\n';

            while (true) {
                // Check if this node has already been copied to the BDSGraph
                auto iter = nodes_created.find(node->id);

                if (iter == nodes_created.end()) {
                    cout << spoa_id << " " << base_index << " ";

                    char base = gfa_graph.get_base(gfa_handle, base_index);

                    cout << base << '\n';

                    auto new_subgraph_handle = subgraph.create_handle(string(1, base));
                    nodes_created.emplace(node->id, new_subgraph_handle);

                    if (base_index > 0) {
                        subgraph.create_edge(previous_subgraph_handle, new_subgraph_handle);
                    }

                    previous_subgraph_handle = new_subgraph_handle;
                }
                else{
                    if (base_index > 0) {
                        subgraph.create_edge(previous_subgraph_handle, iter->second);
                    }

                    previous_subgraph_handle = iter->second;
                }

                subgraph.append_step(subgraph_path_handle, previous_subgraph_handle);

                base_index++;

                // Check if the spoa path has ended
                if (!(node = node->Successor(spoa_id))) {
                    break;
                }
            }
        }
    }
}


void add_alignments_to_poa(
        const HashGraph& gfa_graph,
        Graph& spoa_graph,
        unique_ptr<AlignmentEngine>& alignment_engine,
        array <map <handle_t, uint32_t>, 2>& handle_paths_per_side,
        const vector<edge_t>& biclique){

    handle_paths_per_side[0].clear();
    handle_paths_per_side[1].clear();

    // If the graph already has some sequences in it, then start the id at that number
    uint32_t spoa_id = spoa_graph.sequences().size();

    for (auto& edge: biclique){
        cout << gfa_graph.get_id(edge.first) << "->" << gfa_graph.get_id(edge.second) << " handle = " << as_integer(edge.first) << "->" << as_integer(edge.second) << '\n';
        if (handle_paths_per_side[0].count(edge.first) == 0){
            cout << gfa_graph.get_id(edge.first) << " L=" << gfa_graph.get_length(edge.first) << " H=" << as_integer(edge.first) << '\n';

            handle_paths_per_side[0].emplace(edge.first, spoa_id++);
            auto sequence = gfa_graph.get_sequence(edge.first);

            auto alignment = alignment_engine->Align(sequence, spoa_graph);
            spoa_graph.AddAlignment(alignment, sequence);
        }
        if (handle_paths_per_side[1].count(edge.second) == 0){
            cout << gfa_graph.get_id(edge.second) << " L=" << gfa_graph.get_length(edge.second) << " H=" << as_integer(edge.second) << '\n';

            handle_paths_per_side[1].emplace(edge.second, spoa_id++);

            auto sequence = gfa_graph.get_sequence(edge.second);

            auto alignment = alignment_engine->Align(sequence, spoa_graph);
            spoa_graph.AddAlignment(alignment, sequence);
        }
    }
}


void align_biclique_overlaps(
        size_t i,
        const HashGraph& gfa_graph,
        const Bicliques& bicliques){

    // TODO: switch to fetch_add atomic
    const auto& biclique = bicliques[i];

    if (biclique.empty()){
        return;
    }

    array <map <handle_t, uint32_t>, 2> handle_paths_per_side;

    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 5, -3, -3, -1);

    spoa::Graph spoa_graph{};

    add_alignments_to_poa(gfa_graph, spoa_graph, alignment_engine, handle_paths_per_side, biclique);

    cout << '\n';

    auto consensus = spoa_graph.GenerateConsensus();

//    std::cout << ">Consensus: " << consensus << std::endl;

    spoa::Graph seeded_spoa_graph{};

    auto seeded_alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 6, -2, -4, -1);

    auto alignment = seeded_alignment_engine->Align(consensus, seeded_spoa_graph);
    seeded_spoa_graph.AddAlignment(alignment, consensus);

    // Iterate a second time on alignment, this time with consensus as the seed
    add_alignments_to_poa(gfa_graph, seeded_spoa_graph, alignment_engine, handle_paths_per_side, biclique);

    {
        auto seeded_consensus = seeded_spoa_graph.GenerateConsensus();

        string prefix = "spoa_overlap_" + to_string(i);
        seeded_spoa_graph.PrintDot(prefix + ".dot");

        string command = "dot -Tpng " + prefix + ".dot -o " + prefix + ".png";
        run_command(command);

        auto msa = seeded_spoa_graph.GenerateMultipleSequenceAlignment();

        for (const auto& it : msa) {
            std::cout << it << std::endl;
        }
        std::cout << '\n';
    }

    HashGraph subgraph;

    convert_spoa_to_bdsg(gfa_graph,subgraph,seeded_spoa_graph,handle_paths_per_side);

    if (subgraph.get_node_count() < 200){
        string test_path_prefix = "test_bluntify_subgraph_" + std::to_string(i);
        handle_graph_to_gfa(subgraph, test_path_prefix + ".gfa");
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    unchop(&subgraph);

    if (subgraph.get_node_count() < 200){
        string test_path_prefix = "test_bluntify_subgraph_unchopped_" + std::to_string(i);
        handle_graph_to_gfa(subgraph, test_path_prefix + ".gfa");
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

}


void bluntify(string gfa_path){
    ifstream file(gfa_path);

    HashGraph gfa_graph;
    IncrementalIdMap<string> id_map;
    OverlapMap overlaps;

    gfa_to_handle_graph(gfa_path, gfa_graph, id_map, overlaps);

    {
        size_t id = 1;
        for (auto& item: id_map.names) {
            cout << id++ << " " << item << '\n';
        }
    }

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

    // TODO: delete adjacency components vector if unneeded

    map_splice_sites_by_node(gfa_graph, bicliques, node_to_biclique_edge);

    Duplicator super_duper(
            node_to_biclique_edge,
            bicliques,
            overlaps);

    if (gfa_graph.get_node_count() < 30){
        string test_path_prefix = "test_bluntify_" + std::to_string(0);
        handle_graph_to_gfa(gfa_graph, test_path_prefix + ".gfa");
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    super_duper.duplicate_all_node_termini(gfa_graph);

    // TODO: keep provenance map for node ID, or else maintain ID map while editing?

    if (gfa_graph.get_node_count() < 30){
        string test_path_prefix = "test_bluntify_" + std::to_string(1);
        handle_graph_to_gfa(gfa_graph, test_path_prefix + ".gfa");
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    for (size_t i=0; i<bicliques.size(); i++){
        align_biclique_overlaps(i, gfa_graph, bicliques);
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

