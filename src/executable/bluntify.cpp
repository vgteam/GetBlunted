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


class BicliqueEdgeIndex{
public:
    size_t biclique_index;
    size_t edge_index;

    BicliqueEdgeIndex(size_t biclique, size_t edge);
};


BicliqueEdgeIndex::BicliqueEdgeIndex(size_t biclique, size_t edge):
    biclique_index(biclique),
    edge_index(edge)
{}


class Bicliques{
public:
    vector <vector <edge_t> > bicliques;

    edge_t& operator[](BicliqueEdgeIndex i);
    const edge_t& operator[](BicliqueEdgeIndex i) const;

    vector <edge_t>& operator[](size_t i);
    const vector <edge_t>& operator[](size_t i) const;

    size_t size() const;
};


size_t Bicliques::size() const{
    return bicliques.size();
}


edge_t& Bicliques::operator[](BicliqueEdgeIndex i){
    return bicliques[i.biclique_index][i.edge_index];
}


const edge_t& Bicliques::operator[](BicliqueEdgeIndex i) const{
    return bicliques[i.biclique_index][i.edge_index];
}


vector <edge_t>& Bicliques::operator[](size_t i){
    return bicliques[i];
}


const vector <edge_t>& Bicliques::operator[](size_t i) const{
    return bicliques[i];
}


class OverlapInfo{
public:
    size_t edge_index;
    size_t length;

    OverlapInfo(size_t edge_index, size_t length);
};


OverlapInfo::OverlapInfo(size_t edge_index, size_t length) :
    edge_index(edge_index),
    length(length)
{}


class NodeInfo{
public:
    array <map <size_t, vector <OverlapInfo> >, 2> factored_overlaps;
    const vector <vector <BicliqueEdgeIndex> >& node_to_biclique_edge;
    const Bicliques& bicliques;
    const HandleGraph& gfa_graph;
    const OverlapMap& overlaps;
    const size_t node_id;

    NodeInfo(
            const vector <vector <BicliqueEdgeIndex> >& node_to_biclique_edge,
            const Bicliques& bicliques,
            const HandleGraph& gfa_graph,
            const OverlapMap& overlaps,
            size_t node_id);

    void factor_overlaps_by_biclique_and_side();
    void sort_factored_overlaps();

    size_t get_overlap_length(edge_t edge, bool side);
    void print_stats();
};


NodeInfo::NodeInfo(
        const vector <vector <BicliqueEdgeIndex> >& node_to_biclique_edge,
        const Bicliques& bicliques,
        const HandleGraph& gfa_graph,
        const OverlapMap& overlaps,
        size_t node_id):
    node_to_biclique_edge(node_to_biclique_edge),
    bicliques(bicliques),
    gfa_graph(gfa_graph),
    overlaps(overlaps),
    node_id(node_id)
{
    factor_overlaps_by_biclique_and_side();
    sort_factored_overlaps();
}


void NodeInfo::print_stats() {
    cout << "Node " << node_id << '\n';

    for (size_t side: {0,1}){
        auto biclique_overlaps = factored_overlaps[side];

        cout << "  Side " << side << '\n';

        for (const auto& item: biclique_overlaps){
            auto biclique_index = item.first;
            auto overlap_infos = item.second;

            cout << "    Biclique " << biclique_index << '\n';

            for (const auto& overlap_info: overlap_infos){
                cout << "      " << overlap_info.edge_index << " " << overlap_info.length << '\n';
            }
        }

        cout << '\n';
    }

    cout << '\n';
}


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
        Bicliques& bicliques,
        mutex& biclique_mutex){

    auto& adjacency_component = adjacency_components[i];

    // Skip trivial adjacency components (dead ends)
    if (adjacency_component.size() == 1) {
        return;
    }

    adjacency_component.decompose_into_bipartite_blocks([&](const BipartiteGraph& bipartite_graph){
        biclique_mutex.lock();
        bicliques.bicliques.emplace_back();
        biclique_mutex.unlock();

        vector<bipartition> biclique_cover = BicliqueCover(bipartite_graph).get();

        deduplicate_and_canonicalize_edges_in_bicliques(biclique_cover, bicliques.bicliques.back(), gfa_graph, overlaps);
    });
}


size_t NodeInfo::get_overlap_length(edge_t edge, bool side){
    pair<size_t, size_t> lengths;
    overlaps.at(edge)->second.compute_lengths(lengths);

    size_t length;
    if (side == 0){
        length = lengths.first;
    }
    else{
        length = lengths.second;
    }

    return length;
}


// For one node, make a mapping: (side -> (biclique_index -> (edge_index,length) ) )
void NodeInfo::factor_overlaps_by_biclique_and_side() {

    for (auto& index: node_to_biclique_edge[node_id]) {
        auto edge = bicliques[index];

        auto left_node_id = gfa_graph.get_id(edge.first);
        auto right_node_id = gfa_graph.get_id(edge.second);

        // It's possible that the edge is a self-edge. Add the edge (index) to any side that it matches.
        if (left_node_id == nid_t(node_id)) {
            auto length = get_overlap_length(edge, 0);
            factored_overlaps[0][index.biclique_index].emplace_back(index.edge_index, length);
        }
        if (right_node_id == nid_t(node_id)) {
            auto length = get_overlap_length(edge, 1);
            factored_overlaps[1][index.biclique_index].emplace_back(index.edge_index, length);
        }
    }
}


void NodeInfo::sort_factored_overlaps(){
    for (bool side: {0,1}) {
        auto biclique_edge_indexes = factored_overlaps[side];

        for (auto& biclique: biclique_edge_indexes) {
            auto& biclique_index = biclique.first;
            auto& overlap_infos = biclique.second;

            sort(overlap_infos.begin(), overlap_infos.end(), [&](OverlapInfo a, OverlapInfo b){
                return a.length > b.length;
            });
        }
    }
}


//void make_babies(HandleGraph& gfa_graph, NodeInfo node_info, bool side){
//    for (){
//
//    }
//}


void duplicate_termini(
        const vector <vector <BicliqueEdgeIndex> >& node_to_biclique_edge,
        const Bicliques& bicliques,
        HandleGraph& gfa_graph,
        OverlapMap& overlaps){

    for (size_t node_id=1; node_id<node_to_biclique_edge.size(); node_id++){
        // Factor the overlaps into hierarchy: side -> biclique -> (overlap, length)
        NodeInfo node_info(node_to_biclique_edge, bicliques, gfa_graph, overlaps, node_id);

    }
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
            node_to_biclique_edge[right_node_id].emplace_back(i,j);
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
    Bicliques bicliques;
    mutex biclique_mutex;

    auto size = gfa_graph.get_node_count() + 1;
    vector <vector <BicliqueEdgeIndex> > node_to_biclique_edge(size);

    std::cout << "Total adjacency components:\t" << adjacency_components.size() << '\n' << '\n';

    for (size_t i = 0; i<adjacency_components.size(); i++){
        print_adjacency_components_stats(i,adjacency_components,id_map,gfa_graph);
        compute_all_bicliques(i, gfa_graph, overlaps, adjacency_components, bicliques, biclique_mutex);
    }

    // TODO: delete adjacency components

    map_splice_sites_by_node(gfa_graph, bicliques, node_to_biclique_edge);

    duplicate_termini(node_to_biclique_edge, bicliques, gfa_graph, overlaps);

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

