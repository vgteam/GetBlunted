#include "Splicer.hpp"
#include <queue>

using std::queue;
using std::tie;


namespace bluntifier {


Splicer::Splicer(MutablePathMutableHandleGraph& gfa_graph, vector<vector<SpliceData> >& splice_sites, const size_t node_id):
    gfa_graph(gfa_graph),
    splice_sites(splice_sites),
    node_id(node_id)
{}


/// Gather all the indexes that the node needs to be split at (initially, for duplication purposes)
/// TODO: adapt for multiple bicliques
void Splicer::find_duplication_sites(
        vector <pair <size_t,size_t> >& left_sites,
        vector <pair <size_t,size_t> >& right_sites){

    // Keep track of the middlemost splice site for every biclique, on each side of this node
    map<size_t,size_t> max_left_sites;
    map<size_t,size_t> min_right_sites;

    auto h = gfa_graph.get_handle(node_id, false);
    size_t node_length = gfa_graph.get_length(h);

    for (size_t i=0; i<splice_sites[node_id].size(); i++) {
        auto& site = splice_sites[node_id][i];

        if (site.forward_splice_is_left()) {
            auto coord = site.get_forward_coordinate(node_length);

            auto iter = max_left_sites.find(site.biclique_index);

            if (iter != max_left_sites.end()){
                auto prev_max_index = iter->second;
                auto& prev_max_site = splice_sites[node_id][prev_max_index];
                auto prev_max_coord = prev_max_site.get_forward_coordinate(node_length);

                if (coord > prev_max_coord){
                    iter->second = i;
                }
            }
            else{
                max_left_sites.emplace(site.biclique_index, i);
            }
        }
        else{
            auto coord = site.get_forward_coordinate(node_length);

            auto iter = min_right_sites.find(site.biclique_index);

            if (iter != min_right_sites.end()){
                auto prev_min_index = iter->second;
                auto& prev_min_site = splice_sites[node_id][prev_min_index];
                auto prev_min_coord = prev_min_site.get_forward_coordinate(node_length);

                if (coord < prev_min_coord){
                    iter->second = i;
                }
            }
            else{
                min_right_sites.emplace(site.biclique_index, i);
            }
        }
    }

    cout << "DUPLICATION SITES\n";
    cout << "LEFT sites\n";
    for (auto& item: max_left_sites){
        auto coord = splice_sites[node_id][item.second].get_forward_coordinate(node_length);
        cout << item.first << " " << coord << '\n';
        left_sites.emplace_back(item);
    }
    cout << "RIGHT sites\n";
    for (auto& item: min_right_sites){
        auto coord = splice_sites[node_id][item.second].get_forward_coordinate(node_length);
        cout << item.first << " " << coord << '\n';
        right_sites.emplace_back(item);
    }
    cout << '\n';

    // Return objects should be a sorted vector, not a map, because maps are unweildy and no one likes them
    sort(left_sites.begin(), left_sites.end(), [&](const pair<size_t,size_t>& a, const pair<size_t,size_t>& b){
        auto& a_site = splice_sites[node_id][a.second];
        auto a_value = a_site.get_forward_coordinate(node_length);

        auto& b_site = splice_sites[node_id][b.second];
        auto b_value = b_site.get_forward_coordinate(node_length);

        return a_value < b_value;
    });

    sort(right_sites.begin(), right_sites.end(), [&](const pair<size_t,size_t>& a, const pair<size_t,size_t>& b){
        auto& a_site = splice_sites[node_id][a.second];
        auto a_value = a_site.get_forward_coordinate(node_length);

        auto& b_site = splice_sites[node_id][b.second];
        auto b_value = b_site.get_forward_coordinate(node_length);

        return a_value < b_value;
    });
}


void Splicer::find_overlapping_overlaps(){
    vector<size_t> indexes;

    auto h = gfa_graph.get_handle(node_id, false);
    size_t node_length = gfa_graph.get_length(h);

    for (size_t i=0; i<splice_sites[node_id].size(); i++){
        indexes.push_back(i);
    }

    // Sort the indexes instead of the array itself
    sort(indexes.begin(), indexes.end(), [&](const size_t& a, const size_t& b){
        auto a_value = splice_sites[node_id][a].get_forward_coordinate(node_length);
        auto b_value = splice_sites[node_id][b].get_forward_coordinate(node_length);

        if (a_value == b_value){
            return splice_sites[node_id][a].forward_splice_is_left();
        }
        return a_value < b_value;
    });

    bool right_visited = false;
    queue <size_t> right_visited_queue;
    bool is_left;

    // In the sorted splice sites, find every case where a left site is after a right site or vice versa.
    // These cases are the overlapping overlaps.
    for (auto& i: indexes){
        is_left = splice_sites[node_id][i].forward_splice_is_left();

        if (is_left){
            while (not right_visited_queue.empty()){
                auto i_queue = right_visited_queue.front();
                right_visited_queue.pop();

                auto coord = splice_sites[node_id][i_queue].get_forward_coordinate(node_length);
                cout << "overlapping overlap:\tR " << coord << '\n';
            }

            if (right_visited){
                auto coord = splice_sites[node_id][i].get_forward_coordinate(node_length);
                cout << "overlapping overlap:\tL " << coord << '\n';
            }
        }
        else {
            right_visited_queue.push(i);
            right_visited = true;
        }
    }
}


void Splicer::duplicate_terminus(size_t site_index){
    auto& site = splice_sites[node_id][site_index];

    auto h = gfa_graph.get_handle(node_id, false);
    size_t node_length = gfa_graph.get_length(h);

    auto coord = site.get_forward_coordinate(node_length);

    handle_t left_handle;
    handle_t right_handle;
//    tie (left_handle, right_handle) = gfa_graph.divide_handle(h, coord);

    cout << coord << '\n';

    // Update all the splice sites that just got broken
    for (auto other_site: splice_sites[node_id]){
        auto other_coord = other_site.get_forward_coordinate(node_length);

        if (other_coord >= coord and not other_site.forward_splice_is_left()){
            cout << other_site.get_forward_coordinate(node_length) << '\n' << '\n';
        }
    }
}


void Splicer::duplicate_all_termini(){
    // Pairs of (biclique_index, site_index)
    vector <pair <size_t,size_t> > left_sites;
    vector <pair <size_t,size_t> > right_sites;

    find_duplication_sites(left_sites, right_sites);

    for (auto& item: right_sites){
        duplicate_terminus(item.second);
    }
    for (auto& item: left_sites){
        duplicate_terminus(item.second);
    }
}


}