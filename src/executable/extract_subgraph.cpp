#include <bdsg/hash_graph.hpp>
#include "IncrementalIdMap.hpp"
#include "handle_to_gfa.hpp"
#include "gfa_to_handle.hpp"
#include "OverlapMap.hpp"

#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <queue>
#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>

using bdsg::HashGraph;
using handlegraph::handle_t;
using handlegraph::nid_t;

using std::unordered_set;
using std::unordered_map;
using std::function;
using std::ofstream;
using std::queue;
using std::string;
using std::vector;
using std::runtime_error;


#include "Filesystem.hpp"
#include "utility.hpp"
#include "CLI11.hpp"

using ghc::filesystem::path;
using bluntifier::write_bfs_to_gfa;
using bluntifier::IncrementalIdMap;
using bluntifier::OverlapMap;


int main(int argc, char* argv[]){
    path gfa_path;
    size_t radius;
    string start_name;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input_gfa",
            gfa_path,
            "Path to GFA containing phased non-overlapping segments")
            ->required();

    app.add_option(
            "-r,--radius",
            radius,
            "Distance from start node to traverse")
            ->required();

    app.add_option(
            "-s,--start",
            start_name,
            "Name of node/segment to start at")
            ->required();

    CLI11_PARSE(app, argc, argv);

    bluntifier::write_bfs_to_gfa(gfa_path, start_name, radius);


    return 0;
}
