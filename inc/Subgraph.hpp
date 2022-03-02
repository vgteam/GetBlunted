#ifndef BLUNTIFIER_SUBGRAPH_HPP
#define BLUNTIFIER_SUBGRAPH_HPP

#include "handlegraph/handle_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "utility.hpp"

#include <array>
#include <map>

using bdsg::HashGraph;
using handlegraph::path_handle_t;
using handlegraph::handle_t;

using std::array;
using std::map;


namespace bluntifier {

class PathInfo{
public:
    path_handle_t path_handle;
    uint32_t poa_id;
    bool biclique_side;

    PathInfo(path_handle_t path_handle, uint32_t poa_id, bool biclique_side);
    PathInfo()=default;
};


class Subgraph{
public:
    // The small graph resulting from POA of overlap sequences
    HashGraph graph;

    // In terms of the biclique, record the side of each overlap and the name needed to fetch that path after copying
    array <map <handle_t, PathInfo>, 2> paths_per_handle;

    Subgraph()=default;
};



}

#endif //BLUNTIFIER_SUBGRAPH_HPP
