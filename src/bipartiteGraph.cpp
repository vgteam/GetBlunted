/**
 * \file bipartiteGraph.cpp
 *
 * Implements algorithm for computing the biclique cover of a bipartite graph.
 */
#include "bipartiteGraph.hpp"

namespace bluntifier {

BipartiteGraph::BipartiteGraph(const HandleGraph& graph,
                               const bipartition& partition)
: graph(&graph),
{
    _partition.first.reserve(partition.first.size());
    _partition.second.reserve(partition.second.size());
    // sort to remove system dependent behavior
    sort(_partition.first.begin(), _partition.first.end());
    sort(_partition.second.begin(), _partition.second.end());
    // map the handles back to their index as well
    left_partition_index.reserve(_partition.first.size());
    right_partition_index.reserve(_partition.second.size());
    for (size_t i = 0; i < _partition.first.size(); ++i) {
        left_partition_index[_partition.first[i]] = i;
    }
    for (size_t i = 0; i < _partition.second.size(); ++i) {
        right_partition_index[_partition.second[i]] = i;
    }
}

BipartiteGraph::BipartiteGraph(const HandleGraph& graph,
                               const ordered_bipartition& partition)
: graph(&graph),
  _partition(partition)
{
    // map the handles back to their index as well
    left_partition_index.reserve(_partition.first.size());
    right_partition_index.reserve(_partition.second.size());
    for (size_t i = 0; i < _partition.first.size(); ++i) {
        left_partition_index[_partition.first[i]] = i;
    }
    for (size_t i = 0; i < _partition.second.size(); ++i) {
        right_partition_index[_partition.second[i]] = i;
    }
}

BipartiteGraph::~BipartiteGraph() {
    
}

BipartiteGraph::const_iterator BipartiteGraph::left_begin() const {
    return _partition.first.begin();
}

BipartiteGraph::const_iterator BipartiteGraph::left_end() const {
    return _partition.first.end();
}

BipartiteGraph::const_iterator BipartiteGraph::left_iterator(const handle_t node) const {
    return _partition.first.begin() + left_partition_index.at(node);
}

BipartiteGraph::size_t BipartiteGraph::left_size() const {
    return _partition.first.size();
}

BipartiteGraph::const_iterator BipartiteGraph::right_begin() const {
    return _partition.second.begin();
}

BipartiteGraph::const_iterator BipartiteGraph::right_end() const {
    return _partition.second.end();
}

BipartiteGraph::const_iterator BipartiteGraph::right_iterator(const handle_t node) const {
    return _partition.second.begin() + right_partition_index.at(node);
}

BipartiteGraph::size_t BipartiteGraph::right_size() const {
    return _partition.second.size();
}

bool BipartiteGraph::for_each_adjacent_side(const handle_t& side,
                                            const function<bool(handle_t)>& lambda) const {
    return graph->follow_edges(side, false, [&](const handle_t& neighbor) {
        return lambda(graph->flip(neighbor));
    });
}

const ordered_bipartition& BipartiteGraph::bipartition() const {
    return partition;
}
}