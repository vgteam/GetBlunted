#ifndef BLUNTIFIER_DUPLICATOR_HPP
#define BLUNTIFIER_DUPLICATOR_HPP

#include "BicliqueCover.hpp"
#include "Biclique.hpp"
#include "NodeInfo.hpp"
#include "OverlapMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "duplicate_terminus.hpp"
#include "utility.hpp"

#include "bdsg/hash_graph.hpp"

using handlegraph::MutablePathDeletableHandleGraph;
using handlegraph::as_integer;
using handlegraph::handle_t;
using bdsg::HashGraph;

namespace bluntifier{


class Duplicator{
public:
    /// Attributes ///
    const vector <vector <BicliqueEdgeIndex> >& node_to_biclique_edge;
    Bicliques& bicliques;
    OverlapMap& overlaps;


    /// Methods ///
    Duplicator(
            const vector <vector <BicliqueEdgeIndex> >& node_to_biclique_edge,
            Bicliques& bicliques,
            OverlapMap& overlaps);

    void duplicate_termini(MutablePathDeletableHandleGraph& gfa_graph);

private:
    void repair_edges(
            MutablePathDeletableHandleGraph& gfa_graph,
            const array <map <size_t, handle_t>, 2>& biclique_side_to_child,
            nid_t old_node_id,
            handle_t old_handle,
            handle_t old_handle_flipped);

    void remove_participating_edges(
            const array <deque <size_t>, 2>& sorted_bicliques_per_side,
            MutablePathDeletableHandleGraph& gfa_graph,
            nid_t parent_node);
};


}

#endif //BLUNTIFIER_DUPLICATOR_HPP
