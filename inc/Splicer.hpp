#ifndef BLUNTIFIER_SPLICER_HPP
#define BLUNTIFIER_SPLICER_HPP

#include "handlegraph/handle_graph.hpp"
#include "bdsg/hash_graph.hpp"

#include "Pileup.hpp"

using handlegraph::MutablePathMutableHandleGraph;

using handlegraph::handle_t;
using bluntifier::SpliceData;

using std::numeric_limits;
using std::cout;
using std::map;

namespace bluntifier {

class Splicer {
public:
    /// Attributes ///

    MutablePathMutableHandleGraph& gfa_graph;
    const vector<vector<SpliceData> >& splice_sites;
    const size_t node_id;


    /// Methods ///

    Splicer(MutablePathMutableHandleGraph& gfa_graph, vector<vector<SpliceData> >& splice_sites, size_t node_id);

    void find_duplication_sites(vector<pair<size_t, size_t> >& left_sites, vector<pair<size_t, size_t> >& right_sites);

    void find_overlapping_overlaps();

    void duplicate_terminus(size_t site_index);

    void duplicate_all_termini();
};

}

#endif //BLUNTIFIER_SPLICER_HPP
