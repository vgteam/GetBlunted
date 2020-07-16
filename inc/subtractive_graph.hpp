#ifndef BLUNTIFIER_SUBTRACTIVE_GRAPH_HPP_INCLUDED
#define BLUNTIFIER_SUBTRACTIVE_GRAPH_HPP_INCLUDED

/** \file
 * subtractive_graph.hpp: defines a handle graph implementation of a subgraph
 */

#include <unordered_set>
#include <functional>
#include <limits>
#include "handlegraph/handle_graph.hpp"
#include "handlegraph/expanding_overlay_graph.hpp"

namespace bluntifier {

using handlegraph::nid_t;
using handlegraph::edge_t;
using handlegraph::ExpandingOverlayGraph;
using handlegraph::handle_t;
using std::unordered_set;
using std::function;

using namespace std;

    /**
     * A HandleGraph implementation that acts as a subgraph of some other HandleGraph
     * using a layer of indirection. Edges may be removed relative to the underlying graph
     * but no other subsetting operations are allowed.
     */
    class SubtractiveHandleGraph : public ExpandingOverlayGraph {
    public:
        
        /// Initialize with a super graph and nodes returned by iterators to handles
        /// from the super graph
        SubtractiveHandleGraph(const HandleGraph& super);
        
        /// Remove an edge that is present in the super graph. Removing the same
        /// edge multiple times has no additional effect.
        void subtract_edge(const handle_t& left, const handle_t& right);
        
        //////////////////////////
        /// HandleGraph interface
        //////////////////////////
        
        // Method to check if a node exists by ID
        bool has_node(nid_t node_id) const;
        
        /// Look up the handle for the node with the given ID in the given orientation
        handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;
        
        /// Get the ID from a handle
        nid_t get_id(const handle_t& handle) const;
        
        /// Get the orientation of a handle
        bool get_is_reverse(const handle_t& handle) const;
        
        /// Invert the orientation of a handle (potentially without getting its ID)
        handle_t flip(const handle_t& handle) const;
        
        /// Get the length of a node
        size_t get_length(const handle_t& handle) const;
        
        /// Get the sequence of a node, presented in the handle's local forward
        /// orientation.
        string get_sequence(const handle_t& handle) const;
        
        /// Loop over all the handles to next/previous (right/left) nodes. Passes
        /// them to a callback which returns false to stop iterating and true to
        /// continue. Returns true if we finished and false if we stopped early.
        bool follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
        
        /// Loop over all the nodes in the graph in their local forward
        /// orientations, in their internal stored order. Stop if the iteratee
        /// returns false. Can be told to run in parallel, in which case stopping
        /// after a false return value is on a best-effort basis and iteration
        /// order is not defined.
        bool for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
        
        /// Return the number of nodes in the graph
        size_t get_node_count() const;
        
        /// Return the smallest ID in the graph, or some smaller number if the
        /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
        nid_t min_node_id() const;
        
        /// Return the largest ID in the graph, or some larger number if the
        /// largest ID is unavailable. Return value is unspecified if the graph is empty.
        nid_t max_node_id() const;
        
        ///////////////////////////////////
        /// ExpandingOverlayGraph interface
        ///////////////////////////////////
        
        /**
         * Returns the handle in the underlying graph that corresponds to a handle in the
         * overlay
         */
        handle_t get_underlying_handle(const handle_t& handle) const;
        
    private:
        const HandleGraph* super = nullptr;
        unordered_set<edge_t> subtracted_edges;
        
    };
}

#endif /* BLUNTIFIER_SUBTRACTIVE_GRAPH_HPP_INCLUDED */
