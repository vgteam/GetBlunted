#ifndef BLUNTIFIER_COPY_GRAPH_HPP
#define BLUNTIFIER_COPY_GRAPH_HPP

#ifndef VG_ALGORITHMS_COPY_GRAPH_HPP_INCLUDED
#define VG_ALGORITHMS_COPY_GRAPH_HPP_INCLUDED

/**
 * \file copy_graph.hpp
 *
 * Defines algorithms for copying data between handle graphs
 */

#include <iostream>

#include "handlegraph/handle_graph.hpp"
#include "handlegraph/path_handle_graph.hpp"
#include "handlegraph/mutable_path_mutable_handle_graph.hpp"

using handlegraph::HandleGraph;
using handlegraph::PathHandleGraph;
using handlegraph::MutableHandleGraph;
using handlegraph::MutablePathHandleGraph;
using handlegraph::MutablePathMutableHandleGraph;
using handlegraph::path_handle_t;


namespace bluntifier{

/// Copies the nodes and edges from one graph into another.
void copy_handle_graph(const HandleGraph* from, MutableHandleGraph* into);

/// Copies the nodes, edges, and paths from one graph into another.
void copy_path_handle_graph(const PathHandleGraph* from, MutablePathMutableHandleGraph* into);

/// Copies a path from one graph to another. Nodes and edges to support
/// the path must already exist.
void copy_path(const PathHandleGraph* from, const path_handle_t& path,
               MutablePathHandleGraph* into);

}


#endif


#endif //BLUNTIFIER_COPY_GRAPH_HPP
