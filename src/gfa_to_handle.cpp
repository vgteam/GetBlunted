#include "gfa_to_handle.hpp"

using handlegraph::handle_t;
using handlegraph::path_handle_t;
using bdsg::MutablePathHandleGraph;
using bdsg::HandleGraph;

namespace bluntifier {


nid_t parse_gfa_sequence_id(const string& s, IncrementalIdMap<string>& id_map) {
    nid_t id;

    if (id_map.exists(s)){
        id = id_map.get_id(s);
    }
    else{
        id = id_map.insert(s);
    }

    return id;
}


void validate_gfa_edge(const gfak::edge_elem& e) {
    if (e.source_name.empty()) {
        throw GFAFormatError("Error:[gfa_to_handle_graph] Found edge record with missing source name");
    }
    if (e.sink_name.empty()) {
        throw GFAFormatError("Error:[gfa_to_handle_graph] Found edge record with missing sink name");
    }
}


string process_raw_gfa_path_name(const string& path_name_raw) {
    string processed = path_name_raw;
    processed.erase(remove_if(processed.begin(), processed.end(),
                              [](char c) { return isspace(c); }),
                    processed.end());
    return processed;
}


void gfa_to_handle_graph_in_memory(
        istream& in,
        MutableHandleGraph& graph,
        gfak::GFAKluge& gg,
        IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps) {

    if (!in) {
        throw std::ios_base::failure("Error:[gfa_to_handle_graph] Couldn't open input stream");
    }
    gg.parse_gfa_file(in);

    // create nodes
    for (const auto& seq_record : gg.get_name_to_seq()) {
        graph.create_handle(seq_record.second.sequence, parse_gfa_sequence_id(seq_record.first, id_map));
    }

    // create edges
    size_t num_self_overlap_warnings = 0;
    size_t max_num_self_overlap_warnings = 10;
    for (const auto& links_record : gg.get_seq_to_edges()) {
        for (const auto& edge : links_record.second) {
            validate_gfa_edge(edge);

            const nid_t source_id = parse_gfa_sequence_id(edge.source_name, id_map);
            const nid_t sink_id = parse_gfa_sequence_id(edge.sink_name, id_map);

            if (not graph.has_node(source_id)){
                throw runtime_error("ERROR: gfa link (" + edge.source_name + "->" + edge.sink_name + ") "
                                         "contains non-existent node: " + edge.source_name);
            }

            if (not graph.has_node(sink_id)){
                throw runtime_error("ERROR: gfa link (" + edge.source_name + "->" + edge.sink_name + ") "
                                         "contains non-existent node: " + edge.sink_name);
            }


            // note: we're counting on implementations de-duplicating edges
            handle_t a = graph.get_handle(source_id, not edge.source_orientation_forward);
            handle_t b = graph.get_handle(sink_id, not edge.sink_orientation_forward);
            
            // Update the overlap map
            Alignment alignment(edge.alignment);

            pair<size_t, size_t> lengths;
            alignment.compute_lengths(lengths);

            bool valid = true;
            if (lengths.first > graph.get_length(a)){
                cerr << "WARNING: skipping overlap for which sum of cigar operations is > SOURCE node length: "
                     << edge.source_name << "->" << edge.sink_name << '\n' << '\n';
                valid = false;
            }
            if (lengths.second > graph.get_length(b)){
                cerr << "WARNING: skipping overlap for which sum of cigar operations is > SINK node length: "
                     << edge.source_name << "->" << edge.sink_name << '\n' << '\n';
                valid = false;
            }
            if (b == graph.flip(a)) {
                if (num_self_overlap_warnings++ < max_num_self_overlap_warnings) {
                    cerr << "WARNING: skipping overlap that reverses back onto the same sequence (in de Bruijn graphs, these can be avoided by choosing an odd-numbered overlap length): " << e.source_name << '\n';
                }
                if (num_self_overlap_warnings == max_num_self_overlap_warnings) {
                    cerr << "suppressing further warnings\n";
                }
                valid = false;
            }

            if (valid) {
                graph.create_edge(a, b);
                overlaps.insert(alignment, a, b);
            }
        }
    }
}


void gfa_to_handle_graph_on_disk(
        const string& filename,
        MutableHandleGraph& graph,
        bool try_id_increment_hint,
        gfak::GFAKluge& gg,
        IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps) {

    // adapted from
    // https://github.com/vgteam/odgi/blob/master/src/gfa_to_handle.cpp

    if (try_id_increment_hint) {

        // find the minimum ID
        nid_t min_id = numeric_limits<nid_t>::max();
        gg.for_each_sequence_line_in_file(filename.c_str(), [&](gfak::sequence_elem s) {
            min_id = std::min(min_id, parse_gfa_sequence_id(s.name, id_map));
        });

        if (min_id != numeric_limits<nid_t>::max()) {
            // we found the min, set it as the increment
            graph.set_id_increment(min_id);
        }
    }

    // add in all nodes
    gg.for_each_sequence_line_in_file(filename.c_str(), [&](gfak::sequence_elem s) {

        graph.create_handle(s.sequence, parse_gfa_sequence_id(s.name, id_map));

    });

    // add in all edges
    size_t num_self_overlap_warnings = 0;
    size_t max_num_self_overlap_warnings = 10;
    gg.for_each_edge_line_in_file(const_cast<char*>(filename.c_str()), [&](gfak::edge_elem e) {
        validate_gfa_edge(e);

        const nid_t source_id = parse_gfa_sequence_id(e.source_name, id_map);
        const nid_t sink_id = parse_gfa_sequence_id(e.sink_name, id_map);

        if (not graph.has_node(source_id)){
            throw runtime_error("ERROR: gfa link (" + e.source_name + "->" + e.sink_name + ") "
                                    "contains non-existent node: " + e.source_name);
        }

        if (not graph.has_node(sink_id)){
            throw runtime_error("ERROR: gfa link (" + e.source_name + "->" + e.sink_name + ") "
                                     "contains non-existent node: " + e.sink_name);
        }

        handle_t a = graph.get_handle(source_id, not e.source_orientation_forward);
        handle_t b = graph.get_handle(sink_id, not e.sink_orientation_forward);

        // Update the overlap map
        Alignment alignment(e.alignment);

        pair<size_t, size_t> lengths;
        alignment.compute_lengths(lengths);

        bool valid = true;
        if (lengths.first > graph.get_length(a)){
            cerr << "WARNING: skipping overlap for which sum of cigar operations is > SOURCE node length: "
                 << e.source_name << "->" << e.sink_name << '\n' << '\n';
            valid = false;
        }
        if (lengths.second > graph.get_length(b)){
            cerr << "WARNING: skipping overlap for which sum of cigar operations is > SINK node length: "
                 << e.source_name << "->" << e.sink_name << '\n' << '\n';
            valid = false;
        }
        if (b == graph.flip(a)) {
            if (num_self_overlap_warnings++ < max_num_self_overlap_warnings) {
                cerr << "WARNING: skipping overlap that reverses back onto the same sequence (in de Bruijn graphs, these can be avoided by choosing an odd-numbered overlap length): " << e.source_name << '\n';
            }
            if (num_self_overlap_warnings == max_num_self_overlap_warnings) {
                cerr << "suppressing further warnings\n";
            }
            valid = false;
        }

        if (valid) {
            graph.create_edge(a, b);
            overlaps.insert(alignment, a, b);
        }
    });
}


/// Parse nodes and edges and load them into the given GFAKluge.
/// If the input is a seekable file, filename will be filled in and unseekable will be nullptr.
/// If the input is not a seekable file, filename may be filled in, and unseekable will be set to a stream to read from.
void gfa_to_handle_graph_load_graph(
        const string& filename,
        istream* unseekable,
        MutableHandleGraph& graph,
        bool try_id_increment_hint,
        gfak::GFAKluge& gg,
        IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps) {

    if (graph.get_node_count() > 0) {
        throw invalid_argument("Error:[gfa_to_handle_graph] Must parse GFA into an empty graph");
    }

    if (!unseekable) {
        // Do the from-disk path
        gfa_to_handle_graph_on_disk(filename, graph, try_id_increment_hint, gg, id_map, overlaps);
    } else {
        // Do the path for streams

        if (try_id_increment_hint) {
            // The ID increment hint can't be done.
            cerr
                    << "Warning:[gfa_to_handle_graph] Skipping node ID increment hint because input stream for GFA does not support seeking. "
                    << "If performance suffers, consider using an alternate graph implementation or reading GFA from hard disk."
                    << endl;
        }

        gfa_to_handle_graph_in_memory(*unseekable, graph, gg, id_map, overlaps);
    }
}


/// After the given GFAKluge has been populated with nodes and edges, load path information.
/// If the input is a seekable file, filename will be filled in and unseekable will be nullptr.
/// If the input is not a seekable file, filename may be filled in, and unseekable will be set to a stream to read from.
void gfa_to_handle_graph_add_paths(const string& filename, istream* unseekable, MutablePathHandleGraph& graph,
                                   gfak::GFAKluge& gg, IncrementalIdMap<string>& id_map) {


    if (!unseekable) {
        // Input is from a seekable file on disk.

        // add in all paths
        gg.for_each_path_element_in_file(filename.c_str(), [&](const string& path_name_raw,
                                                               const string& node_id,
                                                               bool is_rev,
                                                               const string& cigar,
                                                               bool is_empty,
                                                               bool is_circular) {
            // remove white space in path name
            // TODO: why?
            string path_name = process_raw_gfa_path_name(path_name_raw);

            // get either the existing path handle or make a new one
            path_handle_t path;
            if (not graph.has_path(path_name)) {
                path = graph.create_path_handle(path_name, is_circular);
            } else {
                path = graph.get_path_handle(path_name);
            }

            // add the step
            handle_t step = graph.get_handle(parse_gfa_sequence_id(node_id, id_map), is_rev);
            graph.append_step(path, step);
        });
    } else {

        // gg will have parsed the GFA file in the non-path part of the algorithm
        // No reading to do.

        // create paths
        for (const auto& path_record : gg.get_name_to_path()) {

            // process this to match the disk backed implementation
            // TODO: why?
            string path_name = process_raw_gfa_path_name(path_record.first);
            path_handle_t path = graph.create_path_handle(path_name);

            for (size_t i = 0; i < path_record.second.segment_names.size(); ++i) {
                handle_t step = graph.get_handle(parse_gfa_sequence_id(path_record.second.segment_names.at(i), id_map),
                                                  not path_record.second.orientations.at(i));
                graph.append_step(path, step);
            }
        }
    }


}


void gfa_to_handle_graph(
        const string& filename,
        MutableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps,
        bool try_from_disk,
        bool try_id_increment_hint) {

    // What stream should we read from (isntead of opening the file), if any?
    istream* unseekable = nullptr;

    // If we open a file, it will live here.
    unique_ptr<ifstream> opened;

    if (filename == "-") {
        // Read from standard input
        unseekable = &cin;
    } else if (not try_from_disk) {
        // The file may be seekable actually, but we don't want to use the
        // seekable-file codepath for some reason.
        opened = make_unique<ifstream>(filename);
        if (not opened) {
            throw std::ios_base::failure("Error:[gfa_to_handle_graph] Couldn't open file " + filename);
        }
        unseekable = opened.get();
    }

    gfak::GFAKluge gg;
    gfa_to_handle_graph_load_graph(filename, unseekable, graph, try_id_increment_hint, gg, id_map, overlaps);
}


void gfa_to_path_handle_graph(
        const string& filename,
        MutablePathMutableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps,
        bool try_from_disk,
        bool try_id_increment_hint) {


    // What stream should we read from (instead of opening the file), if any?
    istream* unseekable = nullptr;

    // If we open a file, it will live here.
    unique_ptr<ifstream> opened;

    if (filename == "-") {
        // Read from standard input
        unseekable = &cin;
    } else if (not try_from_disk) {
        // The file may be seekable actually, but we don't want to use the
        // seekable-file codepath for some reason.
        opened = make_unique<ifstream>(filename);
        if (not opened) {
            throw std::ios_base::failure("Error:[gfa_to_handle_graph] Couldn't open file " + filename);
        }
        unseekable = opened.get();
    }

    gfak::GFAKluge gg;
    gfa_to_handle_graph_load_graph(filename, unseekable, graph, try_id_increment_hint, gg, id_map, overlaps);

    // TODO: Deduplicate everything other than this line somehow.
    gfa_to_handle_graph_add_paths(filename, unseekable, graph, gg, id_map);
}


void gfa_to_path_handle_graph_in_memory(
        istream& in,
        MutablePathMutableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        OverlapMap& overlaps) {
    gfak::GFAKluge gg;
    gfa_to_handle_graph_load_graph("", &in, graph, false, gg, id_map, overlaps);
    gfa_to_handle_graph_add_paths("", &in, graph, gg, id_map);

}

}
