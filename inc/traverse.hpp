#ifndef BLUNTIFIER_TRAVERSE_HPP
#define BLUNTIFIER_TRAVERSE_HPP

#include "handlegraph/handle_graph.hpp"
#include <string>
#include <vector>

using handlegraph::HandleGraph;
using handlegraph::handle_t;
using std::string;
using std::vector;
using std::pair;


namespace bluntifier {


vector<vector<handle_t>> source_sink_paths(const HandleGraph& graph) {
    vector<vector<handle_t>> paths;

    graph.for_each_handle([&](const handle_t& handle) {
        bool no_left_nbrs = graph.follow_edges(handle, true, [](const handle_t& dummy) { return false; });

        if (no_left_nbrs) {
            vector<pair<vector<handle_t>, size_t>> stack;
            stack.emplace_back(vector<handle_t>(1, handle), 0);

            while (!stack.empty()) {
                if (stack.back().second == stack.back().first.size()) {
                    stack.pop_back();
                    continue;
                }

                handle_t here = stack.back().first[stack.back().second];
                ++stack.back().second;
                vector<handle_t> next_handles;

                graph.follow_edges(here, false, [&](const handle_t& next) {
                    next_handles.push_back(next);
                });

                if (next_handles.empty()) {
                    paths.emplace_back();
                    for (const auto& record : stack) {
                        paths.back().push_back(record.first[record.second - 1]);
                    }
                } else {
                    stack.emplace_back(move(next_handles), 0);
                }
            }
        }
    });

    return paths;
}


}

#endif //BLUNTIFIER_TRAVERSE_HPP
