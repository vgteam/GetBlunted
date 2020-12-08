#include "duplicate_terminus.hpp"
#include "traverse.hpp"
#include "utility.hpp"
#include "handle_to_gfa.hpp"
#include "bdsg/hash_graph.hpp"

#include <iostream>

using bluntifier::duplicate_prefix;
using bluntifier::duplicate_suffix;
using bluntifier::handle_graph_to_gfa;
using bluntifier::source_sink_paths;
using bluntifier::run_command;
using bdsg::HashGraph;

using std::string;
using std::deque;
using std::ofstream;


void test_adjacent_same_lengths_prefix(){
    HashGraph graph;
    string original_sequence = "GATTACA";
    auto h = graph.create_handle(original_sequence);

    deque<handle_t> children;
    deque<size_t> sizes;

    sizes.emplace_back(3);
    sizes.emplace_back(2);
    sizes.emplace_back(2);
    sizes.emplace_back(1);

    {
        string test_path_prefix = "test_adjacent_same_lengths_prefix_" + std::to_string(0);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    duplicate_prefix(graph, sizes, children, h);

    {
        string test_path_prefix = "test_adjacent_same_lengths_prefix_" + std::to_string(1);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    vector<vector<handle_t>> paths = source_sink_paths(graph);

    for (auto& item: children) {
        std::cout << graph.get_id(item) << " " << graph.get_sequence(item) << '\n';
    }

    for (auto& path: paths) {
        string sequence;

        for (auto& item: path) {
            sequence += graph.get_sequence(item);
        }

        if (sequence != original_sequence) {
            throw runtime_error("FAIL: traversal sequence does not equal parent node sequence: "
                                + sequence + " != " + original_sequence);
        }
    }
}


void test_single_value_prefix(){
    HashGraph graph;
    string original_sequence = "GATTACA";
    auto h = graph.create_handle(original_sequence);

    deque<handle_t> children;
    deque<size_t> sizes;

    sizes.emplace_back(1);

    {
        string test_path_prefix = "test_single_value_prefix_" + std::to_string(0);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    duplicate_prefix(graph, sizes, children, h);

    {
        string test_path_prefix = "test_single_value_prefix_" + std::to_string(0);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    vector<vector<handle_t>> paths = source_sink_paths(graph);

    for (auto& item: children) {
        std::cout << graph.get_id(item) << " " << graph.get_sequence(item) << '\n';
    }

    for (auto& path: paths) {
        string sequence;

        for (auto& item: path) {
            sequence += graph.get_sequence(item);
        }

        if (sequence != original_sequence) {
            throw runtime_error("FAIL: traversal sequence does not equal parent node sequence: "
                                + sequence + " != " + original_sequence);
        }
    }
}


void test_zero_value_prefix(){
    HashGraph graph;
    string original_sequence = "GATTACA";
    auto h = graph.create_handle(original_sequence);

    deque<handle_t> children;
    deque<size_t> sizes;

    sizes.emplace_back(1);
    sizes.emplace_back(0);

    {
        string test_path_prefix = "test_zero_value_prefix_" + std::to_string(0);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    duplicate_prefix(graph, sizes, children, h);

    {
        string test_path_prefix = "test_zero_value_prefix_" + std::to_string(1);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    vector<vector<handle_t>> paths = source_sink_paths(graph);

    for (auto& item: children) {
        std::cout << graph.get_id(item) << " " << graph.get_sequence(item) << '\n';
    }

    for (auto& path: paths) {
        string sequence;

        for (auto& item: path) {
            sequence += graph.get_sequence(item);
        }

        if (sequence != original_sequence) {
            throw runtime_error("FAIL: traversal sequence does not equal parent node sequence: "
                                + sequence + " != " + original_sequence);
        }
    }
}


void test_full_length_value_prefix(){
    HashGraph graph;
    string original_sequence = "GATTACA";
    auto h = graph.create_handle(original_sequence);

    deque<handle_t> children;
    deque<size_t> sizes;

    sizes.emplace_back(7);

    {
        string test_path_prefix = "test_full_length_value_prefix_" + std::to_string(0);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    duplicate_prefix(graph, sizes, children, h);

    {
        string test_path_prefix = "test_full_length_value_prefix_" + std::to_string(1);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    vector<vector<handle_t>> paths = source_sink_paths(graph);

    for (auto& item: children) {
        std::cout << graph.get_id(item) << " " << graph.get_sequence(item) << '\n';
    }

    for (auto& path: paths) {
        string sequence;

        for (auto& item: path) {
            sequence += graph.get_sequence(item);
        }

        if (sequence != original_sequence) {
            throw runtime_error("FAIL: traversal sequence does not equal parent node sequence: "
                                + sequence + " != " + original_sequence);
        }
    }
}



void test_adjacent_same_lengths_suffix(){
    HashGraph graph;
    string original_sequence = "GATTACA";
    auto h = graph.create_handle(original_sequence);

    deque<handle_t> children;
    deque<size_t> sizes;

    sizes.emplace_back(3);
    sizes.emplace_back(2);
    sizes.emplace_back(2);
    sizes.emplace_back(1);

    {
        string test_path_prefix = "test_adjacent_same_lengths_suffix_" + std::to_string(0);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    duplicate_suffix(graph, sizes, children, h);

    {
        string test_path_prefix = "test_adjacent_same_lengths_suffix_" + std::to_string(1);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    vector<vector<handle_t>> paths = source_sink_paths(graph);

    for (auto& item: children) {
        std::cout << graph.get_id(item) << " " << graph.get_sequence(item) << '\n';
    }

    for (auto& path: paths) {
        string sequence;

        for (auto& item: path) {
            sequence += graph.get_sequence(item);
        }

        if (sequence != original_sequence) {
            throw runtime_error("FAIL: traversal sequence does not equal parent node sequence: "
                                + sequence + " != " + original_sequence);
        }
    }
}


void test_single_value_suffix(){
    HashGraph graph;
    string original_sequence = "GATTACA";
    auto h = graph.create_handle(original_sequence);

    deque<handle_t> children;
    deque<size_t> sizes;

    sizes.emplace_back(1);

    {
        string test_path_prefix = "test_single_value_suffix_" + std::to_string(0);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    duplicate_suffix(graph, sizes, children, h);

    {
        string test_path_prefix = "test_single_value_suffix_" + std::to_string(1);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    vector<vector<handle_t>> paths = source_sink_paths(graph);

    for (auto& item: children) {
        std::cout << graph.get_id(item) << " " << graph.get_sequence(item) << '\n';
    }

    for (auto& path: paths) {
        string sequence;

        for (auto& item: path) {
            sequence += graph.get_sequence(item);
        }

        if (sequence != original_sequence) {
            throw runtime_error("FAIL: traversal sequence does not equal parent node sequence: "
                                + sequence + " != " + original_sequence);
        }
    }
}


void test_zero_value_suffix(){
    HashGraph graph;
    string original_sequence = "GATTACA";
    auto h = graph.create_handle(original_sequence);

    deque<handle_t> children;
    deque<size_t> sizes;

    sizes.emplace_back(1);
    sizes.emplace_back(0);

    {
        string test_path_prefix = "test_zero_value_suffix_" + std::to_string(0);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    duplicate_suffix(graph, sizes, children, h);

    {
        string test_path_prefix = "test_zero_value_suffix_" + std::to_string(1);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    vector<vector<handle_t>> paths = source_sink_paths(graph);

    for (auto& item: children) {
        std::cout << graph.get_id(item) << " " << graph.get_sequence(item) << '\n';
    }

    for (auto& path: paths) {
        string sequence;

        for (auto& item: path) {
            sequence += graph.get_sequence(item);
        }

        if (sequence != original_sequence) {
            throw runtime_error("FAIL: traversal sequence does not equal parent node sequence: "
                                + sequence + " != " + original_sequence);
        }
    }
}


void test_full_length_value_suffix(){
    HashGraph graph;
    string original_sequence = "GATTACA";
    auto h = graph.create_handle(original_sequence);

    deque<handle_t> children;
    deque<size_t> sizes;

    sizes.emplace_back(7);

    {
        string test_path_prefix = "test_full_length_value_suffix_" + std::to_string(0);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    duplicate_suffix(graph, sizes, children, h);

    {
        string test_path_prefix = "test_full_length_value_suffix_" + std::to_string(1);
        ofstream out(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, out);
        string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + test_path_prefix + ".png";
        run_command(command);
    }

    vector<vector<handle_t>> paths = source_sink_paths(graph);

    for (auto& item: children) {
        std::cout << graph.get_id(item) << " " << graph.get_sequence(item) << '\n';
    }

    for (auto& path: paths) {
        string sequence;

        for (auto& item: path) {
            sequence += graph.get_sequence(item);
        }

        if (sequence != original_sequence) {
            throw runtime_error("FAIL: traversal sequence does not equal parent node sequence: "
                                + sequence + " != " + original_sequence);
        }
    }
}


int main(){
    std::cout << "TESTING PREFIXES:\n";

    std::cout << "TESTING: adjacent_same_lengths:\n";
    test_adjacent_same_lengths_prefix();
    std::cout << "PASS\n\n";

    std::cout << "TESTING: single_value:\n";
    test_single_value_prefix();
    std::cout << "PASS\n\n";

    std::cout << "TESTING: zero_value:\n";
    test_zero_value_prefix();
    std::cout << "PASS\n\n";

    std::cout << "TESTING: full_length_value:\n";
    test_full_length_value_prefix();
    std::cout << "PASS\n\n";

    std::cout << "TESTING SUFFIXES:\n";

    std::cout << "TESTING: adjacent_same_lengths:\n";
    test_adjacent_same_lengths_suffix();
    std::cout << "PASS\n\n";

    std::cout << "TESTING: single_value:\n";
    test_single_value_suffix();
    std::cout << "PASS\n\n";

    std::cout << "TESTING: zero_value:\n";
    test_zero_value_suffix();
    std::cout << "PASS\n\n";

    std::cout << "TESTING: full_length_value:\n";
    test_full_length_value_suffix();
    std::cout << "PASS\n\n";

    return 0;
}



