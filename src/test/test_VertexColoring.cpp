#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>

#include "bdsg/hash_graph.hpp"
#include "BicliqueCover.hpp"

using bluntifier::VertexColoring;

using std::vector;
using std::set;
using std::sort;
using std::cerr;
using std::endl;
using std::function;
using std::pair;
using std::max_element;


class TestVertexColoring : public VertexColoring {
public:
    TestVertexColoring(const vector<vector<size_t>>& graph) : VertexColoring(graph) {}
    TestVertexColoring() = default;
    using VertexColoring::maximal_independent_sets;
    using VertexColoring::lawlers_algorithm;
    using VertexColoring::greedy_coloring;
    using VertexColoring::degree_ordering;
    using VertexColoring::least_first_ordering;
    using VertexColoring::random_ordering;
    using VertexColoring::interchange_greedy_coloring;
    using VertexColoring::lower_bound;
};

bool is_valid_coloring(const vector<vector<size_t>>& graph, const vector<size_t>& coloring) {
    for (size_t i = 0; i < graph.size(); ++i) {
        for (auto j : graph[i]) {
            if (coloring[i] == coloring[j]) {
                return false;
            }
        }
    }
    return true;
}

int main(){
    
    
    // test the maximal independent sets algorithm
    {
        vector<vector<size_t>> graph{
            {1, 3},
            {0, 2},
            {1, 3},
            {0, 2}
        };
        vector<size_t> index_translator(graph.size());
        for (size_t i = 0; i < index_translator.size(); ++i) {
            index_translator[i] = i;
        }
        
        TestVertexColoring tester;
        auto ind_sets = tester.maximal_independent_sets(graph, index_translator);
        
        // 1010
        uint16_t set1 = 10;
        // 0101
        uint16_t set2 = 5;
        
        if (ind_sets.size() != 2) {
            return 1;
        }
        if (ind_sets[0] != set1 && ind_sets[1] != set1) {
            return 1;
        }
        if (ind_sets[0] != set2 && ind_sets[1] != set2) {
            return 1;
        }
    }
    // a harder test case
    {
        vector<vector<size_t>> graph{
            {1, 2},
            {0, 2, 3},
            {0, 1, 4},
            {1, 4, 5},
            {2, 3, 5},
            {3, 4}
        };
        vector<size_t> index_translator(graph.size());
        for (size_t i = 0; i < index_translator.size(); ++i) {
            index_translator[i] = i;
        }
        
        TestVertexColoring tester;
        auto ind_sets = tester.maximal_independent_sets(graph, index_translator);
        set<uint16_t> got(ind_sets.begin(), ind_sets.end());
        set<uint16_t> truth{
            ((1 << 0)|(1 << 3)),
            ((1 << 0)|(1 << 4)),
            ((1 << 0)|(1 << 5)),
            ((1 << 1)|(1 << 4)),
            ((1 << 1)|(1 << 5)),
            ((1 << 2)|(1 << 3)),
            ((1 << 2)|(1 << 5))
        };
        
        if (got != truth) {
            return 1;
        }
    }
    
    
    // test the exact vertex cover algorithm
    {
        // prism of two tetrahedron, should have a 4-coloring
        vector<vector<size_t>> graph{
            {1, 2, 3, 4},
            {0, 2, 3, 5},
            {0, 1, 3, 6},
            {0, 1, 2, 7},
            {0, 5, 6, 7},
            {1, 4, 6, 7},
            {2, 4, 5, 7},
            {3, 4, 5, 6}
        };

        TestVertexColoring tester(graph);
        auto coloring = tester.lawlers_algorithm();

        if (coloring.size() != graph.size()) {
            return 1;
        }
        // should be able to find a 4-coloring (colors start at 0)
        if (*max_element(coloring.begin(), coloring.end()) != 3) {
            return 1;
        }
        // is a vertex coloring
        for (size_t i = 0; i < graph.size(); ++i) {
            for (size_t j : graph[i]) {
                if (coloring[i] == coloring[j]) {
                    return 1;
                }
            }
        }
    }
    
    {
        vector<vector<size_t>> graph{
            {1, 2, 3},
            {0, 2},
            {0, 1},
            {0}
        };
        
        TestVertexColoring tester(graph);
        auto ordering = tester.degree_ordering();
        
        if (ordering.size() != 4) {
            return 1;
        }
        if (ordering[0] != 0) {
            return 1;
        }
        if (!((ordering[1] == 1 && ordering[2] == 2) ||
              (ordering[1] == 2 && ordering[2] == 1))) {
            return 1;
        }
        if (ordering[3] != 3) {
            return 1;
        }
    }
    
    // last first should work its way inward into a path
    {
        vector<vector<size_t>> graph{
            {1},
            {0, 2},
            {1, 3},
            {2, 4},
            {3, 5},
            {4, 6},
            {5, 7},
            {6, 8},
            {7, 9},
            {8, 10},
            {9}
        };
        
        TestVertexColoring tester(graph);
        auto ordering = tester.least_first_ordering();
        
        size_t first = 0;
        size_t last = graph.size() - 1;
        for (auto i : ordering) {
            if (i == first) {
                ++first;
            }
            else if (i == last) {
                --last;
            }
            else {
                return 1;
            }
        }
        if (first < last) {
            return 1;
        }
    }
    
    {
        vector<vector<size_t>> graph{
            {1},
            {0, 2},
            {1, 3},
            {2, 4},
            {3, 5},
            {4, 6},
            {5, 7},
            {6, 8},
            {7, 9},
            {8, 10},
            {9}
        };
        
        TestVertexColoring tester(graph);
        uint64_t seed = 0;
        for (size_t i = 0; i < 20; ++i) {
            vector<size_t> order = tester.random_ordering(seed);
            vector<bool> is_included(graph.size(), false);
            for (auto i : order) {
                is_included[i] = true;
            }
            for (auto incl : is_included) {
                if (!incl) {
                    return 1;
                }
            }
        }
        
        {
            vector<vector<size_t>> graph{
                {1, 4},
                {0},
                {3},
                {2, 4},
                {0, 3, 7},
                {6},
                {5, 7},
                {4, 6}
            };
            
            vector<size_t> order{0, 1, 2, 3, 4, 5, 6, 7};
            
            TestVertexColoring tester(graph);
            auto coloring = tester.interchange_greedy_coloring(order);
            if (*max_element(coloring.begin(), coloring.end()) != 1) {
                return 1;
            }
            if (!is_valid_coloring(graph, coloring)) {
                return 1;
            }
        }
    }
    
    cerr << "vertex coloring tests successful" << endl;
    return 0;
}

