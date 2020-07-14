#ifndef BLUNTIFIER_INCREMENTALIDMAP_HPP
#define BLUNTIFIER_INCREMENTALIDMAP_HPP

#include "handlegraph/handle_graph.hpp"
#include <string>
#include <vector>
#include <map>
#include <memory>

namespace bluntifier{

using std::string;
using std::vector;
using std::unordered_map;
using handlegraph::nid_t;
using std::unique_ptr;


class IncrementalIdMap {
public:
    /// Attributes ///

    vector <unique_ptr <string> > names;
    unordered_map <string, int64_t> ids;

    /// Methods ///

    IncrementalIdMap();

    // Add a node ID to the running list, do whatever needs to be done to make sure the mapping is reversible, and then
    // return its incremental ID, based on the number of nodes added so far
    int64_t insert(const string& s);

    // Find the original node ID from its integer ID
    string get_name(int64_t id);
    int64_t get_id(const string& name);

    // Check if key/value has been added already, returns true if it exists
    bool exists(const string& name);
    bool exists(int64_t id);
};


}
#endif //BLUNTIFIER_INCREMENTALID_HPP
