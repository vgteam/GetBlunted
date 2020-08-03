#include "handlegraph/handle_graph.hpp"
#include "IncrementalIdMap.hpp"
#include <vector>
#include <string>

using bluntifier::IncrementalIdMap;
using handlegraph::nid_t;
using std::runtime_error;
using std::exception;
using std::vector;
using std::string;
using std::cerr;


int main(){

    vector <string> names = {"a", "b", "c" , "d"};
    vector <nid_t> ids = {1, 2, 3, 4};

    IncrementalIdMap id_map;

    string name;
    nid_t id = -1;

    // Test forward mapping
    for (size_t i=0; i<names.size(); i++){
        id = id_map.insert(names[i]);
        if (id != ids[i]){
            throw runtime_error("FAIL: incorrect mapping");
        }

        cerr << id << '\n';
    }

    // Test reverse mapping
    for (size_t i=0; i<names.size(); i++){
        id = id_map.get_id(names[i]);
        name = id_map.get_name(id);
        if (name != names[i]){
            throw runtime_error("FAIL: incorrect mapping");
        }
        cerr << id << " " << name << '\n';
    }

    // Test found/missing items
    bool found_id;
    bool found_name;

    found_id = id_map.exists("a");
    cerr << found_id << " " << "a" << '\n';

    if (not found_id){
        throw runtime_error("FAIL: item retrieval");
    }

    found_id = id_map.exists("NO");
    cerr << found_id << " " << "NO" << '\n';

    if (found_id){
        throw runtime_error("FAIL: item retrieval");
    }

    found_name = id_map.exists(4);
    cerr << found_name << " " << 4 << '\n';

    if (not found_name){
        throw runtime_error("FAIL: name retrieval");
    }

    found_name = id_map.exists(5);
    cerr << found_name << " " << 5 << '\n';

    if (found_name){
        throw runtime_error("FAIL: name retrieval");
    }

    found_name = id_map.exists(-1);
    cerr << found_name << " " << -1 << '\n';

    if (found_name){
        throw runtime_error("FAIL: name retrieval");
    }

    bool pass = false;
    // Test duplicate items
    try {
        id = id_map.insert("a");
    }
    catch (exception& e){
        cerr << e.what() << '\n';
        pass = true;
    }

    if (not pass){
        throw runtime_error("FAIL: map allows duplicate entry");
    }

    cerr << "PASS\n";

    return 0;
}
