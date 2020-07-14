#include "handlegraph/handle_graph.hpp"
#include "IncrementalIdMap.hpp"
#include <vector>
#include <string>

using bluntifier::IncrementalIdMap;
using handlegraph::nid_t;
using std::vector;
using std::string;
using std::cerr;


int main(){

    vector <string> names = {"a", "b", "c" , "d"};
    IncrementalIdMap id_map;

    string name;
    nid_t id = -1;

    // Test forward mapping
    for (auto& item: names){
        id = id_map.insert(item);
        cerr << id << '\n';
    }

    // Test reverse mapping
    for (auto& item: names){
        id = id_map.get_id(item);
        name = id_map.get_name(id);

        cerr << id << " " << name << '\n';
    }

    // Test found/missing items
    bool found_id;
    bool found_name;

    found_id = id_map.exists("a");
    cerr << found_id << " " << "a" << '\n';

    found_id = id_map.exists("NO");
    cerr << found_id << " " << "NO" << '\n';

    found_name = id_map.exists(4);
    cerr << found_name << " " << 4 << '\n';

    found_name = id_map.exists(5);
    cerr << found_name << " " << 5 << '\n';

    found_name = id_map.exists(-1);
    cerr << found_name << " " << -1 << '\n';

    // Test duplicate items
    id = id_map.insert("a");

    return 0;
}
