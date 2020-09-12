#include "IncrementalIdMap.hpp"
#include <memory>

using std::runtime_error;
using std::unique_ptr;
using std::make_unique;

namespace bluntifier{

IncrementalIdMap::IncrementalIdMap()=default;


int64_t IncrementalIdMap::insert(const string& s) {
    if (exists(s)){
        throw runtime_error("Error: attempted to insert duplicate key: " + s);
    }

    // Make a copy of the node name string, and allocate a pointer to it
    names.push_back(make_unique<string>(s));

    // Create an integer node ID (starting from 1)
    int64_t id = names.size();

    // Create a reverse mapping by dereferencing the pointer
    ids.insert({*names.back(),id});

    // For convenience, return the ID number that was generated
    return id;
}


bool IncrementalIdMap::exists(const string& name){
    return (ids.find(name) != ids.end());
}


bool IncrementalIdMap::exists(int64_t id){
    return (id >= 0 and id <= names.size());
}


int64_t IncrementalIdMap::get_id(const string& name) const{
    return ids.at(name);
}


string IncrementalIdMap::get_name(int64_t id) const{
    return *names[id-1];
}


}