#ifndef BLUNTIFIER_UTILITY_HPP
#define BLUNTIFIER_UTILITY_HPP

#include <stdexcept>
#include <string>

using std::string;
using std::runtime_error;

namespace bluntifier {

int64_t find_nth_instance(string& s, char c, size_t n);

int64_t find_nth_instance_from_back(string& s, char c, size_t n);

void trim_to_nth_instance(string& s, char c, size_t n);

void trim_to_nth_instance_from_back(string& s, char c, size_t n);

string parent_path(string s);

string parent_path(string s, uint64_t n);

string join_paths(string a, string b);


}

#endif //BLUNTIFIER_UTILITY_HPP
