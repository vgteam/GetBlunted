#include <iostream>
#include "utility.hpp"

namespace bluntifier {


void run_command(string& argument_string){
    int exit_code = system(argument_string.c_str());

    if (exit_code != 0){
        throw runtime_error("ERROR: command failed to run: " + argument_string);
    }
}


int64_t find_nth_instance(string& s, char c, size_t n){
    if (n < 1){
        throw runtime_error("Error: cannot find nth occurrence for n < 1");
    }

    int64_t n_found = 0;
    for (int64_t i = 0; i<s.size(); i++){
        if (s[i] == c){
            n_found++;
        }
        if (n_found == int64_t(n)){
            return i;
        }
    }

    return -1;
}


int64_t find_nth_instance_from_back(string& s, char c, size_t n){
    if (n < 1){
        throw runtime_error("Error: cannot find nth occurrence for n < 1");
    }

    int64_t n_found = 0;
    for (int64_t i = s.size()-1; i>=0; i--){
        if (s[i] == c){
            n_found++;
        }
        if (n_found == int64_t(n)){
            return i;
        }
    }

    return -1;
}


void trim_to_nth_instance(string& s, char c, size_t n){
    int64_t i = find_nth_instance(s, c, n);

    if (i >= 0){
        s.resize(i+1);
    }
    else{
        throw runtime_error("Error: character " + string(1, c) + " not found in string");
    }
}


void trim_to_nth_instance_from_back(string& s, char c, size_t n){
    int64_t i = find_nth_instance_from_back(s, c, n);

    if (i >= 0){
        s.resize(i+1);
    }
    else{
        throw runtime_error("Error: character " + string(1, c) + " not found in string");
    }
}


string parent_path(string s){
    ///
    /// Pretends the string is a path and returns a copy of the path up to the last non-terminal '/' character
    ///

    // If already root, do nothing
    if (s == "/"){
        return "/";
    }

    if (s.back() == '/'){
        s.resize(s.size()-1);
    }

    trim_to_nth_instance_from_back(s, '/', 1);

    return s;
}


string parent_path(string s, uint64_t n){
    // If already root, do nothing
    if (s == "/"){
        return "/";
    }

    for (uint64_t i=0; i<n; i++){
        s = parent_path(s);
    }

    return s;
}


string join_paths(string a, string b){
    ///
    /// Inefficient method of normalizing and joining paths to ensure that when joined there is one adjoining '/'
    ///

    int64_t i = a.size()-1;
    while (a[i] == '/' and i >= 0){
        a.resize(a.size()-1);
        i--;
    }

    while (b[0] == '/' and not b.empty()){
        b = b.substr(1,b.size()-1);
    }

    return a + "/" + b;
}


}

