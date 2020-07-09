#include <iostream>
#include "utility.hpp"

using std::cerr;
using std::stoi;
using bluntifier::join_paths;
using bluntifier::parent_path;
using bluntifier::find_nth_instance;
using bluntifier::trim_to_nth_instance;
using bluntifier::find_nth_instance_from_back;
using bluntifier::trim_to_nth_instance_from_back;


bool test_find_nth_instance(){
    string s = "012";

    cerr << "Testing front, middle, back\n";
    for (auto& item: s){
        auto true_index = stoi(string(1,item));
        if (true_index != find_nth_instance(s, item, 1)){
            return false;
        }
    }

    cerr << "Testing nonexistent\n";
    if (find_nth_instance(s, 'a', 1) != -1){
        return false;
    }

    cerr << "Testing fewer than requested\n";
    if (find_nth_instance(s, '1', 2) != -1){
        return false;
    }

    return true;
}


bool test_find_nth_instance_from_back(){
    string s = "012";

    cerr << "Testing front, middle, back\n";
    for (auto& item: s){
        auto true_index = stoi(string(1,item));
        auto found_index = find_nth_instance_from_back(s, item, 1);
        if (true_index != found_index){
            return false;
        }
    }

    cerr << "Testing nonexistent\n";
    if (find_nth_instance_from_back(s, 'a', 1) != -1){
        return false;
    }

    cerr << "Testing fewer than requested\n";
    if (find_nth_instance_from_back(s, '1', 2) != -1){
        return false;
    }

    return true;
}


bool test_trim_to_nth_instance(){
    string s = "012";

    trim_to_nth_instance(s, '1', 1);

    if (s != "01"){
        return false;
    }

    return true;
}


bool test_trim_to_nth_instance_from_back(){
    string s = "012";

    trim_to_nth_instance_from_back(s, '1', 1);

    if (s != "01"){
        return false;
    }

    return true;
}


bool test_parent_path(){
    string path = "/root/dir1/filename.ext";

    path = parent_path(path);
    if (path != "/root/dir1/"){
        return false;
    }
    path = parent_path(path);
    if (path != "/root/"){
        return false;
    }
    path = parent_path(path);
    if (path != "/"){
        return false;
    }
    path = parent_path(path);
    if (path != "/"){
        return false;
    }

    return true;
}


bool test_join_paths(){
    string a;
    string b;
    string c;

    a = "/a/b//";
    b = "/c/d";

    c = join_paths(a, b);
    if (c != "/a/b/c/d"){
        return false;
    }
    c = join_paths(b, a);
    if (c != "/c/d/a/b//"){
        return false;
    }
    c = join_paths(a, a);
    if (c != "/a/b/a/b//"){
        return false;
    }
    c = join_paths(b, b);
    if (c != "/c/d/c/d"){
        return false;
    }

    return true;
}


int main(){
    bool pass;

    pass = test_find_nth_instance();
    if (not pass){
        throw runtime_error("FAIL: find_nth_instance");
    }
    pass = test_find_nth_instance_from_back();
    if (not pass){
        throw runtime_error("FAIL: find_nth_instance_from_back");
    }
    pass = test_trim_to_nth_instance();
    if (not pass){
        throw runtime_error("FAIL: trim_to_nth_instance");
    }
    pass = test_trim_to_nth_instance_from_back();
    if (not pass){
        throw runtime_error("FAIL: trim_to_nth_instance_from_back");
    }
    pass = test_parent_path();
    if (not pass){
        throw runtime_error("FAIL: parent_path");
    }
    pass = test_join_paths();
    if (not pass){
        throw runtime_error("FAIL: parent_path");
    }

    cerr << "PASS\n";

    return 0;
}

