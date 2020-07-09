#include <experimental/filesystem>
#include <fstream>

#include "gfakluge.hpp"
#include "utility.hpp"

using std::ifstream;
using std::cout;

using bluntifier::parent_path;
using bluntifier::join_paths;
using gfak::GFAKluge;


int main(){

    string script_path = __FILE__;
    string project_directory = parent_path(script_path, 3);

    // Get test VCF path
    string relative_gfa_path = "/data/test_gfa1.gfa";
    string absolute_gfa_path = join_paths(project_directory, relative_gfa_path);

    ifstream file(absolute_gfa_path);

    GFAKluge gfa_reader;

    gfa_reader.parse_gfa_file(file);

    for (auto& element: gfa_reader.get_name_to_seq()){
        cout << element.first << " " << element.second.sequence << '\n';
    }

    return 0;
}

