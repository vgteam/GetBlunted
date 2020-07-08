#include <experimental/filesystem>
#include <fstream>

#include "gfakluge.hpp"
#include "../../build/gfak/include/gfakluge.hpp"

using std::experimental::filesystem::path;
using std::ifstream;
using std::cout;

using gfak::GFAKluge;


int main(){

    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path relative_gfa_path = "/data/test_gfa1.gfa";
    path absolute_gfa_path = project_directory / relative_gfa_path;

    ifstream file(absolute_gfa_path);

    GFAKluge gfa_reader;

    gfa_reader.parse_gfa_file(file);

    for (auto& element: gfa_reader.get_name_to_seq()){
        cout << element.first << " " << element.second.sequence << '\n';
    }

    return 0;
}

