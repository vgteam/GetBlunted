#include "Bluntifier.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"

#include <iostream>
#include <getopt.h>

using bluntifier::Bluntifier;
using ghc::filesystem::path;
using std::ifstream;
using std::cerr;
using std::cout;
using std::endl;

int main(int argc, char **argv){

    path gfa_path;
    path provenance_path;
    size_t n_threads = 1;
    bool verbose = false;

    CLI::App app{"GetBlunted v1.0.0"};

    app.add_option(
            "-i,--input_gfa",
            gfa_path,
            "Path to GFA containing overlaps")
            ->required();

    app.add_option(
            "-p,--provenance",
            provenance_path,
            "Optionally generate a table containing info about the origin of each output node");

    app.add_option(
            "-t,--threads",
            n_threads,
            "Number of threads to use (maximum)");

    app.add_flag(
            "-V,--verbose",
            verbose,
            "Print a timed log showing progress/steps");

    CLI11_PARSE(app, argc, argv);

    // test input for openability
    if (!ifstream(gfa_path)) {
        cerr << "ERROR: could not open input GFA " << gfa_path << endl;
        return 1;
    }
    
    if (verbose) {
        cerr << "[get_blunted] Executing command:";
        for (int i = 0; i < argc; ++i) {
            cerr << ' ' << argv[i];
        }
        cerr << endl;
    }
    
    Bluntifier bluntifier(gfa_path, provenance_path, verbose);
    bluntifier.bluntify(n_threads);

    return 0;
}

