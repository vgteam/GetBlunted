#include "Bluntifier.hpp"

#include <iostream>
#include <getopt.h>

using bluntifier::Bluntifier;
using std::ifstream;
using std::cerr;
using std::cout;
using std::endl;

void print_usage() {
    cerr << "usage: get_blunted [options] overlap_graph.gfa > blunt_graph.gfa" << endl;
    cerr << endl;
    cerr << "options:" << endl;
    cerr << " -p, --provenance FILEPATH track origin of bluntified sequences in a table here" << endl;
    cerr << " -V, --verbose             log progress to stderr during execution" << endl;
    cerr << " -v, --version             print the version to stdout and exit" << endl;
    cerr << " -h, --help                print this help message to stderr and exit" << endl;
}

void print_version() {
    cout << "get_blunted version 0.0.1" << endl;
}

int main(int argc, char **argv){
    
    string provenance_path;
    bool verbose = false;
    
    int c;
    while (true){
        static struct option long_options[] =
        {
            {"provenance", required_argument, 0, 'p'},
            {"verbose", no_argument, 0, 'V'},
            {"version", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            {0,0,0,0}
        };
        
        int option_index = 0;
        c = getopt_long (argc, argv, "p:Vvh",
                         long_options, &option_index);
        if (c == -1){
            break;
        }
        
        switch(c){
            case 'p':
                provenance_path = optarg;
                break;
            case 'V':
                verbose = true;
                break;
            case 'v':
                print_version();
                return 0;
            case 'h':
            case '?':
                print_usage();
                return 0;
            default:
                print_usage();
                return 1;
        }
    }
    
    if (optind >= argc) {
        // no positional arguments
        cerr << "ERROR: overlap GFA file is required" << endl;
        print_usage();
        return 1;
    }
    
    string gfa_path = argv[optind];
    ++optind;
    
    if (optind < argc) {
        // arguments are included past the positional argument
        
        cerr << "ERROR: unused argument(s):" << endl;
        while (optind < argc) {
            cerr << "\t" << argv[optind] << endl;
            ++optind;
        }
        return 1;
    }
    
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
    bluntifier.bluntify();

    return 0;
}

